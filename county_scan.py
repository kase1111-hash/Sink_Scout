#!/usr/bin/env python3
"""
county_scan.py - run sink_scout over a whole folder of USGS 1-m tiles and merge the results into one
county-wide candidate list and map.

Each tile (or each piece of a tile, with --split) is run together with a buffer pulled from its neighbors,
so pits on tile seams are measured whole and catchments do not stop at the tile edge. Only candidates whose
center falls inside the tile (or piece) itself are kept, so nothing is counted twice.

    python county_scan.py TILE_FOLDER --out county --dry-run          # plan, cell counts, time estimate
    python county_scan.py TILE_FOLDER --out county --osm-roads        # 2 m per tile (default), ~2-3 h for 98 tiles
    python county_scan.py TILE_FOLDER --out county --agg 1 --split 3 --osm-roads   # true 1 m, ~9 pieces per tile
    python county_scan.py TILE_FOLDER --out county --merge-only       # rebuild the county files from finished runs

Resume-safe: stop it any time (Ctrl+C) and run the same command again; finished pieces are skipped.

Outputs in --out:
    runs/<tile>[_rRcC]/           one sink_scout run per tile or piece (log in run.log)
    roads/<tile>.geojson          OSM roads, fetched once per tile (--osm-roads)
    county_candidates.geojson     kept pits: points + footprints, EPSG:4326, with tile, layer, county_rank
    county_candidates.csv         the same table, no shapes
    county_candidates.kml         Google Earth: folders for pits, off-channel pits, big-wash pits
    county_chains.kml / .geojson  pit chains whose top pit lies in the kept area
    county_overview.png           hillshade of the whole county with every kept pit, top 100 labeled
    county_density_nocoords.png   pit density in 2-km cells, no coordinates or labels (publishable)
    county_summary.txt            per-tile counts, failures, run times

Layers (no per-tile limit by default; set one with the --top options if a county is too crowded):
    pits        on a channel, no shallow/wide, creek or road flag            --top
    offchannel  pit-shaped but off any channel: shafts, cellars, arches,      --top-offchannel
                and, in claystone (Wasatch) hillsides, many real cave entrances
    bigwash     pit-shaped on a channel draining more than --max-catch-ha     --top-bigwash
    shallow     shallow/wide hollows, mostly ponds, a few caves (only with --keep-shallow)  --top-shallow
    Road-flagged depressions are dropped; shallow/wide ones are dropped unless --keep-shallow.
    (Checked Oct 2026 against one surveyor's list of 229 Mesa County caves: per-tile limits of 50/10/10
    cut caves with a candidate within 100 m of an entrance from 124 to 91; keeping shallow hollows added
    2,437 candidates for 8 more caves within 25 m.)

KEEP THE OUTPUT PRIVATE. It is a map of cave locations. The density figure is the only one meant for print.
"""

import argparse
import csv
import glob
import json
import math
import os
import re
import subprocess
import sys
import time

LAYERS = ("pits", "offchannel", "bigwash", "shallow")
TILE_RE = re.compile(r"x(\d+)y(\d+)_(.+?)\.tif$", re.IGNORECASE)
BYTES_PER_CELL = 180          # measured peak memory of sink_scout per DEM cell
SEC_PER_MCELL = 3.8           # measured on Kase's desktop: 25 M cells in ~95 s


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------------------- plan

def find_tiles(folder):
    import rasterio
    tiles = []
    for p in sorted(glob.glob(os.path.join(folder, "*.tif"))):
        m = TILE_RE.search(os.path.basename(p))
        if not m:
            continue
        with rasterio.open(p) as src:
            b = src.bounds
            crs = src.crs
            res = abs(src.transform.a)
        col, row = int(m.group(1)), int(m.group(2))
        bounds = (b.left, b.bottom, b.right, b.top)
        # USGS 1-m tiles overlap their neighbors by a few meters. Split the ground on the nominal
        # 10-km grid the name describes, so a pit on a seam belongs to exactly one tile.
        nominal = (col * 10000.0, (row - 1) * 10000.0, (col + 1) * 10000.0, row * 10000.0)
        if all(abs(u - v) < 100 for u, v in zip(bounds, nominal)):
            bounds = nominal
        tiles.append(dict(path=p, col=col, row=row, project=m.group(3), res=res,
                          id=f"x{m.group(1)}y{m.group(2)}", bounds=bounds,
                          crs=crs.to_string() if crs else None, crs_obj=crs))
    # two projects can share a tile number (UTM 12 and 13): make ids unique
    seen = {}
    for t in tiles:
        seen.setdefault(t["id"], []).append(t)
    for tid, group in seen.items():
        if len(group) > 1:
            for t in group:
                zone = re.search(r"UTM(\d+)", t["project"])
                t["id"] = f"{tid}_z{zone.group(1)}" if zone else f"{tid}_{t['project'][:12]}"
    return tiles


def overlaps(a, b):
    return a[0] < b[2] and a[2] > b[0] and a[1] < b[3] and a[3] > b[1]


def plan_units(tiles, split, buffer_m):
    units = []
    for t in tiles:
        l, b, r, tp = t["bounds"]
        w, h = (r - l) / split, (tp - b) / split
        same = [o for o in tiles if o["project"] == t["project"] and o["crs"] == t["crs"]]
        for i in range(split):            # i = row from the top
            for j in range(split):
                core = (l + j * w, tp - (i + 1) * h, l + (j + 1) * w, tp - i * h)
                buf = (core[0] - buffer_m, core[1] - buffer_m, core[2] + buffer_m, core[3] + buffer_m)
                srcs = [o for o in same if overlaps(o["bounds"], buf)]
                srcs.sort(key=lambda o: o is not t)          # the tile itself first (sets the grid)
                name = t["id"] if split == 1 else f"{t['id']}_r{i}c{j}"
                units.append(dict(name=name, tile=t, core=core, buf=buf, srcs=srcs))
    return units


def buf_to_bbox_arg(unit):
    """lat1,lon1,lat2,lon2 for sink_scout --bbox. sink_scout snaps it back to the grid in the DEM CRS."""
    from rasterio.warp import transform_bounds
    from rasterio.crs import CRS
    w, s, e, n = transform_bounds(unit["tile"]["crs_obj"], CRS.from_epsg(4326), *unit["buf"], densify_pts=21)
    return f"{s:.6f},{w:.6f},{n:.6f},{e:.6f}"


def estimate(units, agg):
    # the lat/lon box sink_scout gets back is a little larger than the buffered rectangle: ~4% each way
    cells = [((u["buf"][2] - u["buf"][0]) * 1.04 / (agg * u["tile"]["res"]))
             * ((u["buf"][3] - u["buf"][1]) * 1.04 / (agg * u["tile"]["res"])) for u in units]
    return max(cells), sum(cells)


# ----------------------------------------------------------------------------- run

def fetch_roads(tile, buffer_m, roads_dir, sink_scout_dir):
    """One Overpass request per tile, covering the tile plus its buffer; every piece of the tile reuses it."""
    path = os.path.join(roads_dir, f"{tile['id']}.geojson")
    if os.path.exists(path):
        return path
    sys.path.insert(0, sink_scout_dir)
    import sink_scout
    from rasterio.transform import from_bounds
    l, b, r, t = tile["bounds"]
    bb = (l - buffer_m, b - buffer_m, r + buffer_m, t + buffer_m)
    sink_scout.fetch_osm_roads(tile["crs_obj"], from_bounds(*bb, 100, 100), (100, 100), path)
    time.sleep(2)     # be polite to the Overpass servers
    return path if os.path.exists(path) else None


def previous_run(out):
    """(set of source tile names, agg) a finished run used, read from its DONE file or, for runs made by
    the first version of this script, from the command line at the top of run.log."""
    done = os.path.join(out, "DONE")
    try:
        info = json.load(open(done))
        if isinstance(info, dict):
            return set(info["srcs"]), int(info["agg"])
    except Exception:
        pass
    try:
        first = open(os.path.join(out, "run.log")).readline()
    except OSError:
        return None, None
    import shlex
    try:
        toks = [t.strip('"') for t in shlex.split(first, posix=False)]
    except ValueError:
        toks = [t.strip('"') for t in first.split()]
    srcs = {os.path.basename(t.replace("\\", "/")) for t in toks if t.lower().endswith(".tif")}
    agg = int(toks[toks.index("--agg") + 1]) if "--agg" in toks else 1
    return srcs, agg


def needs_run(unit, args):
    """A run is redone if it never finished, or if the tiles around it have changed since (a tile added
    next to it) or the cell size has, because its buffer then held different ground."""
    out = os.path.join(args.out, "runs", unit["name"])
    if not os.path.exists(os.path.join(out, "DONE")):
        return True
    srcs, agg = previous_run(out)
    want = {os.path.basename(s["path"]) for s in unit["srcs"]}
    return srcs != want or agg != args.agg


def run_unit(unit, args, sink_scout_path, roads_path):
    out = os.path.join(args.out, "runs", unit["name"])
    done = os.path.join(out, "DONE")
    if not needs_run(unit, args):
        return "skip", 0.0
    if os.path.exists(done):
        os.remove(done)
    os.makedirs(out, exist_ok=True)
    cmd = [args.python or sys.executable, sink_scout_path] + [s["path"] for s in unit["srcs"]]
    if len(unit["srcs"]) > 1:
        cmd.append("--mosaic")
    cmd += ["--bbox", buf_to_bbox_arg(unit), "--agg", str(args.agg), "--min-ratio", f"{args.min_ratio:g}",
            "--max-catch-ha", f"{args.max_catch_ha:g}", "--thumbs", "0", "--no-overview", "--out", out]
    if roads_path:
        cmd += ["--roads", roads_path]
    if args.extra:
        cmd += args.extra.split()
    t0 = time.time()
    with open(os.path.join(out, "run.log"), "w") as lf:
        lf.write(" ".join(f'"{c}"' if " " in c else c for c in cmd) + "\n\n")
        lf.flush()
        rc = subprocess.call(cmd, stdout=lf, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    if rc == 0:
        with open(done, "w") as f:
            json.dump({"secs": round(dt, 1), "agg": args.agg,
                       "srcs": sorted(os.path.basename(s["path"]) for s in unit["srcs"])}, f)
        return "ok", dt
    return f"failed (exit {rc}, see {os.path.join(out, 'run.log')})", dt


# ----------------------------------------------------------------------------- merge

def layer_of(p, max_catch_ac, keep_shallow=False):
    flags = p.get("flags") or ""
    if "road" in flags:
        return None
    if "shallow/wide" in flags:
        return "shallow" if keep_shallow else None
    if "off-channel" in flags:
        return "offchannel"
    if "creek" in flags:
        return "bigwash"
    return "pits"


def merge(units, args):
    from rasterio.warp import transform as rtransform
    from rasterio.crs import CRS
    wgs = CRS.from_epsg(4326)
    caps = {"pits": args.top, "offchannel": args.top_offchannel, "bigwash": args.top_bigwash,
            "shallow": args.top_shallow}
    by_tile, chains, polys, per_tile = {}, [], {}, {}
    missing = 0
    for u in units:
        d = os.path.join(args.out, "runs", u["name"])
        tid = u["tile"]["id"]
        st = per_tile.setdefault(tid, {"pieces": 0, "done": 0, "raw": 0, "in_core": 0})
        st["pieces"] += 1
        if not os.path.exists(os.path.join(d, "DONE")):
            missing += 1
            continue
        st["done"] += 1
        cg = os.path.join(d, "candidates.geojson")
        if not os.path.exists(cg):
            continue            # ran fine, nothing passed the filters
        gj = json.load(open(cg))
        l, b, r, t = u["core"]
        local_polys = {}
        for f in gj["features"]:
            p = f["properties"]
            if f["geometry"]["type"] != "Point":
                local_polys[p["rank"]] = f["geometry"]
                continue
            st["raw"] += 1
            if not (l <= p["x"] < r and b <= p["y"] < t):
                continue
            st["in_core"] += 1
            lay = layer_of(p, args.max_catch_ha * 2.471, args.keep_shallow)
            if lay is None:
                continue
            q = dict(p)
            q.update(tile=tid, piece=u["name"], layer=lay, run_rank=p["rank"])
            by_tile.setdefault((tid, lay), []).append(q)
            polys[(u["name"], p["rank"])] = None
        for key in list(polys):
            if key[0] == u["name"] and key[1] in local_polys:
                polys[key] = local_polys[key[1]]
        chg = os.path.join(d, "chains.geojson")
        if os.path.exists(chg):
            for f in json.load(open(chg))["features"]:
                lon0, lat0 = f["geometry"]["coordinates"][0]
                xs, ys = rtransform(wgs, u["tile"]["crs_obj"], [lon0], [lat0])
                if l <= xs[0] < r and b <= ys[0] < t:
                    f["properties"].update(tile=tid, piece=u["name"])
                    chains.append(f)
    # The same pit can come back from two runs: where two USGS projects (UTM zones) overlap, and on a
    # seam if a tile's edges are not on the nominal grid. Measured from different grids it lands within
    # a few meters of itself (all such pairs in the Mesa County run were under 5 m), while separate pits
    # found by one run were never closer than 8.6 m. So: one pit per 6 m across different runs, never
    # merging two candidates from the same run, keeping the higher score.
    allq = sorted((q for lst in by_tile.values() for q in lst), key=lambda q: -q["score"])
    grid, keep_ids, n_dup = {}, set(), 0
    gsz = 0.0002
    for q in allq:
        gi, gj = int(q["lat"] / gsz), int(q["lon"] / gsz)
        dup = False
        for a in (gi - 1, gi, gi + 1):
            for b in (gj - 1, gj, gj + 1):
                for o in grid.get((a, b), ()):
                    if o["piece"] != q["piece"] and math.hypot(
                            (o["lat"] - q["lat"]) * 110574,
                            (o["lon"] - q["lon"]) * 111320 * math.cos(math.radians(q["lat"]))) < 6.0:
                        dup = True
        if dup:
            n_dup += 1
            continue
        grid.setdefault((gi, gj), []).append(q)
        keep_ids.add(id(q))
    for key in list(by_tile):
        by_tile[key] = [q for q in by_tile[key] if id(q) in keep_ids]
    uniq_chains, cgrid = [], {}
    for c in chains:
        lo0, la0 = c["geometry"]["coordinates"][0]
        gi, gj = int(la0 / gsz), int(lo0 / gsz)
        near = [o for a in (gi - 1, gi, gi + 1) for b in (gj - 1, gj, gj + 1) for o in cgrid.get((a, b), ())]
        if any(o["properties"]["piece"] != c["properties"]["piece"] and math.hypot(
                (o["geometry"]["coordinates"][0][1] - la0) * 110574,
                (o["geometry"]["coordinates"][0][0] - lo0) * 111320 * math.cos(math.radians(la0))) < 6.0
               for o in near):
            continue
        cgrid.setdefault((gi, gj), []).append(c)
        uniq_chains.append(c)
    chains = uniq_chains
    if n_dup:
        log(f"{n_dup} duplicate candidate(s) removed (the same pit found by two runs)")
    args._n_dup = n_dup
    kept = []
    for (tid, lay), lst in by_tile.items():
        lst.sort(key=lambda q: -q["score"])
        cap = caps[lay] if caps[lay] and caps[lay] > 0 else None       # 0 = no limit
        kept += lst[:cap]
        per_tile[tid][lay] = len(lst[:cap])
        per_tile[tid][lay + "_avail"] = len(lst)
    for lay in caps:
        sub = sorted((q for q in kept if q["layer"] == lay), key=lambda q: -q["score"])
        for i, q in enumerate(sub, 1):
            q["county_rank"] = i
    kept.sort(key=lambda q: (LAYERS.index(q["layer"]), q["county_rank"]))
    write_county(kept, polys, chains, per_tile, missing, args)
    return kept, chains, per_tile, missing


def write_county(kept, polys, chains, per_tile, missing, args):
    o = args.out
    feats = []
    for q in kept:
        props = {k: v for k, v in q.items() if k not in ("label", "parent", "x", "y")}
        feats.append({"type": "Feature", "properties": props,
                      "geometry": {"type": "Point", "coordinates": [q["lon"], q["lat"]]}})
        g = polys.get((q["piece"], q["run_rank"]))
        if g:
            feats.append({"type": "Feature", "properties": {"layer": q["layer"], "county_rank": q["county_rank"],
                                                            "tile": q["tile"]}, "geometry": g})
    json.dump({"type": "FeatureCollection", "features": feats}, open(os.path.join(o, "county_candidates.geojson"), "w"))
    cols = ["layer", "county_rank", "tile", "piece", "run_rank", "score", "depth_ft", "diam_ft", "ratio", "elong",
            "catch_ac", "road_m", "chain_n", "chain_len_m", "min_z_m", "flags", "lat", "lon"]
    with open(os.path.join(o, "county_candidates.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for q in kept:
            w.writerow([round(q[c], 5) if isinstance(q.get(c), float) else q.get(c, "") for c in cols])
    json.dump({"type": "FeatureCollection", "features": chains}, open(os.path.join(o, "county_chains.geojson"), "w"))

    # KML: one folder per layer, colored by score (yellow = best)
    def color(score, smax):
        t = max(0.0, min(1.0, score / smax)) if smax > 0 else 0
        r_, g_, b_ = 255, int(60 + 175 * t), int(40 * (1 - t))
        return f"ff{b_:02x}{g_:02x}{r_:02x}"
    titles = {"pits": "Pits on washes", "offchannel": "Off-channel pits (shafts, cellars, arches, alcoves)",
              "bigwash": "Pits on big washes (creek flag)",
              "shallow": "Shallow, wide hollows (mostly ponds; a few caves)"}
    k = ['<?xml version="1.0" encoding="UTF-8"?>', '<kml xmlns="http://www.opengis.net/kml/2.2"><Document>',
         '<name>sink_scout county scan</name>']
    for lay in LAYERS:
        sub = [q for q in kept if q["layer"] == lay]
        if not sub:
            continue
        smax = max(q["score"] for q in sub)
        k.append(f"<Folder><name>{titles[lay]} ({len(sub)})</name>")
        for q in sub:
            desc = (f"score {q['score']:.2f}<br/>depth {q['depth_ft']:.0f} ft, diameter {q['diam_ft']:.0f} ft, "
                    f"d/w {q['ratio']:.2f}<br/>catchment {q['catch_ac']:.0f} ac<br/>chain of {q['chain_n']}"
                    f"<br/>floor {q['min_z_m'] * 3.2808:,.0f} ft<br/>flags {q['flags'] or 'none'}"
                    f"<br/>tile {q['tile']}, run {q['piece']} rank {q['run_rank']}")
            k.append(f"<Placemark><name>{lay[0].upper()}{q['county_rank']}</name><description><![CDATA[{desc}]]>"
                     f"</description><Style><IconStyle><color>{color(q['score'], smax)}</color><scale>0.8</scale>"
                     f"<Icon><href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href></Icon>"
                     f"</IconStyle></Style><Point><coordinates>{q['lon']:.7f},{q['lat']:.7f},0</coordinates>"
                     f"</Point></Placemark>")
        k.append("</Folder>")
    if chains:
        k.append(f"<Folder><name>Chains ({len(chains)})</name>")
        for c in chains:
            p = c["properties"]
            k.append(f"<Placemark><name>{p['n_pits']} pits, {p['stem_ft']:,} ft</name><Style><LineStyle>"
                     f"<color>ff00d7ff</color><width>3</width></LineStyle></Style><LineString><tessellate>1"
                     f"</tessellate><coordinates>" + " ".join(f"{a:.7f},{b:.7f},0" for a, b in c["geometry"]["coordinates"])
                     + "</coordinates></LineString></Placemark>")
        k.append("</Folder>")
    k.append("</Document></kml>")
    open(os.path.join(o, "county_candidates.kml"), "w").write("\n".join(k))

    with open(os.path.join(o, "county_summary.txt"), "w") as f:
        n = {lay: sum(1 for q in kept if q["layer"] == lay) for lay in LAYERS}
        f.write(f"county scan: {len(per_tile)} tiles, {sum(v['pieces'] for v in per_tile.values())} runs "
                f"({missing} not finished)\nkept: {n['pits']} pits, {n['offchannel']} off-channel, "
                f"{n['bigwash']} big-wash, " + (f"{n['shallow']} shallow, " if n['shallow'] else "") + f"{len(chains)} chains; "
                f"{getattr(args, '_n_dup', 0)} duplicates removed\n\n")
        f.write(f"{'tile':<16}{'runs':>6}{'raw':>7}{'inside':>8}{'pits':>12}{'offch':>10}{'bigwash':>10}\n")
        for tid in sorted(per_tile):
            v = per_tile[tid]
            f.write(f"{tid:<16}{v['done']:>3}/{v['pieces']:<2}{v['raw']:>7}{v['in_core']:>8}"
                    f"{v.get('pits', 0):>5}/{v.get('pits_avail', 0):<6}{v.get('offchannel', 0):>4}/"
                    f"{v.get('offchannel_avail', 0):<5}{v.get('bigwash', 0):>4}/{v.get('bigwash_avail', 0):<5}\n")
    log(f"county files written: {n['pits']} pits, {n['offchannel']} off-channel, {n['bigwash']} big-wash, "
        f"{len(chains)} chains")


# ----------------------------------------------------------------------------- maps

def county_hillshade(tiles, res):
    """Decimated mosaic of every tile in one CRS (the most common one), for the overview maps."""
    import numpy as np
    import rasterio
    from rasterio.warp import reproject, transform_bounds, Resampling
    from rasterio.transform import from_origin
    from collections import Counter
    crs = Counter(t["crs"] for t in tiles).most_common(1)[0][0]
    crs_obj = next(t["crs_obj"] for t in tiles if t["crs"] == crs)
    bxs = [transform_bounds(t["crs_obj"], crs_obj, *t["bounds"]) for t in tiles]
    L, B = min(b[0] for b in bxs), min(b[1] for b in bxs)
    R, T = max(b[2] for b in bxs), max(b[3] for b in bxs)
    W, H = int(math.ceil((R - L) / res)), int(math.ceil((T - B) / res))
    dst_t = from_origin(L, T, res, res)
    z = np.full((H, W), np.nan, dtype=np.float32)
    for t in tiles:
        with rasterio.open(t["path"]) as src:
            f = max(1, int(res / abs(src.transform.a)))
            h, w = max(1, src.height // f), max(1, src.width // f)
            a = src.read(1, out_shape=(h, w), resampling=Resampling.average, masked=True).astype(np.float32)
            a = a.filled(np.nan)
            a[(a < -1000) | (a > 9000)] = np.nan
            st = src.transform * src.transform.scale(src.width / w, src.height / h)
            tmp = np.full((H, W), np.nan, dtype=np.float32)
            reproject(a, tmp, src_transform=st, src_crs=src.crs, dst_transform=dst_t, dst_crs=crs_obj,
                      src_nodata=np.nan, dst_nodata=np.nan, resampling=Resampling.bilinear)
            ok = np.isfinite(tmp)
            z[ok] = tmp[ok]
    gy, gx = np.gradient(np.nan_to_num(z, nan=float(np.nanmean(z))), res)
    az, alt = math.radians(315), math.radians(45)
    slope = np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
    hs[~np.isfinite(z)] = np.nan
    return hs, (L, R, B, T), crs_obj


def write_maps(kept, chains, tiles, args):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from rasterio.warp import transform as rtransform
    from rasterio.crs import CRS
    log("building county hillshade for the maps...")
    hs, ext, crs_obj = county_hillshade(tiles, args.map_res)
    wgs = CRS.from_epsg(4326)
    L, R, B, T = ext
    aspect = (T - B) / (R - L)
    pits = [q for q in kept if q["layer"] == "pits"]

    def xy(lst):
        if not lst:
            return np.array([]), np.array([])
        xs, ys = rtransform(wgs, crs_obj, [q["lon"] for q in lst], [q["lat"] for q in lst])
        return np.array(xs), np.array(ys)

    # 1. private overview: every kept pit, top 100 labeled
    fig, ax = plt.subplots(figsize=(16, 16 * aspect + 1))
    ax.imshow(hs, extent=ext, cmap="gray", vmin=-0.2, vmax=1.0, interpolation="bilinear")
    for c in chains:
        xs, ys = rtransform(wgs, crs_obj, *zip(*c["geometry"]["coordinates"]))
        ax.plot(xs, ys, color="orange", lw=1.2)
    for lay, mk, col in (("offchannel", "s", "tab:cyan"), ("bigwash", "^", "tab:blue")):
        x, y = xy([q for q in kept if q["layer"] == lay])
        ax.scatter(x, y, s=10, marker=mk, c=col, edgecolors="none", alpha=0.8, label=lay)
    x, y = xy(pits)
    if len(x):
        sc = ax.scatter(x, y, s=14, c=[q["score"] for q in pits], cmap="plasma", edgecolors="white",
                        linewidths=0.3, label="pits on washes")
        plt.colorbar(sc, ax=ax, fraction=0.025, pad=0.01, label="score")
        for q, xi, yi in list(zip(pits, x, y))[:100]:
            ax.annotate(str(q["county_rank"]), (xi, yi), fontsize=6, xytext=(2, 2), textcoords="offset points",
                        color="black", weight="bold")
    ax.set_xlim(L, R); ax.set_ylim(B, T)
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title(f"sink_scout county scan: {len(pits)} pits on washes, "
                 f"{sum(q['layer'] == 'offchannel' for q in kept)} off-channel, "
                 f"{sum(q['layer'] == 'bigwash' for q in kept)} on big washes, {len(chains)} chains")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "county_overview.png"), dpi=args.map_dpi)
    plt.close(fig)

    # 2. publishable density: counts of pits per cell, no axes, no labels, no points
    cell = args.density_km * 1000
    nx, ny = int(math.ceil((R - L) / cell)), int(math.ceil((T - B) / cell))
    grid = np.zeros((ny, nx))
    dx, dy = xy([q for q in kept if q["layer"] in ("pits", "bigwash")])
    for xi, yi in zip(dx, dy):
        grid[min(ny - 1, int((T - yi) // cell)), min(nx - 1, int((xi - L) // cell))] += 1
    fig, ax = plt.subplots(figsize=(12, 12 * aspect + 0.6))
    ax.imshow(hs, extent=ext, cmap="gray", vmin=-0.2, vmax=1.0, interpolation="bilinear")
    g = np.ma.masked_where(grid < 1, grid)
    im = ax.imshow(g, extent=(L, L + nx * cell, T - ny * cell, T), cmap="inferno_r", alpha=0.65,
                   interpolation="nearest", vmin=1)
    plt.colorbar(im, ax=ax, fraction=0.025, pad=0.01, label=f"candidate pits per {args.density_km:g}-km square")
    sb = 10000
    ax.plot([L + 0.04 * (R - L), L + 0.04 * (R - L) + sb], [B + 0.04 * (T - B)] * 2, color="black", lw=3)
    ax.text(L + 0.04 * (R - L), B + 0.06 * (T - B), "10 km   north up", fontsize=9)
    ax.set_xlim(L, R); ax.set_ylim(B, T); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "county_density_nocoords.png"), dpi=args.map_dpi)
    plt.close(fig)
    log("maps written: county_overview.png (private), county_density_nocoords.png (publishable)")


# ----------------------------------------------------------------------------- main

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("folder", help="folder holding the USGS_one_meter_x..y.._*.tif tiles")
    p.add_argument("--out", default="county")
    p.add_argument("--sink-scout", help="path to sink_scout.py (default: next to this script)")
    p.add_argument("--python", help="python executable for the runs (default: this one)")
    p.add_argument("--agg", type=int, default=2, help="cell size in meters for the runs (default 2)")
    p.add_argument("--split", type=int, default=1,
                   help="split each tile into N x N pieces (use 3 with --agg 1 so each run fits in memory)")
    p.add_argument("--buffer", type=float, default=500.0, help="meters borrowed from neighbors around each piece")
    p.add_argument("--min-ratio", type=float, default=0.2)
    p.add_argument("--max-catch-ha", type=float, default=300.0)
    p.add_argument("--osm-roads", action="store_true", help="fetch OSM roads once per tile and use them in every run")
    p.add_argument("--extra", help='extra sink_scout options, quoted, e.g. --extra "--min-depth 2"')
    p.add_argument("--top", type=int, default=0, help="pits on washes kept per tile (default 0 = all)")
    p.add_argument("--top-offchannel", type=int, default=0, help="off-channel pits kept per tile (default 0 = all)")
    p.add_argument("--top-bigwash", type=int, default=0, help="big-wash pits kept per tile (default 0 = all)")
    p.add_argument("--keep-shallow", action="store_true",
                   help="keep shallow/wide hollows as their own layer instead of dropping them (mostly ponds)")
    p.add_argument("--top-shallow", type=int, default=0, help="shallow hollows kept per tile (default 0 = all)")
    p.add_argument("--only", help="comma-separated tile ids to run (e.g. x70y433,x71y433) - for a test")
    p.add_argument("--dry-run", action="store_true", help="print the plan and estimates, run nothing")
    p.add_argument("--merge-only", action="store_true", help="skip the runs, rebuild the county files and maps")
    p.add_argument("--no-maps", action="store_true")
    p.add_argument("--map-res", type=float, default=30.0, help="hillshade cell for the county maps (m, default 30)")
    p.add_argument("--map-dpi", type=int, default=200)
    p.add_argument("--density-km", type=float, default=2.0, help="cell size of the publishable density map")
    args = p.parse_args()

    sink_scout_path = os.path.abspath(args.sink_scout or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                                       "sink_scout.py"))
    if not os.path.exists(sink_scout_path):
        sys.exit(f"sink_scout.py not found at {sink_scout_path}; pass --sink-scout")
    tiles = find_tiles(args.folder)
    if not tiles:
        sys.exit(f"no USGS_one_meter_x..y.._*.tif tiles in {args.folder}")
    projects = sorted({t["project"] for t in tiles})
    units = plan_units(tiles, args.split, args.buffer)
    if args.only:
        want = set(args.only.split(","))
        units = [u for u in units if u["tile"]["id"] in want]
    biggest, total = estimate(units, args.agg)
    log(f"{len(tiles)} tiles in {len(projects)} project(s): " + ", ".join(
        f"{pr} ({sum(t['project'] == pr for t in tiles)})" for pr in projects))
    log(f"{len(units)} runs at {args.agg * tiles[0]['res']:g} m cells, {args.buffer:g} m buffer; largest run {biggest / 1e6:.0f} M cells "
        f"(~{biggest * BYTES_PER_CELL / 1e9:.1f} GB peak); total {total / 1e9:.2f} G cells, "
        f"~{total / 1e6 * SEC_PER_MCELL / 3600:.1f} h at the measured speed")
    if biggest > 30e6:
        log("WARNING: runs this large have not been tried on your machine; consider --split "
            f"{args.split + 1} or --agg {args.agg + 1}")
    if len(projects) > 1:
        log("note: tiles from different projects (UTM zones) are not mosaicked together; "
            "pits on the zone line are measured from one side only")
    todo = [u for u in units if needs_run(u, args)]
    redo = [u for u in todo if os.path.exists(os.path.join(args.out, "runs", u["name"], "DONE"))]
    if len(todo) < len(units):
        log(f"{len(units) - len(todo)} run(s) already finished; {len(todo)} to do"
            + (f", of which {len(redo)} are reruns because the tiles around them changed" if redo else ""))
    if args.dry_run:
        for u in redo[:10]:
            print(f"  rerun {u['name']}: sources now {', '.join(s['id'] for s in u['srcs'])}")
        for u in units[:5]:
            print(f"  {u['name']}: {len(u['srcs'])} source tile(s), --bbox {buf_to_bbox_arg(u)}")
        if len(units) > 5:
            print(f"  ... {len(units) - 5} more")
        return

    os.makedirs(os.path.join(args.out, "runs"), exist_ok=True)
    if not args.merge_only:
        roads_dir = os.path.join(args.out, "roads")
        os.makedirs(roads_dir, exist_ok=True)
        t_start, done_secs, n_done = time.time(), 0.0, 0
        failures = []
        for k, u in enumerate(units, 1):
            roads = None
            if args.osm_roads:
                try:
                    roads = fetch_roads(u["tile"], args.buffer, roads_dir, os.path.dirname(sink_scout_path))
                except Exception as ex:
                    log(f"  roads for {u['tile']['id']} failed ({ex}); running without")
            status, dt = run_unit(u, args, sink_scout_path, roads)
            if status == "skip":
                continue
            n_done += 1
            done_secs += dt
            left = sum(1 for v in units[k:] if needs_run(v, args))
            eta = done_secs / n_done * left / 3600
            log(f"[{k}/{len(units)}] {u['name']}: {status} in {dt:.0f} s; {left} left, ~{eta:.1f} h")
            if status != "ok":
                failures.append(u["name"])
        if failures:
            log(f"{len(failures)} run(s) failed: {', '.join(failures[:20])} - fix and rerun the same command")
    kept, chains, per_tile, missing = merge(units, args)
    if missing:
        log(f"{missing} run(s) not finished yet; the county files cover the rest")
    if not args.no_maps and kept:
        write_maps(kept, chains, [t for t in tiles if any(u["tile"] is t for u in units)], args)


if __name__ == "__main__":
    main()
