#!/usr/bin/env python3
"""
field_kit.py - turn a sink_scout or county_scan candidate list into what you take into the field.

    python field_kit.py county/county_candidates.geojson --tiles tiles --chains county/county_chains.geojson --out field
    python field_kit.py county/county_candidates.geojson --tiles tiles --chains county/county_chains.geojson \
        --truth ground_truth.csv --tag land=sma.shp:ADMIN_AGEN --tag geology=geology.shp:UNIT_NAME \
        --sheets top25 --out field

It reads the candidate list and the same elevation tiles the scan used. Nothing here reruns the scan, so it
works on any finished run and can be run again whenever a verdict or a map layer is added.

Outputs in --out:
    field_candidates.kmz       Google Earth: every candidate with its lidar close-up in the balloon, colored by
                               verdict where one is known; chains with their profiles
    field_candidates.gpx       waypoints for a GPS or phone, and each chain as a track
    field_candidates.csv       the list with close-up names, tags and verdicts
    closeups/<ID>.png          hillshade with 1-m contours, the pit's footprint and its neighbours, 300 m across;
                               it shows the ground under snow, trees and shadow, where photographs cannot
    profiles/<chain>.png       the ground down each chain, the pit floors, the straight line joining them, and
                               how much ground lies above that line between pits (relative heights only)
    profiles/chains.csv        per chain: pits, wash length, floor drop, least and median ground above the line
    sheets.pdf                 with --sheets: one page per pit for the clipboard: close-up, numbers, profile,
                               and the checks from the safety chapter with boxes to tick
    truth_report.txt           with --truth: how often each kind of candidate turned out to be a pipe, by layer
                               and by rank, and which independently known caves the scan found or missed
    ground_truth_worklist.csv  every candidate, with the verdicts already given filled in; add verdicts as you
                               check them and use it as the next --truth file. Verdicts follow the pit, not its
                               rank, so they survive reruns.

Verdict words (first word counts): pipe | known | talus | mine | cellar | arch | slot | plunge | pothole |
pond | culvert | road | logs | artifact | unknown.  Use "known" only for caves known before or apart from the
program (on the list, surveyed by someone, found on foot): those are the ones that measure what it misses.

KEEP THE OUTPUTS PRIVATE. Everything here except the chain profiles carries cave locations.
"""

import argparse
import csv
import glob
import json
import math
import os
import re
import statistics
import sys
import time
import zipfile

import numpy as np

M2FT = 3.28084
PIPE = {"pipe", "known", "walked", "surveyed", "cave"}
REAL = {"talus", "mine", "shaft", "cellar", "arch", "alcove", "slot", "plunge", "pothole", "sinkhole"}
FALSE = {"pond", "reservoir", "culvert", "road", "logs", "log", "dam", "artifact", "building", "none", "nothing"}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def vclass(verdict):
    w = re.split(r"[\s,;:/()]+", (verdict or "").strip().lower())
    w = w[0] if w and w[0] else ""
    if w in PIPE:
        return "pipe"
    if w in REAL:
        return "real"
    if w in FALSE:
        return "false"
    return "unknown"


# ----------------------------------------------------------------------------- inputs

def load_candidates(path):
    gj = json.load(open(path))
    cands, polys = [], {}
    for f in gj["features"]:
        p = f["properties"]
        g = f["geometry"]
        if g["type"] == "Point":
            c = dict(p)
            if "layer" in p and p.get("county_rank") is not None:
                c["id"] = f"{p['layer'][0].upper()}{p['county_rank']}"
                c["rank_key"] = p["county_rank"]
            else:
                c["layer"] = "pits" if p.get("on_channel", True) else "offchannel"
                c["id"] = f"#{p['rank']}"
                c["rank_key"] = p["rank"]
            c["lon"], c["lat"] = g["coordinates"][0], g["coordinates"][1]
            cands.append(c)
        elif g["type"] == "Polygon":
            if "layer" in p and p.get("county_rank") is not None:
                key = f"{p['layer'][0].upper()}{p['county_rank']}"
            elif p.get("rank") is not None:
                key = f"#{p['rank']}"
            else:
                continue
            polys[key] = g["coordinates"]
    order = {"pits": 0, "offchannel": 1, "bigwash": 2}
    cands.sort(key=lambda c: (order.get(c["layer"], 3), c["rank_key"]))
    return cands, polys


def load_chains(path, cands):
    """Chains, each with its full member list. A county chain can include pits the county list dropped (shallow
    ones, or ones past a tile's cap); those are read from the run's own candidate file in runs/<piece>/, so a
    profile shows every pit in the chain. Members also in the list keep their list id (P30); the others are
    labeled by their rank in that run (r12)."""
    if not path or not os.path.exists(path):
        return []
    feats = [f for f in json.load(open(path))["features"] if f["geometry"]["type"] == "LineString"]
    feats.sort(key=lambda f: -f["properties"].get("stem_ft", 0))
    runs_dir = os.path.join(os.path.dirname(os.path.abspath(path)), "runs")
    by_piece = {(c.get("piece"), c.get("run_rank")): c for c in cands if c.get("piece")}
    by_rank = {c.get("rank"): c for c in cands if not c.get("piece")}
    cache = {}

    def member(c=None, rank=None, q=None, lonlat=None):
        if c is not None:
            return dict(id=c["id"], lat=c["lat"], lon=c["lon"], floor=c["min_z_m"],
                        rim=c.get("spill_z_m", c["min_z_m"]), r=c.get("diam_m", 10.0) / 2.0, cand=c)
        return dict(id=f"r{rank}", lat=lonlat[1], lon=lonlat[0], floor=q["min_z_m"],
                    rim=q.get("spill_z_m", q["min_z_m"]), r=q.get("diam_m", 10.0) / 2.0, cand=None)

    chains, n_from_runs = [], 0
    for i, f in enumerate(feats, 1):
        p = f["properties"]
        members = []
        piece = p.get("piece")
        for r in p.get("ranks") or []:
            c = by_piece.get((piece, r)) if piece else by_rank.get(r)
            if c is not None:
                members.append(member(c=c))
                continue
            if piece:
                if piece not in cache:
                    rf = os.path.join(runs_dir, piece, "candidates.geojson")
                    cache[piece] = ({g["properties"]["rank"]: g for g in json.load(open(rf))["features"]
                                     if g["geometry"]["type"] == "Point"} if os.path.exists(rf) else {})
                g = cache[piece].get(r)
                if g is not None:
                    members.append(member(rank=r, q=g["properties"], lonlat=g["geometry"]["coordinates"]))
                    n_from_runs += 1
        chains.append(dict(cid=f"C{i}", coords=f["geometry"]["coordinates"], props=p, pits=[], members=members))
    if n_from_runs:
        log(f"{n_from_runs} chain member(s) not in the list were read from the run folders")
    return chains


class Tiles:
    """Every elevation file in a folder, opened once; windows are read across tile seams within one CRS."""

    def __init__(self, folder, crs_override=None, sink_scout=None):
        import rasterio
        self.rasterio = rasterio
        self.ss = sink_scout
        paths = sorted(glob.glob(os.path.join(folder, "*.tif")) + glob.glob(os.path.join(folder, "*.img")))
        if not paths:
            sys.exit(f"no .tif or .img elevation files in {folder}")
        self.items = []
        override = sink_scout.crs_from_arg(crs_override) if crs_override else None
        for p in paths:
            ds = rasterio.open(p)
            crs = override or ds.crs
            if crs is None:
                continue
            self.items.append(dict(path=p, ds=ds, crs=crs, key=crs.to_string(), b=ds.bounds,
                                   f=sink_scout.linear_factor(crs)))
        log(f"{len(self.items)} elevation tiles indexed")
        self._tf = {}

    def _to(self, key, crs, lons, lats):
        from rasterio.warp import transform as rtransform
        return rtransform("EPSG:4326", crs, list(lons), list(lats))

    def group_for(self, lon, lat, half_m):
        """The CRS group with a tile whose data best surrounds the point (furthest from that tile's edge)."""
        best, best_d = None, -1e18
        xy = {}
        for it in self.items:
            if it["key"] not in xy:
                xs, ys = self._to(it["key"], it["crs"], [lon], [lat])
                xy[it["key"]] = (xs[0], ys[0])
            x, y = xy[it["key"]]
            b = it["b"]
            d = min(x - b.left, b.right - x, y - b.bottom, b.top - y) * it["f"]
            if d > best_d:
                best, best_d = it, d
        if best is None or best_d < -half_m:
            return None
        return best

    def window(self, lon, lat, half_m, zfac=1.0):
        """(z in meters, transform, crs, cell_m, center_xy) for a square of +-half_m around the point, or None."""
        from rasterio.merge import merge
        it = self.group_for(lon, lat, half_m)
        if it is None:
            return None
        x, y = self._to(it["key"], it["crs"], [lon], [lat])
        x, y = x[0], y[0]
        h = half_m / it["f"]
        bounds = (x - h, y - h, x + h, y + h)
        srcs = [o["ds"] for o in self.items if o["key"] == it["key"] and
                o["b"].left < bounds[2] and o["b"].right > bounds[0] and o["b"].bottom < bounds[3] and o["b"].top > bounds[1]]
        if not srcs:
            return None
        arr, tr = merge(srcs, bounds=bounds, nodata=float("nan"))
        z = arr[0].astype(np.float64)
        nd = srcs[0].nodata
        if nd is not None and np.isfinite(nd):
            z[z == nd] = np.nan
        z[(z < -1000) | (z > 1e5)] = np.nan
        if not np.isfinite(z).any():
            return None
        return z / zfac, tr, it["crs"], abs(tr.a) * it["f"], (x, y), it["f"]


def detect_zfac(tiles, cands, n=25):
    """1.0 if the tiles are in meters like the candidate list, 3.28 if they are in feet."""
    ratios = []
    for c in cands[:: max(1, len(cands) // n)][:n]:
        w = tiles.window(c["lon"], c["lat"], 5.0)
        if w is None or not c.get("min_z_m"):
            continue
        z = w[0]
        v = np.nanmin(z) if np.isfinite(z).any() else None
        if v:
            ratios.append(v / c["min_z_m"])
    if not ratios:
        return 1.0
    r = statistics.median(ratios)
    return M2FT if abs(r - M2FT) < 0.3 else 1.0


# ----------------------------------------------------------------------------- drawing

def hillshade(z, cell, az=315.0, alt=45.0):
    zz = np.where(np.isfinite(z), z, np.nanmean(z))
    gy, gx = np.gradient(zz, cell)
    slope = np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    a, b = math.radians(az), math.radians(alt)
    hs = np.sin(b) * np.cos(slope) + np.cos(b) * np.sin(slope) * np.cos(a - aspect)
    hs[~np.isfinite(z)] = np.nan
    return hs


def save_closeup(path_noext, c, tiles, polys, cands, size_m, zoom_m, zfac):
    """Two panels: the ground around the pit, and a close zoom at full resolution where the texture of the floor
    shows: boulders (a talus pit), straight ridges across the wash (fallen logs), a clean rim (a collapse).
    Written as JPEG to keep a county's worth of them small; PNG if JPEG is not available."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.6, 4.2))
    if not draw_closeup(a1, c, tiles, polys, cands, size_m, zfac, title=False):
        plt.close(fig)
        return None
    a1.add_patch(Rectangle((-zoom_m / 2, -zoom_m / 2), zoom_m, zoom_m, fill=False, ec="yellow", lw=0.9))
    draw_closeup(a2, c, tiles, polys, cands, zoom_m, zfac, title=False, step=0.5, neighbors=False)
    a1.set_title(f"{size_m:.0f} m around (yellow box = zoom)", fontsize=7.5)
    a2.set_title(f"{zoom_m:.0f} m zoom: boulders, logs, rims", fontsize=7.5)
    line2 = c.get("flags") or "no flags"
    if c.get("depth_open_ft") is not None and c["depth_ft"] - c["depth_open_ft"] >= 2:
        line2 += f"  |  {c['depth_open_ft']:.0f} ft without thin ridges"
    if c.get("verdict"):
        line2 += f"  |  verdict: {c['verdict']}"
    fig.suptitle(f"{c['id']}   {c['depth_ft']:.0f} ft deep, {c['diam_ft']:.0f} ft across, {c['catch_ac']:.0f} ac "
                 f"draining in\n{line2}", fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    try:
        fig.savefig(path_noext + ".jpg", dpi=100, pil_kwargs={"quality": 82})
        name = os.path.basename(path_noext) + ".jpg"
    except Exception:
        fig.savefig(path_noext + ".png", dpi=100)
        name = os.path.basename(path_noext) + ".png"
    plt.close(fig)
    return name


def draw_closeup(ax, c, tiles, polys, cands, size_m, zfac, title=True, step=None, neighbors=True):
    from rasterio.warp import transform as rtransform
    w = tiles.window(c["lon"], c["lat"], size_m / 2.0, zfac)
    if w is None:
        ax.text(0.5, 0.5, "no elevation data here", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return False
    z, tr, crs, cell, (cx, cy), f = w
    H, W = z.shape
    x0 = (tr.c - cx) * f
    y1 = (tr.f - cy) * f
    ext = (x0, x0 + W * cell, y1 - H * cell, y1)
    ax.imshow(hillshade(z, cell), cmap="gray", vmin=-0.1, vmax=1.0, extent=ext, interpolation="bilinear")
    zr = float(np.nanmax(z) - np.nanmin(z))
    if step is None:
        step = 1.0 if zr <= 40 else (2.0 if zr <= 100 else 5.0)
    lv = np.arange(math.floor(np.nanmin(z)), math.ceil(np.nanmax(z)) + step, step)
    if 2 <= len(lv) <= 400:
        xs = np.linspace(ext[0] + cell / 2, ext[1] - cell / 2, W)
        ys = np.linspace(ext[2] + cell / 2, ext[3] - cell / 2, H)
        zz = np.where(np.isfinite(z), z, np.nan)[::-1]
        thin = [l for i, l in enumerate(lv) if round(l / step) % 5]
        bold = [l for i, l in enumerate(lv) if round(l / step) % 5 == 0]
        if thin:
            ax.contour(xs, ys, zz, levels=thin, colors="#3a6ea5", linewidths=0.35, alpha=0.6)
        if bold:
            ax.contour(xs, ys, zz, levels=bold, colors="#1d3f66", linewidths=0.8, alpha=0.8)
    ring = polys.get(c["id"])
    if ring:
        for r in ring:
            px, py = rtransform("EPSG:4326", crs, [p[0] for p in r], [p[1] for p in r])
            ax.plot([(u - cx) * f for u in px], [(v - cy) * f for v in py], color="red", lw=1.4)
    half = size_m / 2.0
    near = [o for o in cands if o is not c and abs(o["lat"] - c["lat"]) < half / 1e5 * 1.2
            and abs(o["lon"] - c["lon"]) < half / 7.5e4 * 1.2] if neighbors else []
    if near:
        ox, oy = rtransform("EPSG:4326", crs, [o["lon"] for o in near], [o["lat"] for o in near])
        for o, u, v in zip(near, ox, oy):
            du, dv = (u - cx) * f, (v - cy) * f
            if abs(du) < half and abs(dv) < half:
                ax.plot(du, dv, "o", ms=4, mfc="none", mec="orange", mew=1.2)
                ax.annotate(o["id"], (du, dv), xytext=(3, 3), textcoords="offset points", fontsize=6.5,
                            color="orange", weight="bold")
    ax.plot(0, 0, "x", color="yellow", ms=8 if size_m >= 150 else 6, mew=2)
    sb = 50.0 if size_m >= 150 else (20.0 if size_m >= 100 else 10.0)
    ax.plot([-half * 0.9, -half * 0.9 + sb], [-half * 0.9] * 2, color="white", lw=3)
    ax.plot([-half * 0.9, -half * 0.9 + sb], [-half * 0.9] * 2, color="black", lw=1.2)
    ax.text(-half * 0.9, -half * 0.84, f"{sb:.0f} m ({sb * M2FT:.0f} ft)   north up   contours {step:g} m",
            fontsize=6.5, color="black", bbox=dict(fc="white", ec="none", alpha=0.7, pad=1))
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        line2 = c.get("flags") or "no flags"
        if c.get("verdict"):
            line2 += f"  |  verdict: {c['verdict']}"
        ax.set_title(f"{c['id']}   {c['depth_ft']:.0f} ft deep, {c['diam_ft']:.0f} ft across, "
                     f"{c['catch_ac']:.0f} ac\n{line2}", fontsize=8)
    return True


def profile_data(ch, tiles, cands, zfac):
    """Ground elevations down the chain line, and the candidates that sit on it."""
    from rasterio.warp import transform as rtransform
    lons = [p[0] for p in ch["coords"]]
    lats = [p[1] for p in ch["coords"]]
    clon, clat = (min(lons) + max(lons)) / 2.0, (min(lats) + max(lats)) / 2.0
    span = max(1.0, max(abs(max(lats) - min(lats)) * 110574,
                        abs(max(lons) - min(lons)) * 111320 * math.cos(math.radians(clat))))
    w = tiles.window(clon, clat, span / 2.0 + 40.0, zfac)
    if w is None:
        return None
    z, tr, crs, cell, (cx, cy), f = w
    xs, ys = rtransform("EPSG:4326", crs, lons, lats)
    inv = ~tr
    g = []
    for x, y in zip(xs, ys):
        cc, rr = inv * (x, y)
        rr, cc = int(rr), int(cc)
        g.append(z[rr, cc] if 0 <= rr < z.shape[0] and 0 <= cc < z.shape[1] else np.nan)
    s = [0.0]
    for i in range(1, len(xs)):
        s.append(s[-1] + math.hypot(xs[i] - xs[i - 1], ys[i] - ys[i - 1]) * f)
    s = np.array(s)
    g = np.array(g)
    X, Y = np.array(xs), np.array(ys)
    pits = []
    if ch.get("members"):
        # the chain's own pits, placed where the line passes their floors
        mx, my = rtransform("EPSG:4326", crs, [m["lon"] for m in ch["members"]], [m["lat"] for m in ch["members"]])
        for m, u, v in zip(ch["members"], mx, my):
            i = int(np.argmin(np.hypot(X - u, Y - v)))
            pits.append(dict(m, s=float(s[i])))
    else:
        # no member list (an older chain file): candidates within 15 m of the line
        near = [o for o in cands if min(lats) - 0.0003 < o["lat"] < max(lats) + 0.0003
                and min(lons) - 0.0004 < o["lon"] < max(lons) + 0.0004]
        if near:
            ox, oy = rtransform("EPSG:4326", crs, [o["lon"] for o in near], [o["lat"] for o in near])
            for o, u, v in zip(near, ox, oy):
                d = np.hypot(X - u, Y - v) * f
                i = int(np.argmin(d))
                if d[i] <= 15.0:
                    pits.append(dict(id=o["id"], cand=o, s=float(s[i]), floor=o["min_z_m"],
                                     rim=o.get("spill_z_m", o["min_z_m"]), r=o.get("diam_m", 10.0) / 2.0))
    pits.sort(key=lambda p: p["s"])
    return s, g, pits


def cover_stats(s, g, pits):
    gaps = []
    for a, b in zip(pits, pits[1:]):
        lo, hi = a["s"] + a["r"], b["s"] - b["r"]
        sel = (s > lo) & (s < hi) & np.isfinite(g)
        if sel.sum() < 2:
            continue
        line = a["floor"] + (b["floor"] - a["floor"]) * (s[sel] - a["s"]) / max(1e-6, b["s"] - a["s"])
        cov = g[sel] - line
        gaps.append(dict(a=a["id"], b=b["id"], min=float(np.min(cov)), med=float(np.median(cov)),
                         length=float(hi - lo)))
    return gaps


def draw_profile(ax, ch, prof, title=True):
    s, g, pits = prof
    ref = pits[0]["rim"] if pits else float(np.nanmax(g))
    ax.plot(s * M2FT, (g - ref) * M2FT, color="#555", lw=1.2, label="ground down the wash")
    if pits:
        ps = np.array([p["s"] for p in pits]) * M2FT
        pf = (np.array([p["floor"] for p in pits]) - ref) * M2FT
        ax.plot(ps, pf, "--", color="#b5651d", lw=1.2, label="line joining pit floors")
        ax.plot(ps, pf, "v", color="#b5651d", ms=7)
        for p, x, y in zip(pits, ps, pf):
            # the wash line skirts each pit's lowest point; draw the drop from the wash to the floor
            gi = int(np.argmin(np.abs(s - p["s"])))
            top = (np.nanmax(g[max(0, gi - 3):gi + 4]) - ref) * M2FT if np.isfinite(g[max(0, gi - 3):gi + 4]).any() else y
            ax.plot([x, x], [top, y], ":", color="#7a3e00", lw=1.0)
        for p, x, y in zip(pits, ps, pf):
            ax.annotate(p["id"], (x, y), xytext=(0, -12), textcoords="offset points", fontsize=7,
                        ha="center", color="#7a3e00" if p.get("cand") is not None else "#999")
        for gp, a, b in zip(cover_stats(s, g, pits), pits, pits[1:]):
            sel = (s > a["s"] + a["r"]) & (s < b["s"] - b["r"])
            line = a["floor"] + (b["floor"] - a["floor"]) * (s[sel] - a["s"]) / max(1e-6, b["s"] - a["s"])
            ax.fill_between(s[sel] * M2FT, (line - ref) * M2FT, (g[sel] - ref) * M2FT, color="#e8c99b", alpha=0.6)
            xm = (a["s"] + b["s"]) / 2 * M2FT
            ax.annotate(f"{gp['min'] * M2FT:.0f}-{gp['med'] * M2FT:.0f} ft", (xm, (a["floor"] + b["floor"]) / 2 * M2FT
                        - ref * M2FT), xytext=(0, 6), textcoords="offset points", fontsize=6.5, ha="center")
    ax.set_xlabel("feet along the wash from the top pit")
    ax.set_ylabel("feet below the top pit's rim")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper right")
    if title:
        p = ch["props"]
        ax.set_title(f"{ch['cid']}: {p.get('n_pits', len(pits))} pits, {p.get('stem_ft', 0):,} ft along the wash, "
                     f"{p.get('floor_drop_ft', 0):.0f} ft floor to floor\nshaded: ground above the line joining "
                     f"pit floors (least-median), not roof thickness. Preliminary, not a survey.", fontsize=8)


# ----------------------------------------------------------------------------- tags and truth

def parse_tag(spec):
    if "=" not in spec:
        sys.exit(f"--tag wants NAME=path[:FIELD], got {spec}")
    name, rest = spec.split("=", 1)
    field = None
    if ":" in rest:
        head, tail = rest.rsplit(":", 1)
        if tail and not any(ch in tail for ch in "\\/") and not tail.lower().endswith((".shp", ".json", ".geojson")):
            rest, field = head, tail
    return name.strip(), rest.strip(), field


def tag_points(cands, spec, ss):
    from matplotlib.path import Path
    from rasterio.warp import transform as rtransform
    name, path, field = parse_tag(spec)
    geoms, vcrs = ss.read_vector(path)
    if not geoms:
        log(f"tag {name}: no features in {path}")
        return name
    fields = list(geoms[0].get("properties", {}).keys())
    if not field or field not in fields:
        sys.exit(f"--tag {name}: name the attribute to use after a colon, e.g. {name}={path}:{fields[0] if fields else 'FIELD'}\n"
                 f"fields in {os.path.basename(path)}: {', '.join(fields)}")
    if vcrs is None:
        log(f"tag {name}: {path} has no .prj; assuming longitude/latitude")
        px, py = [c["lon"] for c in cands], [c["lat"] for c in cands]
    else:
        px, py = rtransform("EPSG:4326", vcrs, [c["lon"] for c in cands], [c["lat"] for c in cands])
    P = np.column_stack([px, py])
    vals = [None] * len(cands)
    for g in geoms:
        t = g.get("type")
        if t == "Polygon":
            polys = [g["coordinates"]]
        elif t == "MultiPolygon":
            polys = g["coordinates"]
        else:
            continue
        rings = [np.asarray(r, dtype=float)[:, :2] for poly in polys for r in poly if len(r) >= 3]
        if not rings:
            continue
        allr = np.vstack(rings)
        x0, y0 = allr.min(axis=0)
        x1, y1 = allr.max(axis=0)
        sel = np.where((P[:, 0] >= x0) & (P[:, 0] <= x1) & (P[:, 1] >= y0) & (P[:, 1] <= y1))[0]
        sel = [i for i in sel if vals[i] is None]
        if not sel:
            continue
        cnt = np.zeros(len(sel), dtype=int)
        for r in rings:
            cnt += Path(r).contains_points(P[sel]).astype(int)
        for i, k in zip(sel, cnt):
            if k % 2 == 1:
                vals[i] = str(g["properties"].get(field, ""))
    for c, v in zip(cands, vals):
        c[f"tag_{name}"] = v if v is not None else "(none)"
    log(f"tag {name}: {sum(v is not None for v in vals)} of {len(cands)} candidates fall inside a polygon")
    return name


def read_truth(path):
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            r = {k.strip().lower(): (v or "").strip() for k, v in r.items() if k}
            try:
                la, lo = float(r.get("lat", "")), float(r.get("lon", ""))
            except ValueError:
                continue
            if not r.get("verdict"):
                continue
            rows.append(dict(lat=la, lon=lo, verdict=r["verdict"], name=r.get("name", ""),
                             notes=r.get("notes", ""), date=r.get("date", ""), id=r.get("id", "")))
    return rows


def match_truth(cands, rows, match_m):
    def dist(a, b):
        return math.hypot((a["lat"] - b["lat"]) * 110574,
                          (a["lon"] - b["lon"]) * 111320 * math.cos(math.radians(a["lat"])))
    for r in rows:
        best = min(cands, key=lambda c: dist(c, r)) if cands else None
        d = dist(best, r) if best else 1e9
        r["match"], r["match_m"] = (best if d <= match_m else None), d
        r["nearest"] = best
        if r["match"] is not None:
            c = r["match"]
            if not c.get("verdict") or vclass(c["verdict"]) == "unknown":
                c["verdict"], c["truth_name"], c["truth_notes"], c["truth_date"] = r["verdict"], r["name"], r["notes"], r["date"]


def truth_report(cands, rows, out, match_m):
    L = []
    L.append(f"ground truth: {len(rows)} verdict(s); matched to a candidate within {match_m:g} m: "
             f"{sum(1 for r in rows if r['match'] is not None)}\n")
    for lay in ("pits", "offchannel", "bigwash"):
        sub = [c for c in cands if c["layer"] == lay and c.get("verdict") and vclass(c["verdict"]) != "unknown"]
        if not sub:
            continue
        n = len(sub)
        k = {cl: sum(1 for c in sub if vclass(c["verdict"]) == cl) for cl in ("pipe", "real", "false")}
        L.append(f"{lay}: {n} checked: pipe {k['pipe']} ({100 * k['pipe'] / n:.0f}%), other real hole {k['real']}, "
                 f"false {k['false']} ({100 * k['false'] / n:.0f}%)")
        words = {}
        for c in sub:
            w = re.split(r"[\s,;:/()]+", c["verdict"].lower())[0]
            words[w] = words.get(w, 0) + 1
        L.append("   " + ", ".join(f"{w} {n_}" for w, n_ in sorted(words.items(), key=lambda t: -t[1])))
        if lay == "pits":
            for lo, hi in ((1, 25), (26, 50), (51, 100), (101, 200), (201, 10 ** 9)):
                band = [c for c in sub if lo <= c["rank_key"] <= hi]
                if band:
                    p = sum(1 for c in band if vclass(c["verdict"]) == "pipe")
                    rl = sum(1 for c in band if vclass(c["verdict"]) in ("pipe", "real"))
                    L.append(f"   rank {lo}-{hi if hi < 10 ** 9 else 'end'}: {len(band)} checked, "
                             f"pipe {100 * p / len(band):.0f}%, real hole {100 * rl / len(band):.0f}%")
        L.append("")
    known = [r for r in rows if re.split(r"[\s,;:/()]+", r["verdict"].lower())[0] == "known"]
    if known:
        found = [r for r in known if r["match"] is not None]
        L.append(f"KNOWN CAVES (known apart from the program): found {len(found)} of {len(known)}")
        for r in known:
            if r["match"] is not None:
                c = r["match"]
                L.append(f"   found  {r['name'] or '(unnamed)'}: {c['id']} ({c['layer']}), {r['match_m']:.0f} m from the "
                         f"recorded point")
            else:
                nb = r["nearest"]
                L.append(f"   missed {r['name'] or '(unnamed)'}: nearest candidate {nb['id'] if nb else '-'} "
                         f"{r['match_m']:.0f} m away; run sink_scout with --check {r['lat']:.5f},{r['lon']:.5f} "
                         f"on that ground to see why")
    open(os.path.join(out, "truth_report.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


def write_worklist(cands, rows, out, tag_names):
    cols = ["id", "layer", "rank", "score", "depth_ft", "diam_ft", "catch_ac", "chain_n", "flags"] + \
           [f"tag_{t}" for t in tag_names] + ["lat", "lon", "verdict", "name", "notes", "date"]
    with open(os.path.join(out, "ground_truth_worklist.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for c in cands:
            w.writerow([c["id"], c["layer"], c["rank_key"], round(c.get("score", 0), 3), round(c["depth_ft"], 1),
                        round(c["diam_ft"], 1), round(c["catch_ac"], 1), c.get("chain_n", 1), c.get("flags", "")]
                       + [c.get(f"tag_{t}", "") for t in tag_names]
                       + [f"{c['lat']:.7f}", f"{c['lon']:.7f}", c.get("verdict", ""), c.get("truth_name", ""),
                          c.get("truth_notes", ""), c.get("truth_date", "")])
        for r in rows:
            if r["match"] is None:
                w.writerow(["", "(no candidate)", "", "", "", "", "", "", ""] + ["" for _ in tag_names]
                           + [f"{r['lat']:.7f}", f"{r['lon']:.7f}", r["verdict"], r["name"], r["notes"], r["date"]])


# ----------------------------------------------------------------------------- outputs

def kml_color(c, smax):
    cl = vclass(c.get("verdict"))
    if cl == "pipe":
        return "ff3cd23c"
    if cl == "real":
        return "ffd2c83c"
    if cl == "false":
        return "ff3c3cdc"
    t = max(0.0, min(1.0, c.get("score", 0) / smax)) if smax > 0 else 0
    return f"ff{int(40 * (1 - t)):02x}{int(60 + 175 * t):02x}ff"


def balloon(c, tag_names, closeup):
    rows = [("depth", f"{c['depth_ft']:.0f} ft"), ("across", f"{c['diam_ft']:.0f} ft"),
            ("depth/width", f"{c['ratio']:.2f}"), ("catchment", f"{c['catch_ac']:.0f} ac"),
            ("floor", f"{c['min_z_m'] * M2FT:,.0f} ft")]
    if c.get("depth_open_ft") is not None:
        rows.append(("depth without thin ridges", f"{c['depth_open_ft']:.0f} ft"))
    if c.get("chain_cid"):
        rows.append(("chain", f"{c['chain_cid']} ({c.get('chain_n', '?')} pits)"))
    if c.get("road_m") is not None:
        rows.append(("nearest road/track", f"{c['road_m'] * M2FT:,.0f} ft"))
    if (c.get("nsib") or 1) >= 5:
        rows.append(("note", f"{c['nsib']} low spots in one hollow: the boulder-pit signature"))
    for t in tag_names:
        rows.append((t, c.get(f"tag_{t}", "")))
    rows.append(("flags", c.get("flags") or "none"))
    if c.get("verdict"):
        rows.append(("verdict", f"{c['verdict']} {c.get('truth_notes', '')}".strip()))
    img = f'<img src="closeups/{closeup}" width="600"/><br/>' if closeup else ""
    tbl = "".join(f"<tr><td><b>{k}</b></td><td>{v}</td></tr>" for k, v in rows)
    return f"<![CDATA[{img}<table>{tbl}</table>]]>"


def write_kmz(cands, chains, out, tag_names, closeups, profiles):
    titles = {"pits": "Pits on washes", "offchannel": "Off-channel pits", "bigwash": "Pits on big washes"}
    k = ['<?xml version="1.0" encoding="UTF-8"?>', '<kml xmlns="http://www.opengis.net/kml/2.2"><Document>',
         "<name>sink_scout field kit</name>"]
    for lay in ("pits", "offchannel", "bigwash"):
        sub = [c for c in cands if c["layer"] == lay]
        if not sub:
            continue
        smax = max(c.get("score", 0) for c in sub)
        k.append(f"<Folder><name>{titles[lay]} ({len(sub)})</name>")
        for c in sub:
            k.append(f"<Placemark><name>{c['id']}</name><description>{balloon(c, tag_names, closeups.get(c['id']))}"
                     f"</description><Style><IconStyle><color>{kml_color(c, smax)}</color><scale>0.9</scale><Icon>"
                     f"<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href></Icon></IconStyle>"
                     f"</Style><Point><coordinates>{c['lon']:.7f},{c['lat']:.7f},0</coordinates></Point></Placemark>")
        k.append("</Folder>")
    if chains:
        k.append(f"<Folder><name>Chains ({len(chains)})</name>")
        for ch in chains:
            p = ch["props"]
            img = f'<img src="profiles/{profiles[ch["cid"]]}" width="520"/><br/>' if ch["cid"] in profiles else ""
            desc = f"<![CDATA[{img}pits: {' '.join(ch['pits']) or '-'}<br/>{p.get('stem_ft', 0):,} ft along the wash]]>"
            k.append(f"<Placemark><name>{ch['cid']}: {p.get('n_pits', '?')} pits, {p.get('stem_ft', 0):,} ft</name>"
                     f"<description>{desc}</description><Style><LineStyle><color>ff00d7ff</color><width>3</width>"
                     f"</LineStyle></Style><LineString><tessellate>1</tessellate><coordinates>"
                     + " ".join(f"{a:.7f},{b:.7f},0" for a, b in ch["coords"]) + "</coordinates></LineString></Placemark>")
        k.append("</Folder>")
    k.append("</Document></kml>")
    path = os.path.join(out, "field_candidates.kmz")
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("doc.kml", "\n".join(k))
        for name in closeups.values():
            z.write(os.path.join(out, "closeups", name), f"closeups/{name}")
        for name in profiles.values():
            z.write(os.path.join(out, "profiles", name), f"profiles/{name}")
    return path


def write_gpx(cands, chains, out, tag_names):
    def esc(s):
        return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    L = ['<?xml version="1.0" encoding="UTF-8"?>',
         '<gpx version="1.1" creator="sink_scout field_kit" xmlns="http://www.topografix.com/GPX/1/1">']
    for c in cands:
        d = (f"{c['depth_ft']:.0f} ft deep, {c['diam_ft']:.0f} ft across, {c['catch_ac']:.0f} ac"
             + "".join(f", {c.get('tag_' + t)}" for t in tag_names if c.get("tag_" + t) not in (None, "(none)"))
             + (f", verdict {c['verdict']}" if c.get("verdict") else ""))
        L.append(f'<wpt lat="{c["lat"]:.7f}" lon="{c["lon"]:.7f}"><ele>{c.get("spill_z_m", c["min_z_m"]):.1f}</ele>'
                 f"<name>{esc(c['id'])}</name><desc>{esc(d)}</desc><sym>Flag, Blue</sym></wpt>")
    for ch in chains:
        L.append(f"<trk><name>{esc(ch['cid'])} {ch['props'].get('n_pits', '')} pits</name><trkseg>"
                 + "".join(f'<trkpt lat="{b:.7f}" lon="{a:.7f}"/>' for a, b in ch["coords"]) + "</trkseg></trk>")
    L.append("</gpx>")
    open(os.path.join(out, "field_candidates.gpx"), "w").write("\n".join(L))


def write_csv(cands, out, tag_names, closeups):
    cols = ["id", "layer", "rank", "score", "depth_ft", "diam_ft", "ratio", "catch_ac", "chain", "chain_n",
            "floor_ft", "road_ft", "depth_open_ft", "flags"] + [f"tag_{t}" for t in tag_names] + \
           ["verdict", "lat", "lon", "closeup"]
    with open(os.path.join(out, "field_candidates.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for c in cands:
            w.writerow([c["id"], c["layer"], c["rank_key"], round(c.get("score", 0), 3), round(c["depth_ft"], 1),
                        round(c["diam_ft"], 1), round(c["ratio"], 3), round(c["catch_ac"], 1), c.get("chain_cid", ""),
                        c.get("chain_n", 1), round(c["min_z_m"] * M2FT), round(c["road_m"] * M2FT) if c.get("road_m") is not None else "",
                        round(c["depth_open_ft"], 1) if c.get("depth_open_ft") is not None else "", c.get("flags", "")]
                       + [c.get(f"tag_{t}", "") for t in tag_names]
                       + [c.get("verdict", ""), f"{c['lat']:.7f}", f"{c['lon']:.7f}", closeups.get(c["id"], "")])


CHECKLIST = [
    "BEFORE    [  ] rain/snow in the last 14 days checked: ____________    [  ] someone knows where and when: ______________",
    "          [  ] grotto cave-rescue number carried          [  ] three lights each, helmet, mask           ",
    "AT MOUTH  floor:  [  ] dry, cracked into plates    [  ] hard, a heel leaves an edge    [  ] dark, sheen, mud on stones -> STOP",
    "          fall cones in the first 30 ft:  [  ] dusty    [  ] sharp, no dust -> roof is working",
    "PARTY     in: ______________________    at the mouth (does not come in): ______________________",
    "FOUND     [  ] pit only    [  ] entrance    [  ] skylight(s) ____    [  ] choke    passage walked ______ ft  (pace / tape)",
    "          [  ] rodent sign / middens    [  ] bats    [  ] water or recent flow    photos # __________",
    "VERDICT   [  ] pipe   [  ] talus   [  ] pond   [  ] culvert   [  ] logs   [  ] mine / cellar   [  ] other ____________",
    "NOTES     ____________________________________________________________________________________________",
    "          ____________________________________________________________________________________________",
]


def write_sheets(sel, chains_by_id, tiles, polys, cands, out, size_m, zfac, tag_names, profiles_cache):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from rasterio.warp import transform as rtransform
    path = os.path.join(out, "sheets.pdf")
    with PdfPages(path) as pdf:
        for c in sel:
            fig = plt.figure(figsize=(8.5, 11))
            fig.text(0.06, 0.965, f"{c['id']}   {c['layer']}   rank {c['rank_key']}", fontsize=15, weight="bold")
            fig.text(0.06, 0.945, "PRIVATE: a cave location. For the grotto and the people walking it, not for posting.",
                     fontsize=8, color="#a00")
            ax = fig.add_axes([0.06, 0.53, 0.46, 0.40])
            draw_closeup(ax, c, tiles, polys, cands, size_m, zfac, title=False)
            axz = fig.add_axes([0.60, 0.53, 0.30, 0.17])
            draw_closeup(axz, c, tiles, polys, cands, 60.0, zfac, title=False, step=0.5, neighbors=False)
            axz.set_title("60 m zoom", fontsize=7)
            it = tiles.group_for(c["lon"], c["lat"], 10)
            utm = ""
            if it is not None:
                x, y = rtransform("EPSG:4326", it["crs"], [c["lon"]], [c["lat"]])
                epsg = it["crs"].to_epsg()
                utm = f"{x[0]:,.0f} E  {y[0]:,.0f} N  ({'EPSG:' + str(epsg) if epsg else 'tile CRS'})"
            lines = [f"lat, lon   {c['lat']:.6f}, {c['lon']:.6f}", f"grid      {utm}", "",
                     f"floor     {c['min_z_m'] * M2FT:,.0f} ft     rim {c.get('spill_z_m', c['min_z_m']) * M2FT:,.0f} ft",
                     f"depth     {c['depth_ft']:.0f} ft    across {c['diam_ft']:.0f} ft    d/w {c['ratio']:.2f}",
                     f"catchment {c['catch_ac']:.0f} ac"]
            if c.get("depth_open_ft") is not None:
                lines.append(f"depth without thin ridges {c['depth_open_ft']:.0f} ft")
            if c.get("chain_cid"):
                lines.append(f"chain     {c['chain_cid']}, {c.get('chain_n', '?')} pits")
            if c.get("road_m") is not None:
                lines.append(f"road/track {c['road_m'] * M2FT:,.0f} ft away (straight line)")
            for t in tag_names:
                lines.append(f"{t:<9} {c.get('tag_' + t, '')}")
            lines.append(f"flags     {c.get('flags') or 'none'}")
            if (c.get("nsib") or 1) >= 5:
                lines.append(f"note      {c['nsib']} low spots in one hollow (boulder pit?)")
            if c.get("verdict"):
                lines.append(f"verdict   {c['verdict']}  {c.get('truth_notes', '')}")
            fig.text(0.56, 0.92, "\n".join(lines), fontsize=7.6, family="monospace", va="top")
            ch = chains_by_id.get(c.get("chain_cid"))
            if ch is not None and profiles_cache.get(ch["cid"]) is not None:
                axp = fig.add_axes([0.08, 0.30, 0.86, 0.19])
                draw_profile(axp, ch, profiles_cache[ch["cid"]], title=False)
                axp.set_title(f"{ch['cid']}: ground down the wash; shaded = ground above the line joining pit floors "
                              f"(not roof thickness)", fontsize=7.5)
            fig.text(0.06, 0.255, "FIELD CHECKS", fontsize=10, weight="bold")
            fig.text(0.06, 0.24, "\n\n".join(CHECKLIST), fontsize=7.2, family="monospace", va="top")
            fig.text(0.06, 0.02, f"date ____________   weather ____________   sink_scout field kit, "
                                 f"printed {time.strftime('%Y-%m-%d')}", fontsize=7)
            pdf.savefig(fig)
            plt.close(fig)
    return path


# ----------------------------------------------------------------------------- main

def select_sheets(spec, cands, chains):
    if not spec:
        return []
    by_id = {c["id"]: c for c in cands}
    out = []
    for tok in [t.strip() for t in spec.split(",") if t.strip()]:
        m = re.fullmatch(r"top(\d+)", tok, re.IGNORECASE)
        if m:
            out += [c for c in cands if c["layer"] == "pits"][:int(m.group(1))]
        elif tok.upper().startswith("C") and tok[1:].isdigit():
            ch = next((h for h in chains if h["cid"] == tok.upper()), None)
            out += [by_id[i] for i in (ch["pits"] if ch else []) if i in by_id]   # members in the list
        elif tok in by_id:
            out.append(by_id[tok])
        else:
            log(f"--sheets: no candidate {tok}")
    seen, uniq = set(), []
    for c in out:
        if c["id"] not in seen:
            seen.add(c["id"])
            uniq.append(c)
    return uniq


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("candidates", help="county_candidates.geojson or a run's candidates.geojson")
    p.add_argument("--tiles", required=True, help="folder with the elevation tiles the scan used")
    p.add_argument("--chains", help="county_chains.geojson or a run's chains.geojson (profiles, tracks)")
    p.add_argument("--out", default="field")
    p.add_argument("--crs", help="CRS override for tiles with bad metadata (e.g. mesa)")
    p.add_argument("--truth", help="CSV of verdicts: lat, lon, verdict[, name, notes, date]")
    p.add_argument("--match-m", type=float, default=40.0, help="how close a verdict must be to a candidate (m)")
    p.add_argument("--tag", action="append", default=[], help="NAME=polygons.shp:FIELD, e.g. land=sma.shp:ADMIN_AGEN "
                   "(repeatable): label each candidate with the polygon it falls in")
    p.add_argument("--layers", default="pits,offchannel,bigwash", help="which layers get close-ups and waypoints")
    p.add_argument("--closeup-m", type=float, default=300.0, help="close-up window in meters (default 300)")
    p.add_argument("--zoom-m", type=float, default=60.0, help="zoom panel in meters (default 60)")
    p.add_argument("--no-closeups", action="store_true", help="skip the close-ups (fast; KMZ without pictures)")
    p.add_argument("--sheets", help="field sheets: top25, a chain (C3), or ids (P31,P36,O4), comma-separated")
    args = p.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    try:
        import sink_scout as ss
    except ImportError:
        sys.exit("field_kit.py needs sink_scout.py in the same folder")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(args.out, exist_ok=True)
    cands, polys = load_candidates(args.candidates)
    want = set(args.layers.split(","))
    cands = [c for c in cands if c["layer"] in want]
    log(f"{len(cands)} candidates from {os.path.basename(args.candidates)}")
    chains = load_chains(args.chains, cands)
    tiles = Tiles(args.tiles, args.crs, ss)
    zfac = detect_zfac(tiles, cands)
    if zfac != 1.0:
        log("tiles are in feet; converting")

    tag_names = [tag_points(cands, spec, ss) for spec in args.tag]
    rows = []
    if args.truth:
        rows = read_truth(args.truth)
        match_truth(cands, rows, args.match_m)
        log(f"{len(rows)} verdicts read, {sum(1 for r in rows if r['match'] is not None)} matched")

    # chains: profiles, membership
    prof_cache, profiles = {}, {}
    if chains:
        os.makedirs(os.path.join(args.out, "profiles"), exist_ok=True)
        rows_c = []
        for ch in chains:
            pr = profile_data(ch, tiles, cands, zfac)
            prof_cache[ch["cid"]] = pr
            if pr is None:
                continue
            s, g, pits = pr
            ch["pits"] = [q["id"] for q in pits]
            for q in pits:
                if q.get("cand") is not None:
                    q["cand"].setdefault("chain_cid", ch["cid"])
            gaps = cover_stats(s, g, pits)
            fig, ax = plt.subplots(figsize=(7.5, 3.4))
            draw_profile(ax, ch, pr)
            fig.tight_layout()
            name = f"{ch['cid']}.png"
            fig.savefig(os.path.join(args.out, "profiles", name), dpi=110)
            plt.close(fig)
            profiles[ch["cid"]] = name
            cp = ch["props"]
            rows_c.append([ch["cid"], cp.get("n_pits", len(pits)), " ".join(ch["pits"]), cp.get("stem_ft", ""),
                           cp.get("floor_drop_ft", ""), cp.get("grade_pct", ""),
                           round(min(gp["min"] for gp in gaps) * M2FT, 1) if gaps else "",
                           round(statistics.median([gp["med"] for gp in gaps]) * M2FT, 1) if gaps else ""])
        with open(os.path.join(args.out, "profiles", "chains.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["chain", "n_pits", "pits_on_line", "stem_ft", "floor_drop_ft", "grade_pct",
                        "least_cover_ft", "median_cover_ft"])
            w.writerows(rows_c)
        log(f"{len(profiles)} chain profiles drawn")

    closeups = {}
    if not args.no_closeups:
        os.makedirs(os.path.join(args.out, "closeups"), exist_ok=True)
        t0 = time.time()
        for i, c in enumerate(cands, 1):
            name = save_closeup(os.path.join(args.out, "closeups", c["id"].replace("#", "r")), c, tiles, polys,
                                cands, args.closeup_m, args.zoom_m, zfac)
            if name:
                closeups[c["id"]] = name
            if i % 100 == 0:
                log(f"  {i}/{len(cands)} close-ups ({(time.time() - t0) / i:.2f} s each)")
        log(f"{len(closeups)} close-ups drawn")

    kmz = write_kmz(cands, chains, args.out, tag_names, closeups, profiles)
    write_gpx(cands, chains, args.out, tag_names)
    write_csv(cands, args.out, tag_names, closeups)
    write_worklist(cands, rows, args.out, tag_names)
    if rows:
        truth_report(cands, rows, args.out, args.match_m)
    for t in tag_names:
        pits = [c for c in cands if c["layer"] == "pits"]
        cnt = {}
        for c in pits:
            cnt[c.get(f"tag_{t}")] = cnt.get(c.get(f"tag_{t}"), 0) + 1
        log(f"pits by {t}: " + ", ".join(f"{k} {v}" for k, v in sorted(cnt.items(), key=lambda kv: -kv[1])[:10]))
    sel = select_sheets(args.sheets, cands, chains)
    if sel:
        sp = write_sheets(sel, {h["cid"]: h for h in chains}, tiles, polys, cands, args.out, args.closeup_m, zfac,
                          tag_names, prof_cache)
        log(f"{len(sel)} field sheet(s) -> {sp}")
    log(f"done -> {kmz}")


if __name__ == "__main__":
    main()
