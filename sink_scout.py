#!/usr/bin/env python3
"""
sink_scout.py — find piping-cave collapse sinkholes in bare-earth lidar DEMs.

WHAT IT LOOKS FOR
    A piping cave (pseudokarst pipe in colluvium / bentonitic fill) betrays itself
    at the surface as a chain of steep-sided collapse sinkholes sitting ON the axis
    of a wash, with real catchment draining into them. A normal wash has no closed
    depressions at all, so anything that survives the filters below is interesting
    by definition.

METHOD
    1. Depression fill (grayscale reconstruction, Vincent 1993) of the DEM.
       fill - dem  = depth of every closed depression, pixel by pixel.
    2. Flats produced by the fill are given a micro-gradient toward their outlet
       (geodesic distance inside the flat), then D8 flow direction + accumulation
       are computed on the conditioned surface, so drainage runs THROUGH each basin
       and the accumulation inside a basin equals the catchment feeding it.
    3. Footprints (fill depth > --floor) are broken into individual pits: h-minima
       markers find every pit with its own depth >= --min-depth, including pits
       hiding inside a larger basin (a wash dammed by a road fill, two sinks
       sharing a pour point); watershed splits them and each pit's spill level is
       the saddle where it would overflow, so its depth is its own, not the dam's.
       A flood-fill from the bottom then trims each pit to its rim (where the wall
       slope drops below --rim-slope), so a sink does not inherit the footprint of
       the wash floor around it.
    4. Each pit is measured: depth, area, equivalent diameter, depth/width ratio
       ("pit shape"), elongation, catchment (ha / acres), whether a channel
       (acc >= --channel-ha) enters it, distance to the nearest road, geology
       membership, and chain membership (other pits within --chain-dist).
    5. Thin-dam test: every ridge narrower than --dam-width (default 4 m) is shaved off
       around each pit and the pit is measured again. A hollow held by a fallen log
       across a ravine, a fence line, a thin berm or a lidar bridge over a narrow
       channel drains away; a collapse pit, with ground all round it, keeps its depth.
       Pits that lose more than half their depth are flagged thin-dam and kept out of
       chains, so a staircase of logs down a forested valley cannot pose as a pipe.
    6. Hard filters (depth, area, edge contact) then a weighted score: depth, pit
       shape, log catchment and chain bonus add; off-channel, shallow/wide, creek, road,
       thin-dam, off-geology and nodata-edge multiply the score down. Ranked and written out.

    Scale: ~25 s and ~1.6 GB for 9 million cells (3000 x 3000). A 4-ft county tile
    runs in seconds; mosaic 10-20 tiles at a time, or use --agg 2 for very large areas.

INPUT
    One or more projected DEMs (GeoTIFF, ERDAS .img, anything GDAL reads).
    Elevation and the CRS must be in linear units. If the file's CRS tag is wrong
    (Mesa County tiles say UTM 13 but are UTM 12 NAD83 in US survey feet) override
    it:  --crs mesa   (or any EPSG code / proj string).  --z-units ft if the
    elevations are feet (default for --crs mesa).

OUTPUT  (in --out)
    candidates.shp/.shx/.dbf/.prj   ranked points at the deepest cell (DEM CRS)
    depressions.shp/...             footprint polygons, same attributes
    candidates.csv                  all attributes, lat/lon, feet and meters
    offchannel_pits.csv             pit-shaped holes with no drainage: not pipes, but
                                    cellar holes / prospects / cisterns show up here
    candidates.geojson              EPSG:4326, points + footprints
    candidates.kml                  Google Earth placemarks + outlines
    chains.csv / .kml / .geojson    pit chains linked along one drainage: a line down the wash from
                                    the top pit to the bottom pit with the wash length and drop between them
    refine_commands.txt             after a coarse --agg sweep: one full-resolution --bbox command per
                                    chain and per top candidate (coarse cells find the fields and the
                                    long chains; 1-m cells are what separate small pits from ponds).
                                    Catchments, the creek flag and chains are computed inside the box,
                                    so widen --refine-margin (default 300 m) for pits on big washes.
    overview.png                    hillshade, channels, ranked candidates
    thumbs/rank_###.png             close-ups of the top --thumbs candidates
    dem.tif, fill_depth.tif, flowacc_ha.tif, labels.tif   (--save-rasters; dem.tif is the ground as
                                    processed, the file to send when a single pit needs a second opinion)

EXAMPLES
    python sink_scout.py --selftest --out selftest_out
    python sink_scout.py tile.img --crs mesa --out park --check 39.08,-108.54
        (--check reports the rank of a known sinkhole, or why it was dropped)
    python sink_scout.py *.img --crs mesa --mosaic --osm-roads --out gunnison_sweep --thumbs 40
    python sink_scout.py tile_1m.tif --bbox 39.075,-108.545,39.085,-108.530 --out park_1m
        (--bbox clips first, so a 1-m tile can run at full resolution over one hillside)
    python sink_scout.py dem.tif --geology co_geol.shp --geo-field UNIT --geo-units TRc,Jm --out run1

REQUIREMENTS
    pip install numpy scipy scikit-image rasterio pyshp matplotlib
"""

import argparse
import csv
import json
import math
import os
import sys
import time
import urllib.parse
import urllib.request

import numpy as np
from scipy import ndimage as ndi

# ----------------------------------------------------------------------------- utils

S8 = np.ones((3, 3), bool)
OFFS = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
DI = np.array([o[0] for o in OFFS], np.int8)
DJ = np.array([o[1] for o in OFFS], np.int8)
DIST = np.array([math.hypot(*o) for o in OFFS])
M2FT = 1.0 / 0.3048
HA2AC = 2.4710538
CRS_SHORTCUTS = {
    # Mesa County CO 2016 QL2 deliverables: UTM 12 NAD83, US survey feet, elevations in feet
    "mesa": "+proj=utm +zone=12 +datum=NAD83 +units=us-ft +no_defs",
}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def need(modname, pipname=None):
    try:
        return __import__(modname)
    except ImportError:
        sys.exit(f"Missing dependency '{modname}'. Install with: pip install {pipname or modname}")


def shifted(arr, di, dj, fill):
    """arr shifted so that out[i,j] = arr[i+di, j+dj]; out-of-bounds filled."""
    H, W = arr.shape
    out = np.full_like(arr, fill)
    si = slice(max(0, di), H + min(0, di))
    sj = slice(max(0, dj), W + min(0, dj))
    ti = slice(max(0, -di), H + min(0, -di))
    tj = slice(max(0, -dj), W + min(0, -dj))
    out[ti, tj] = arr[si, sj]
    return out


# ----------------------------------------------------------------------------- loading


def crs_from_arg(text):
    from rasterio.crs import CRS
    return CRS.from_user_input(CRS_SHORTCUTS.get(text.lower(), text))


def linear_factor(crs):
    """meters per CRS linear unit."""
    if crs.is_geographic:
        sys.exit("DEM is in geographic (lat/lon) coordinates. Reproject to UTM first, or pass --crs.")
    try:
        return crs.linear_units_factor[1]
    except Exception:
        name = (crs.linear_units or "metre").lower()
        if "survey" in name:
            return 0.3048006096012192
        if "foot" in name or "feet" in name or name == "ft":
            return 0.3048
        return 1.0


def block_aggregate(arr, valid, n):
    """Block-mean downsample by integer factor n, ignoring nodata."""
    H, W = arr.shape
    H2, W2 = H // n, W // n
    a = np.where(valid, arr, 0.0)[: H2 * n, : W2 * n].reshape(H2, n, W2, n)
    v = valid[: H2 * n, : W2 * n].reshape(H2, n, W2, n)
    cnt = v.sum(axis=(1, 3))
    s = a.sum(axis=(1, 3))
    out = np.where(cnt > 0, s / np.maximum(cnt, 1), np.nan)
    return out, cnt > (n * n) // 2


def parse_bbox(text):
    """'lat1,lon1,lat2,lon2' (two opposite corners, any order) -> (south, west, north, east)."""
    try:
        a, b, c, d = [float(v) for v in text.split(",")]
    except ValueError:
        sys.exit("--bbox wants lat1,lon1,lat2,lon2 in decimal degrees, e.g. --bbox 39.075,-108.545,39.085,-108.530")
    return min(a, c), min(b, d), max(a, c), max(b, d)


def load_dem(paths, crs_override=None, mosaic=False, nodata=None, zmin=-1000.0, z_units="m", agg=1, bbox=None):
    rasterio = need("rasterio")
    from rasterio.transform import Affine
    from rasterio.crs import CRS
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds, Window

    srcs = [rasterio.open(p) for p in paths]
    src_nodata = srcs[0].nodata
    crs = crs_from_arg(crs_override) if crs_override else srcs[0].crs
    if crs is None:
        sys.exit("DEM carries no CRS. Pass --crs (e.g. --crs EPSG:26912 or --crs mesa).")

    bounds = None
    if bbox:
        from rasterio.windows import bounds as win_bounds
        s_, w_, n_, e_ = parse_bbox(bbox)
        raw = transform_bounds(CRS.from_epsg(4326), crs, w_, s_, e_, n_)  # (left, bottom, right, top) in DEM CRS
        # snap to the first source's pixel grid so merge/window reads need no resampling
        t0 = srcs[0].transform
        win0 = from_bounds(*raw, transform=t0)
        win0 = Window(math.floor(win0.col_off), math.floor(win0.row_off),
                      math.ceil(win0.col_off + win0.width) - math.floor(win0.col_off),
                      math.ceil(win0.row_off + win0.height) - math.floor(win0.row_off))
        bounds = win_bounds(win0, t0)
        log(f"clipping to --bbox: lat {s_:.5f}..{n_:.5f}, lon {w_:.5f}..{e_:.5f}")

    if len(srcs) > 1:
        if not mosaic:
            sys.exit("Several DEMs given: add --mosaic to merge them, or run one at a time.")
        from rasterio.merge import merge
        # Merge into NaN for float DEMs. Some tiles declare their nodata as the float32 minimum
        # (-3.4028234663852886e+38); handed to merge as the destination nodata, that value fails
        # rasterio's representability check, the destination is left zero-filled and nothing is
        # copied in: a "100% valid" DEM of z = 0. NaN sidesteps the check, and merge still masks each
        # source by its own nodata. Integer DEMs keep their declared nodata when it is in range.
        dt0 = np.dtype(srcs[0].dtypes[0])
        if np.issubdtype(dt0, np.floating):
            merge_nd = float("nan")
        elif src_nodata is not None and np.isfinite(src_nodata) and \
                np.iinfo(dt0).min <= src_nodata <= np.iinfo(dt0).max:
            merge_nd = src_nodata
        else:
            merge_nd = None
        kw = dict(bounds=bounds) if bounds else {}
        if merge_nd is not None:
            kw["nodata"] = merge_nd
        arr, transform = merge(srcs, **kw)
        arr = arr[0]
        src_nodata = merge_nd
    else:
        src = srcs[0]
        if bounds:
            win = from_bounds(*bounds, transform=src.transform)
            try:
                win = win.intersection(Window(0, 0, src.width, src.height)).round_offsets().round_lengths()
            except Exception:
                sys.exit("--bbox does not overlap this DEM")
            if win.width <= 0 or win.height <= 0:
                sys.exit("--bbox does not overlap this DEM")
            arr = src.read(1, window=win)
            transform = src.window_transform(win)
        else:
            arr = src.read(1)
            transform = src.transform
    for s in srcs:
        s.close()
    if nodata is None:
        nodata = src_nodata

    arr = arr.astype(np.float64)
    valid = np.isfinite(arr) & (arr > zmin)
    if nodata is not None and np.isfinite(nodata):
        valid &= arr != nodata
    if agg > 1:
        arr, valid = block_aggregate(arr, valid, agg)
        transform = transform * Affine.scale(agg, agg)
    if z_units == "ft":
        arr = arr * 0.3048006096012192
    arr[~valid] = np.nan

    f = linear_factor(crs)
    cell_x, cell_y = abs(transform.a) * f, abs(transform.e) * f
    if abs(cell_x - cell_y) > 0.01 * cell_x:
        log(f"warning: non-square cells {cell_x:.3f} x {cell_y:.3f} m; using the mean")
    cell = 0.5 * (cell_x + cell_y)
    log(f"DEM {arr.shape[1]} x {arr.shape[0]} cells, {cell:.2f} m/cell, "
        f"{valid.mean()*100:.1f}% valid, z {np.nanmin(arr):.1f}..{np.nanmax(arr):.1f} m, CRS: {crs.to_string()[:60]}")
    return arr, valid, transform, crs, cell, f


# ----------------------------------------------------------------------------- hydrology


def fill_depressions(dem, valid):
    """Priority-flood equivalent via grayscale reconstruction. Returns filled surface (nodata -> very low)."""
    from skimage.morphology import reconstruction

    LOW = -1.0e9
    z = np.where(valid, dem, LOW)
    inv_p = np.pad(~valid, 1, constant_values=True)
    adj_inv = ndi.binary_dilation(inv_p, structure=S8)[1:-1, 1:-1]
    border = valid & adj_inv  # cells that can drain off the grid / into nodata
    seed = np.where(border | ~valid, z, np.nanmax(dem) + 10.0)
    filled = reconstruction(seed, z, method="erosion")
    return filled


def resolve_flats(filled, valid, eps):
    """Add eps * (geodesic distance to outlet + 1) inside flats so every cell has a strictly lower neighbor.
    eps is capped so the added gradient can never exceed the smallest real elevation step in the DEM."""
    z = np.where(valid, filled, -np.inf)
    minnb = np.full_like(z, np.inf)
    minpos = np.inf
    with np.errstate(invalid="ignore"):
        for di, dj in OFFS:
            nz = shifted(z, di, dj, -np.inf)
            np.minimum(minnb, nz, out=minnb)
            d = z - nz
            pos = d[valid & np.isfinite(d) & (d > 0)]
            if pos.size:
                minpos = min(minpos, float(pos.min()))
    flat = valid & ~(minnb < z)
    zeps = z.copy()
    n_flat = int(flat.sum())
    if n_flat == 0:
        return zeps, 0
    nonflat = valid & ~flat
    exits = np.zeros_like(flat)
    for di, dj in OFFS:
        nz = shifted(z, di, dj, -np.inf)
        nnf = shifted(nonflat, di, dj, False)
        exits |= flat & nnf & (nz <= z)

    from skimage.graph import MCP_Geometric
    cost = np.where(flat, 1.0, 10.0)
    mcp = MCP_Geometric(cost, fully_connected=True)
    starts = np.argwhere(exits)
    if len(starts) == 0:
        return zeps, n_flat
    try:
        dist, _ = mcp.find_costs(starts, max_step_cost=5.0)
    except TypeError:  # scikit-image < 0.26
        dist, _ = mcp.find_costs(starts, max_cost=5.0)
    reach = flat & np.isfinite(dist)
    dmax = float(dist[reach].max()) if reach.any() else 0.0
    if np.isfinite(minpos):
        eps = min(eps, 0.5 * minpos / (dmax + 2.0))
    zeps[reach] += eps * (dist[reach] + 1.0)
    stuck = int((flat & ~reach).sum())
    return zeps, stuck


def d8_receivers(zeps, valid):
    """Steepest-descent receiver (flat index) for every valid cell; -1 = drains off grid / into nodata."""
    H, W = zeps.shape
    best = np.full((H, W), -np.inf)
    bk = np.zeros((H, W), np.int8)
    with np.errstate(invalid="ignore"):
      for k, (di, dj) in enumerate(OFFS):
        slope = (zeps - shifted(zeps, di, dj, -np.inf)) / DIST[k]
        upd = slope > best
        best[upd] = slope[upd]
        bk[upd] = k
    itype = np.int32 if H * W < 2**31 - 1 else np.int64
    ii, jj = np.indices((H, W), dtype=itype)
    ri = ii + DI[bk]
    rj = jj + DJ[bk]
    inside = (ri >= 0) & (ri < H) & (rj >= 0) & (rj < W)
    rec = np.full((H, W), -1, itype)
    ok = valid & inside
    rec[ok] = ri[ok] * W + rj[ok]
    rv = np.zeros((H, W), bool)
    rv[ok] = valid.ravel()[rec[ok]]
    rec[ok & ~rv] = -1
    bad = ok & rv & ~(best > 0)
    rec[bad] = -1
    return rec, int(bad.sum())


def flow_accumulation(rec, valid):
    """Vectorised Kahn topological accumulation (cells). No sort, no python per-cell loop."""
    n = rec.size
    recf = rec.ravel()
    has = recf >= 0
    indeg = np.bincount(recf[has], minlength=n).astype(np.int32)
    acc = valid.ravel().astype(np.float64)
    frontier = np.flatnonzero((indeg == 0) & valid.ravel())
    rounds = 0
    while frontier.size:
        rounds += 1
        r = recf[frontier]
        m = r >= 0
        r = r[m]
        if r.size == 0:
            break
        np.add.at(acc, r, acc[frontier][m])
        np.subtract.at(indeg, r, 1)
        frontier = np.unique(r[indeg[r] == 0])
    leftover = int(((indeg > 0) & valid.ravel()).sum())
    if leftover:
        log(f"warning: {leftover} cells never resolved in flow routing (cycles?)")
    return acc.reshape(rec.shape), rounds


# ----------------------------------------------------------------------------- vectors in


def read_vector(path):
    """Return (list of GeoJSON-like geometries, CRS or None). Shapefile via pyshp, else GeoJSON."""
    from rasterio.crs import CRS
    ext = os.path.splitext(path)[1].lower()
    geoms, crs = [], None
    if ext == ".shp":
        shapefile = need("shapefile", "pyshp")
        r = shapefile.Reader(path)
        for shp, rec in zip(r.iterShapes(), r.iterRecords()):
            gj = shp.__geo_interface__
            gj["properties"] = rec.as_dict()
            geoms.append(gj)
        prj = os.path.splitext(path)[0] + ".prj"
        if os.path.exists(prj):
            with open(prj) as f:
                crs = CRS.from_wkt(f.read())
    else:
        with open(path) as f:
            gj = json.load(f)
        feats = gj.get("features", [gj] if gj.get("type") != "FeatureCollection" else [])
        for ft in feats:
            g = ft.get("geometry", ft)
            if g:
                g = dict(g)
                g["properties"] = ft.get("properties", {})
                geoms.append(g)
        crs = CRS.from_epsg(4326)
    return geoms, crs


def reproject_geom(geom, src_crs, dst_crs):
    from rasterio.warp import transform as rtransform
    if src_crs is None or src_crs == dst_crs:
        return geom

    def conv(coords):
        if not coords:
            return coords
        if isinstance(coords[0], (int, float)):
            xs, ys = rtransform(src_crs, dst_crs, [coords[0]], [coords[1]])
            return [xs[0], ys[0]]
        if isinstance(coords[0][0], (int, float)):
            xs = [c[0] for c in coords]
            ys = [c[1] for c in coords]
            tx, ty = rtransform(src_crs, dst_crs, xs, ys)
            return [[x, y] for x, y in zip(tx, ty)]
        return [conv(c) for c in coords]

    g = dict(geom)
    g["coordinates"] = conv(geom["coordinates"])
    return g


def rasterize_distance(geoms, shape, transform, cell):
    """Rasterize geometries (DEM CRS) and return distance-to-nearest in meters."""
    from rasterio import features
    shapes = [(g, 1) for g in geoms if g.get("type") in
              ("LineString", "MultiLineString", "Polygon", "MultiPolygon")]
    if not shapes:
        return None
    mask = features.rasterize(shapes, out_shape=shape, transform=transform, fill=0,
                              all_touched=True, dtype="uint8").astype(bool)
    if not mask.any():
        return None
    return ndi.distance_transform_edt(~mask) * cell


def fetch_osm_roads(crs, transform, shape, cache_path):
    """Download every OSM highway= way in the DEM bbox (Overpass) and return geometries in DEM CRS."""
    from rasterio.crs import CRS
    from rasterio.warp import transform_bounds
    from rasterio.transform import array_bounds
    if os.path.exists(cache_path):
        log(f"using cached OSM roads {cache_path}")
        geoms, vcrs = read_vector(cache_path)
        return [reproject_geom(g, vcrs, crs) for g in geoms]
    b = array_bounds(shape[0], shape[1], transform)
    w, s, e, n = transform_bounds(crs, CRS.from_epsg(4326), *b)
    pad = 0.002
    q = f'[out:json][timeout:90];way["highway"]({s-pad},{w-pad},{n+pad},{e+pad});out geom;'
    body = urllib.parse.urlencode({"data": q}).encode()
    headers = {"User-Agent": "sink_scout/1.0 (piping-cave sinkhole detector; python-urllib)",
               "Accept": "application/json", "Content-Type": "application/x-www-form-urlencoded"}
    log("fetching roads from OpenStreetMap (Overpass)...")
    data = None
    for host in ("https://overpass-api.de/api/interpreter", "https://overpass.kumi.systems/api/interpreter",
                 "https://lz4.overpass-api.de/api/interpreter"):
        try:
            req = urllib.request.Request(host, data=body, headers=headers, method="POST")
            with urllib.request.urlopen(req, timeout=180) as resp:
                data = json.load(resp)
            break
        except Exception as ex:
            log(f"  {host.split('/')[2]}: {ex}")
    if data is None:
        log("warning: OSM fetch failed on every mirror; continuing without roads")
        return []
    feats = []
    for el in data.get("elements", []):
        pts = el.get("geometry")
        if not pts or len(pts) < 2:
            continue
        feats.append({"type": "Feature",
                      "properties": {"highway": el.get("tags", {}).get("highway", ""),
                                     "name": el.get("tags", {}).get("name", "")},
                      "geometry": {"type": "LineString",
                                   "coordinates": [[p["lon"], p["lat"]] for p in pts]}})
    with open(cache_path, "w") as f:
        json.dump({"type": "FeatureCollection", "features": feats}, f)
    log(f"OSM: {len(feats)} road/track ways, cached to {cache_path}")
    return [reproject_geom(f["geometry"], CRS.from_epsg(4326), crs) for f in feats]


# ----------------------------------------------------------------------------- measurement


def flood_to_rim(zs, mask, rc, min_z, spill, cell, slope_thresh):
    """Raise a water level from the bottom of a pit and stop where the wall flattens out: that is the rim.
    Keeps a sink from inheriting the footprint of the wash floor or dammed basin it sits in.
    Returns (component mask at the rim, rim level)."""
    total = spill - min_z
    fallback = mask & (zs < spill)
    if total <= 0 or slope_thresh <= 0:
        return fallback, spill
    step = min(0.5, max(0.1, total / 40.0))
    levels = np.arange(min_z + step, spill, step)
    levels = np.append(levels, spill)
    prev_R, prev_level, prev_comp, prev_A = 0.0, min_z, None, 0
    for lev in levels:
        lab_, _ = ndi.label(mask & (zs < lev), structure=S8)
        li = lab_[rc]
        if li == 0:
            continue
        comp = lab_ == li
        A = int(comp.sum())
        R = math.sqrt(A * cell * cell / math.pi)
        dR = R - prev_R
        if prev_comp is not None and prev_A >= 6 and dR >= 2 * cell and (lev - prev_level) / dR < slope_thresh:
            # rim lies between prev_level (still in the pit) and lev (flooded out): bisect to pin it down
            lo, hi, lo_comp, lo_R = prev_level, lev, prev_comp, prev_R
            for _ in range(4):
                mid = 0.5 * (lo + hi)
                lab_m, _ = ndi.label(mask & (zs < mid), structure=S8)
                comp_m = lab_m == lab_m[rc]
                R_m = math.sqrt(int(comp_m.sum()) * cell * cell / math.pi)
                if R_m - lo_R >= cell and (mid - lo) / (R_m - lo_R) < slope_thresh:
                    hi = mid
                else:
                    lo, lo_comp, lo_R = mid, comp_m, R_m
            return lo_comp, lo
        prev_R, prev_level, prev_comp, prev_A = R, lev, comp, A
    return (prev_comp if prev_comp is not None else fallback), spill


def split_nested(dem, valid, filled, lab, foot, cell, args):
    """Turn fill footprints into individual pits.
    1. h-minima markers (dynamic >= --min-depth) find every pit, including pits hiding inside a larger
       basin (a road-dammed wash, two sinks sharing a pour point).
    2. Footprints holding several markers are split by watershed; each piece's spill level is the saddle
       where it would overflow into its neighbor, so its depth is its OWN depth, not the dam's.
    3. Each piece is trimmed to its rim with flood_to_rim.
    Returns (reg, spill_lut, parent_lut, nsib_lut): region label image and per-region lookup arrays."""
    from skimage.morphology import h_minima
    from skimage.segmentation import watershed

    n = int(lab.max())
    zlow = np.where(valid, dem, -1.0e9)
    hm = (h_minima(zlow, args.min_depth * 0.999) > 0) & foot
    mk, nm = ndi.label(hm, structure=S8)
    empty = (np.zeros_like(lab), np.zeros(1), np.zeros(1, np.int32), np.zeros(1, np.int32))
    if nm == 0:
        return empty
    mk_parent = ndi.maximum(lab, mk, np.arange(1, nm + 1)).astype(np.int64)
    counts = np.bincount(mk_parent, minlength=n + 1)
    objs = ndi.find_objects(lab)

    reg = np.zeros_like(lab)
    spill, parent, nsib = [0.0], [0], [0]
    nxt = 1
    n_split = 0
    for p in np.flatnonzero(counts >= 1):
        sl = objs[p - 1]
        sub = lab[sl] == p
        zs = zlow[sl]
        p_spill = float(filled[sl][sub].max())
        if counts[p] == 1:
            pieces = [(sub, p_spill)]
        else:
            n_split += 1
            ws = watershed(zs, markers=np.where(sub, mk[sl], 0), mask=sub, connectivity=2)
            ids = np.unique(ws[ws > 0])
            saddle = {int(i): p_spill for i in ids}
            for di, dj in OFFS:
                nws = shifted(ws, di, dj, 0)
                nz = shifted(zs, di, dj, np.inf)
                pair = (ws > 0) & (nws > 0) & (nws != ws)
                if not pair.any():
                    continue
                hgt = np.maximum(zs, nz)[pair]
                for i, h in zip(ws[pair], hgt):
                    if h < saddle[int(i)]:
                        saddle[int(i)] = float(h)
            pieces = [((ws == i) & (zs < saddle[int(i)]), saddle[int(i)]) for i in ids]
        for m, sp in pieces:
            if not m.any():
                continue
            zm = np.where(m, zs, np.inf)
            rc = np.unravel_index(int(np.argmin(zm)), zm.shape)
            comp, rim = flood_to_rim(zs, m, rc, float(zm[rc]), sp, cell, args.rim_slope)
            if not comp.any():
                continue
            view = reg[sl]
            view[comp & (view == 0)] = nxt
            spill.append(float(rim))
            parent.append(int(p))
            nsib.append(len(pieces))
            nxt += 1
    if n_split:
        log(f"{n_split} compound footprints split into component pits (h-minima / watershed)")
    if nxt == 1:
        return empty
    return reg, np.array(spill), np.array(parent, np.int32), np.array(nsib, np.int32)


def measure_depressions(dem, valid, filled, acc_ha, transform, cell, args, road_dist=None, geo_dist=None):
    from skimage.measure import regionprops_table

    depth = np.where(valid, filled - np.where(valid, dem, 0.0), 0.0)
    foot = depth > args.floor
    lab, n = ndi.label(foot, structure=S8)
    log(f"{n} closed depressions deeper than the {args.floor} m floor")
    if n == 0:
        return lab, lab, [], depth

    reg, spill_lut, parent_lut, nsib_lut = split_nested(dem, valid, filled, lab, foot, cell, args)
    nr = int(reg.max())
    if nr == 0:
        log(f"none reach the {args.min_depth} m minimum depth")
        return lab, reg, [], depth
    idx = np.arange(1, nr + 1)
    # per-region depth relative to the region's own spill level
    dreg = np.where(reg > 0, spill_lut[reg] - np.where(valid, dem, 0.0), 0.0)
    area_m2 = ndi.sum(reg > 0, reg, idx) * cell * cell
    maxd = ndi.maximum(dreg, reg, idx)
    H, W = reg.shape
    objs = ndi.find_objects(reg)
    edge = np.array([(o[0].start == 0 or o[1].start == 0 or o[0].stop == H or o[1].stop == W)
                     if o is not None else False for o in objs])
    keep = (maxd >= args.min_depth) & (maxd <= args.max_depth) & \
           (area_m2 >= args.min_area) & (area_m2 <= args.max_area)
    if not args.keep_edge:
        keep &= ~edge
    klabels = idx[keep]
    log(f"{len(klabels)} pass depth {args.min_depth}-{args.max_depth} m, area {args.min_area:g}-{args.max_area:g} m2"
        + ("" if args.keep_edge else ", not touching the raster edge"))
    if len(klabels) == 0:
        return lab, np.zeros_like(reg), [], depth

    lut = np.zeros(nr + 1, np.int32)
    lut[klabels] = np.arange(1, len(klabels) + 1)
    lab2 = lut[reg]
    k = np.arange(1, len(klabels) + 1)

    demf = np.where(valid, dem, np.inf)
    meand = ndi.mean(dreg, lab2, k)
    vol = ndi.sum(dreg, lab2, k) * cell * cell
    minpos = ndi.minimum_position(demf, lab2, k)
    com = ndi.center_of_mass(lab2 > 0, lab2, k)
    maxacc = ndi.maximum(acc_ha, lab2, k)
    channel = acc_ha >= args.channel_ha
    channel_d = ndi.binary_dilation(channel, structure=S8)
    onchan = ndi.maximum(channel_d.astype(np.uint8), lab2, k) > 0
    nod = ndi.binary_dilation(~valid, structure=S8)
    near_nodata = ndi.maximum(nod.astype(np.uint8), lab2, k) > 0
    props = regionprops_table(lab2, properties=("label", "axis_major_length", "axis_minor_length"))
    elong = {l: (a / b if b > 0 else 9.9) for l, a, b in
             zip(props["label"], props["axis_major_length"], props["axis_minor_length"])}
    rdist = ndi.minimum(road_dist, lab2, k) if road_dist is not None else None
    gdist = ndi.minimum(geo_dist, lab2, k) if geo_dist is not None else None

    cands = []
    for i, lk in enumerate(k):
        src = int(klabels[i])
        r, c = minpos[i]
        x, y = transform * (c + 0.5, r + 0.5)
        cy, cx = com[i]
        xc, yc = transform * (cx + 0.5, cy + 0.5)
        a = float(area_m2[src - 1])
        d = float(maxd[src - 1])
        diam = math.sqrt(4 * a / math.pi)
        cands.append(dict(
            label=int(lk), src_label=src, parent=int(parent_lut[src]), nsib=int(nsib_lut[src]),
            row=int(r), col=int(c), x=float(x), y=float(y), xc=float(xc), yc=float(yc),
            depth_m=d, mean_depth_m=float(meand[i]), area_m2=a, diam_m=float(diam),
            volume_m3=float(vol[i]), ratio=float(d / diam) if diam > 0 else 0.0,
            elong=float(elong.get(lk, 1.0)),
            min_z_m=float(dem[r, c]), spill_z_m=float(spill_lut[src]),
            catch_ha=float(maxacc[i]), on_channel=bool(onchan[i]),
            near_nodata=bool(near_nodata[i]), touches_edge=bool(edge[src - 1]),
            road_m=float(rdist[i]) if rdist is not None else None,
            geo_m=float(gdist[i]) if gdist is not None else None,
        ))
    return lab, lab2, cands, depth


def thin_dam_test(dem, valid, lab2, cands, cell, args):
    """How much of each pit's depth survives when every ridge narrower than --dam-width is shaved off.

    A fallen tree across a ravine is often kept as ground in a bare-earth model, and so are fence lines,
    thin berms, and the 'bridges' interpolation can throw across a narrow channel under canopy. Each one
    dams the wash and leaves a closed hollow behind it, measured from the top of the barrier. A grey-scale
    opening with a flat k x k window removes raised features narrower than k cells without touching pits
    (opening only ever lowers ground, and keeps local minima). Each pit is refilled inside a window around
    its footprint, with the window's edge free to drain: a hollow held by a thin barrier drains through
    the gap, a collapse pit keeps nearly all its depth, and a pit separated from the next one by a narrow
    bridge stays closed because the two merge into one hole. Sets depth_open_m, dam_frac and thin_dam."""
    for c in cands:
        c["depth_open_m"], c["dam_frac"], c["thin_dam"] = None, None, False
    if not cands or args.dam_width <= 0:
        return
    k = int(round(args.dam_width / cell)) + 1
    if k % 2 == 0:
        k += 1
    if cell > args.dam_width / 1.5:
        log(f"thin-dam test skipped: {cell:.1f} m cells are too coarse to see a {args.dam_width:g} m ridge")
        return
    objs = ndi.find_objects(lab2)
    H, W = dem.shape
    pad = max(10, 3 * k)
    n_flag = 0
    for c in cands:
        o = objs[c["label"] - 1] if c["label"] - 1 < len(objs) else None
        if o is None:
            continue
        r0, r1 = max(0, o[0].start - pad), min(H, o[0].stop + pad)
        c0, c1 = max(0, o[1].start - pad), min(W, o[1].stop + pad)
        v = valid[r0:r1, c0:c1]
        z = np.where(v, dem[r0:r1, c0:c1], np.nan)
        if not v.any():
            continue
        zhi = np.where(v, z, float(np.nanmax(z)))   # nodata high for the opening, restored after
        opened = ndi.grey_opening(zhi, size=(k, k))
        refill = fill_depressions(np.where(v, opened, np.nan), v)
        rr, cc = c["row"] - r0, c["col"] - c0
        d_open = max(0.0, float(refill[rr, cc] - opened[rr, cc]))
        c["depth_open_m"] = d_open
        c["dam_frac"] = d_open / c["depth_m"] if c["depth_m"] > 0 else 1.0
        if c["dam_frac"] < 0.5:
            c["thin_dam"] = True
            n_flag += 1
    log(f"thin-dam test ({args.dam_width:g} m): {n_flag} of {len(cands)} pits lose more than half their depth "
        f"when ridges thinner than that are removed")


def trace_chains(cands, lab2, rec, dem, cell, transform, args):
    """Link pits that lie on the same drainage line. From each eligible pit the D8 flow path is followed
    downstream (on the conditioned surface, so it leaves the pit at its spill point and runs down the wash)
    until it enters another eligible pit within --chain-dist meters of wash. Each pit therefore has at most one
    downstream neighbour; the links form trees, one per pipe system. For each tree the main stem is the longest
    upstream-to-root path: its along-wash length, straight distance, and floor-to-floor / rim-to-rim elevation
    drop describe the wash between the pits, not the cave under it. Returns a list of chain dicts and fills chain_* on cands.
    Eligibility: ratio >= --chain-ratio, catchment <= --max-catch-ha, not within --road-buffer of a road."""
    for c in cands:
        c["chain_id"], c["chain_n"], c["chain_len_m"] = 0, 1, 0.0
        c["eligible"] = (c["ratio"] >= args.chain_ratio and c["catch_ha"] <= args.max_catch_ha
                         and not (c["road_m"] is not None and c["road_m"] <= args.road_buffer)
                         and not c.get("thin_dam", False))
    elig = [c for c in cands if c["eligible"]]
    if len(elig) < 2:
        return []
    H, W = lab2.shape
    by_label = {c["label"]: c for c in elig}
    # nearest-footprint lookup within a tolerance, so a path that skirts a pit by a few cells still counts as a hit
    tol = max(1, int(round(10.0 / cell)))
    dist_c, (ri, ci) = ndi.distance_transform_edt(lab2 == 0, return_indices=True)
    near = np.where(dist_c <= tol, lab2[ri, ci], 0).astype(np.int32)
    del dist_c, ri, ci
    recf = rec.ravel()
    SQ2 = math.sqrt(2.0)
    links = {}
    for c in elig:
        lab = c["label"]
        cur = c["row"] * W + c["col"]
        dist = 0.0
        path = [(c["row"], c["col"])]
        guard = int(args.chain_dist / cell) * 4 + 100
        while guard > 0:
            guard -= 1
            nxt = int(recf[cur])
            if nxt < 0:
                break
            pr, pc = divmod(cur, W)
            nr, nc = divmod(nxt, W)
            dist += cell * (SQ2 if (nr != pr and nc != pc) else 1.0)
            path.append((nr, nc))
            cur = nxt
            h = int(near[nr, nc])
            if h and h != lab and h in by_label:
                links[lab] = (h, dist, path)
                break
            if dist > args.chain_dist:
                break
    if not links:
        return []
    # components (forest): union-find over links
    parent = {c["label"]: c["label"] for c in elig}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, (b, _, _) in links.items():
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    groups = {}
    for c in elig:
        groups.setdefault(find(c["label"]), []).append(c["label"])

    chains = []
    cid = 0
    for members in groups.values():
        if len(members) < 2:
            continue
        cid += 1
        # distance from each member down to the root (the member with no downstream link inside the tree)
        def to_root(l):
            d, p, seen = 0.0, l, set()
            while p in links and p not in seen:
                seen.add(p)
                nxt, dd, _ = links[p]
                d += dd
                p = nxt
            return d, p
        dists = {l: to_root(l) for l in members}
        root = next(iter(dists.values()))[1]
        top = max(members, key=lambda l: dists[l][0])
        stem_len = dists[top][0]
        # stem path cells, top -> root. Each link's path starts on its pit's floor and stops where it reaches
        # the next pit's edge; the next path starts on that pit's floor. Keep every floor cell, so the line
        # (and any profile drawn along it) goes down into each pit instead of skipping across it.
        cells, p, seen, floors = [], top, set(), set()
        while p in links and p not in seen:
            seen.add(p)
            nxt, _, path = links[p]
            floors.add(path[0])
            cells.extend(path if not cells or cells[-1] != path[0] else path[1:])
            p = nxt
        ct, cr = by_label[top], by_label[root]
        if cells and cells[-1] != (cr["row"], cr["col"]):
            cells.append((cr["row"], cr["col"]))
        floors.add((cr["row"], cr["col"]))
        floor_drop = float(dem[ct["row"], ct["col"]] - dem[cr["row"], cr["col"]])
        rim_drop = float(ct["spill_z_m"] - cr["spill_z_m"])
        straight = math.hypot(ct["x"] - cr["x"], ct["y"] - cr["y"]) * (cell / abs(transform.a))
        chains.append(dict(
            chain_id=cid, n=len(members), members=members, top=top, root=root,
            stem_m=stem_len, straight_m=straight, floor_drop_m=floor_drop, rim_drop_m=rim_drop,
            grade_pct=(100.0 * rim_drop / stem_len) if stem_len > 0 else 0.0,
            cells=[q for i, q in enumerate(cells) if i % max(1, len(cells) // 400) == 0 or q in floors
                   or i == len(cells) - 1] if cells else [],
        ))
        for l in members:
            by_label[l]["chain_id"], by_label[l]["chain_n"], by_label[l]["chain_len_m"] = cid, len(members), float(stem_len)
    chains.sort(key=lambda ch: -ch["stem_m"])
    for i, ch in enumerate(chains, 1):
        ch["chain_id"] = i
        for l in ch["members"]:
            by_label[l]["chain_id"] = i
    return chains


def write_chains(chains, cands, crs, transform, outdir):
    """chains.csv, chains.kml, chains.geojson: one line per pipe system along the wash from its top pit to its root."""
    from rasterio.crs import CRS
    from rasterio.warp import transform as rtransform
    wgs = CRS.from_epsg(4326)
    by_label = {c["label"]: c for c in cands}
    with open(os.path.join(outdir, "chains.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["chain_id", "n_pits", "ranks", "top_rank", "bottom_rank", "stem_ft", "stem_m", "straight_ft",
                     "floor_drop_ft", "rim_drop_ft", "grade_pct", "top_lat", "top_lon", "bottom_lat", "bottom_lon"])
        for ch in chains:
            ranks = sorted(by_label[l]["rank"] for l in ch["members"])
            t, r = by_label[ch["top"]], by_label[ch["root"]]
            wr.writerow([ch["chain_id"], ch["n"], " ".join(str(x) for x in ranks), t["rank"], r["rank"],
                         round(ch["stem_m"] * M2FT, 0), round(ch["stem_m"], 1), round(ch["straight_m"] * M2FT, 0),
                         round(ch["floor_drop_m"] * M2FT, 1), round(ch["rim_drop_m"] * M2FT, 1),
                         round(ch["grade_pct"], 1), t["lat"], t["lon"], r["lat"], r["lon"]])
    feats = []
    kml = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>sink_scout chains</name>',
           '<Style id="chain"><LineStyle><color>ff00d7ff</color><width>3</width></LineStyle></Style>']
    for ch in chains:
        if not ch["cells"]:
            continue
        xs = [transform.c + (cc + 0.5) * transform.a for _, cc in ch["cells"]]
        ys = [transform.f + (rr + 0.5) * transform.e for rr, _ in ch["cells"]]
        lons, lats = rtransform(crs, wgs, xs, ys)
        coords = [[lo, la] for lo, la in zip(lons, lats)]
        ranks = sorted(by_label[l]["rank"] for l in ch["members"])
        name = (f"Chain {ch['chain_id']}: {ch['n']} pits, {ch['stem_m'] * M2FT:,.0f} ft along the wash, "
                f"{ch['floor_drop_m'] * M2FT:.0f} ft drop")
        desc = (f"pits (ranks): {' '.join(str(x) for x in ranks)}<br/>stem {ch['stem_m'] * M2FT:,.0f} ft "
                f"({ch['stem_m']:,.0f} m), straight {ch['straight_m'] * M2FT:,.0f} ft<br/>"
                f"floor-to-floor drop {ch['floor_drop_m'] * M2FT:.0f} ft, rim-to-rim {ch['rim_drop_m'] * M2FT:.0f} ft, "
                f"wash grade {ch['grade_pct']:.1f}%<br/>preliminary: wash length between the top and bottom pits, "
                f"not surveyed cave length")
        feats.append({"type": "Feature", "properties": {"chain_id": ch["chain_id"], "n_pits": ch["n"],
                      "stem_ft": round(ch["stem_m"] * M2FT), "floor_drop_ft": round(ch["floor_drop_m"] * M2FT, 1),
                      "rim_drop_ft": round(ch["rim_drop_m"] * M2FT, 1), "grade_pct": round(ch["grade_pct"], 1),
                      "ranks": ranks}, "geometry": {"type": "LineString", "coordinates": coords}})
        kml.append(f'<Placemark><name>{name}</name><styleUrl>#chain</styleUrl><description><![CDATA[{desc}]]>'
                   f'</description><LineString><tessellate>1</tessellate><coordinates>'
                   + " ".join(f"{lo:.7f},{la:.7f},0" for lo, la in coords) + '</coordinates></LineString></Placemark>')
    kml.append("</Document></kml>")
    with open(os.path.join(outdir, "chains.kml"), "w") as f:
        f.write("\n".join(kml))
    with open(os.path.join(outdir, "chains.geojson"), "w") as f:
        json.dump({"type": "FeatureCollection", "features": feats}, f)


def write_refine_commands(cands, chains, args):
    """After a coarse (--agg) sweep, write one full-resolution --bbox command per chain and per top
    unchained candidate, so the fields the sweep found can be rerun at 1 m without working out boxes by hand.
    The margin matters: flow accumulation only sees the box, so catchment, the on-channel test, the creek
    flag and chain tracing are all box-limited. 300 m is enough to re-measure a pit; use 1000 m or more
    when the verdict depends on the wash above it."""
    top_n = args.refine_top
    margin_m = args.refine_margin
    by_label = {c["label"]: c for c in cands}
    lat0 = cands[0]["lat"]
    dlat = margin_m / 111320.0
    dlon = margin_m / (111320.0 * math.cos(math.radians(lat0)))
    base = [sys.executable.split(os.sep)[-1] if sys.executable else "python", os.path.basename(sys.argv[0])] + list(args.dem)
    if len(args.dem) > 1:
        base.append("--mosaic")
    if args.crs:
        base += ["--crs", args.crs]
    if args.z_units == "ft" and not args.crs:
        base += ["--z-units", "ft"]
    base += ["--min-ratio", f"{args.min_ratio:g}", "--thumbs", "20"]
    if args.osm_roads:
        base.append("--osm-roads")
    if args.roads:
        base += ["--roads", args.roads]
    lines = []
    seen = set()
    for ch in chains:
        lats = [by_label[l]["lat"] for l in ch["members"]]
        lons = [by_label[l]["lon"] for l in ch["members"]]
        seen.update(ch["members"])
        box = f"{min(lats) - dlat:.5f},{min(lons) - dlon:.5f},{max(lats) + dlat:.5f},{max(lons) + dlon:.5f}"
        lines.append(" ".join(base + ["--bbox", box, "--out", f"{args.out}_refine_chain{ch['chain_id']}"]))
    n = 0
    for c in cands:
        if c["label"] in seen or n >= top_n:
            continue
        n += 1
        box = f"{c['lat'] - dlat:.5f},{c['lon'] - dlon:.5f},{c['lat'] + dlat:.5f},{c['lon'] + dlon:.5f}"
        lines.append(" ".join(base + ["--bbox", box, "--check", f"{c['lat']:.5f},{c['lon']:.5f}",
                                      "--out", f"{args.out}_refine_rank{c['rank']}"]))
    with open(os.path.join(args.out, "refine_commands.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return len(lines)


def print_chains(chains, cands, n_show=12):
    if not chains:
        return
    by_label = {c["label"]: c for c in cands}
    print(f"\n{len(chains)} pit chain(s) along a shared drainage (preliminary lengths = wash distance top pit to "
          f"bottom pit, not surveyed length):")
    print(f"{'chain':>5} {'pits':>4} {'stem_ft':>8} {'floor_drop':>10} {'rim_drop':>8} {'grade':>6}  ranks")
    for ch in chains[:n_show]:
        ranks = sorted(by_label[l]["rank"] for l in ch["members"])
        print(f"{ch['chain_id']:>5} {ch['n']:>4} {ch['stem_m'] * M2FT:>8,.0f} {ch['floor_drop_m'] * M2FT:>10.0f} "
              f"{ch['rim_drop_m'] * M2FT:>8.0f} {ch['grade_pct']:>5.1f}%  {' '.join(str(x) for x in ranks)}")
    print("  floor_drop = top pit floor to bottom pit floor (proxy for surveyed depth); rim_drop = rim to rim "
          "(wash grade)")
    if len(chains) > n_show:
        print(f"... {len(chains) - n_show} more in chains.csv")


def score_candidates(cands, args, have_roads, have_geo):
    for c in cands:
        s = 0.0
        s += min(c["depth_m"] / 4.0, 1.5)                       # 4 m deep -> 1
        s += min(c["ratio"] / 0.3, 1.5)                          # depth/width 0.3 -> 1
        s += float(np.clip(math.log10(c["catch_ha"] / args.channel_ha + 1.0), 0, 1.5))
        s += 0.75 * min(c["chain_n"] - 1, 2)
        s -= 0.5 * float(np.clip((c["elong"] - 2.0) / 2.0, 0, 1))
        mult = 1.0
        notes = []
        if not c["on_channel"]:
            mult *= 0.5
            notes.append("off-channel")
        if c["ratio"] < args.min_ratio:
            mult *= 0.4
            notes.append("shallow/wide")
        if c["catch_ha"] > args.max_catch_ha:
            mult *= 0.5
            notes.append("creek")
        if have_roads and c["road_m"] is not None and c["road_m"] <= args.road_buffer:
            mult *= 0.15
            notes.append("road")
        if c.get("thin_dam"):
            mult *= 0.4
            notes.append("thin-dam")
        if have_geo and c["geo_m"] is not None and c["geo_m"] > args.geo_buffer:
            mult *= 0.6
            notes.append("off-geology")
        if c["near_nodata"]:
            mult *= 0.7
            notes.append("nodata-edge")
        if c.get("nsib", 1) > 1:
            notes.append(f"nested({c['nsib']})")
        c["score"] = round(max(s, 0.05) * mult, 3)
        c["flags"] = ",".join(notes)
    cands.sort(key=lambda c: -c["score"])
    for i, c in enumerate(cands, 1):
        c["rank"] = i


# ----------------------------------------------------------------------------- outputs

FIELDS = [  # (shapefile name, csv name, type, decimals)
    ("rank", "rank", "N", 0), ("score", "score", "F", 3),
    ("depth_m", "depth_m", "F", 2), ("depth_ft", "depth_ft", "F", 1),
    ("diam_m", "diam_m", "F", 1), ("diam_ft", "diam_ft", "F", 1),
    ("area_m2", "area_m2", "F", 1), ("vol_m3", "volume_m3", "F", 1),
    ("ratio", "ratio", "F", 3), ("elong", "elong", "F", 2),
    ("catch_ha", "catch_ha", "F", 2), ("catch_ac", "catch_ac", "F", 1),
    ("on_chan", "on_channel", "N", 0), ("road_m", "road_m", "F", 1), ("geo_m", "geo_m", "F", 1),
    ("chain_id", "chain_id", "N", 0), ("chain_n", "chain_n", "N", 0), ("chain_m", "chain_len_m", "F", 1),
    ("min_z_m", "min_z_m", "F", 2), ("spill_z_m", "spill_z_m", "F", 2),
    ("dep_opn_ft", "depth_open_ft", "F", 1), ("dam_frac", "dam_frac", "F", 2),
    ("lat", "lat", "F", 7), ("lon", "lon", "F", 7), ("x", "x", "F", 2), ("y", "y", "F", 2),
    ("flags", "flags", "C", 0), ("label", "label", "N", 0),
    ("parent", "parent", "N", 0), ("nsib", "nsib", "N", 0),
]


def finish_fields(cands, crs):
    from rasterio.crs import CRS
    from rasterio.warp import transform as rtransform
    if not cands:
        return
    lons, lats = rtransform(crs, CRS.from_epsg(4326), [c["x"] for c in cands], [c["y"] for c in cands])
    for c, lo, la in zip(cands, lons, lats):
        c["lat"], c["lon"] = round(la, 7), round(lo, 7)
        c["depth_ft"] = c["depth_m"] * M2FT
        c["diam_ft"] = c["diam_m"] * M2FT
        c["catch_ac"] = c["catch_ha"] * HA2AC
        c["depth_open_ft"] = c["depth_open_m"] * M2FT if c.get("depth_open_m") is not None else None


def ring_area(ring):
    a = 0.0
    for (x1, y1), (x2, y2) in zip(ring, ring[1:] + ring[:1]):
        a += x1 * y2 - x2 * y1
    return 0.5 * a


def footprint_polygons(lab2, transform):
    """{label: [ring, ring, ...]} in DEM CRS, exterior clockwise, holes counter-clockwise (shapefile convention)."""
    from rasterio import features
    out = {}
    for geom, val in features.shapes(lab2.astype(np.int32), mask=lab2 > 0, transform=transform, connectivity=8):
        rings = []
        for i, ring in enumerate(geom["coordinates"]):
            pts = [tuple(p) for p in ring]
            a = ring_area(pts)
            if (i == 0 and a > 0) or (i > 0 and a < 0):
                pts = pts[::-1]
            rings.append(pts)
        out.setdefault(int(val), []).extend(rings)
    return out


def write_prj(path, crs):
    try:
        from rasterio.enums import WktVersion
        wkt = crs.to_wkt(version=WktVersion.WKT1_ESRI)
    except Exception:
        wkt = crs.to_wkt()
    with open(path, "w") as f:
        f.write(wkt)


def record_for(c):
    rec = []
    for shp_name, key, typ, dec in FIELDS:
        v = c.get(key)
        if typ == "N":
            rec.append(int(v) if v is not None else 0)
        elif typ == "F":
            rec.append(round(float(v), dec) if v is not None else -1.0)
        else:
            rec.append(str(v) if v is not None else "")
    return rec


def write_shapefiles(cands, polys, crs, outdir):
    shapefile = need("shapefile", "pyshp")
    for kind in ("candidates", "depressions"):
        path = os.path.join(outdir, kind)
        w = shapefile.Writer(path, shapeType=shapefile.POINT if kind == "candidates" else shapefile.POLYGON)
        for shp_name, key, typ, dec in FIELDS:
            if typ == "C":
                w.field(shp_name, "C", 40)
            elif typ == "N":
                w.field(shp_name, "N", 10, 0)
            else:
                w.field(shp_name, "F", 16, dec)
        for c in cands:
            if kind == "candidates":
                w.point(c["x"], c["y"])
            else:
                rings = polys.get(c["label"])
                if not rings:
                    continue
                w.poly([[list(p) for p in r] for r in rings])
            w.record(*record_for(c))
        w.close()
        write_prj(path + ".prj", crs)


def write_csv(cands, outdir):
    with open(os.path.join(outdir, "candidates.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow([key for _, key, _, _ in FIELDS])
        for c in cands:
            wr.writerow(record_for(c))
    # pit-shaped holes with no drainage into them: cellar holes, prospects, cisterns, talus voids, artifacts.
    # Not pipes, but a by-product list worth having for anyone hunting old sites.
    poi = sorted([c for c in cands if not c["on_channel"] and c["ratio"] >= 0.25], key=lambda c: -c["depth_m"])
    with open(os.path.join(outdir, "offchannel_pits.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow([key for _, key, _, _ in FIELDS])
        for c in poi:
            wr.writerow(record_for(c))
    return len(poi)


def write_geojson_kml(cands, polys, crs, outdir):
    from rasterio.crs import CRS
    from rasterio.warp import transform as rtransform
    wgs = CRS.from_epsg(4326)
    feats = []
    kml = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>sink_scout candidates</name>',
           '<Style id="top"><IconStyle><color>ff0000ff</color><scale>1.3</scale></IconStyle>'
           '<LineStyle><color>ff0000ff</color><width>2.5</width></LineStyle><PolyStyle><fill>0</fill></PolyStyle></Style>',
           '<Style id="mid"><IconStyle><color>ff00a5ff</color><scale>1.0</scale></IconStyle>'
           '<LineStyle><color>ff00a5ff</color><width>2</width></LineStyle><PolyStyle><fill>0</fill></PolyStyle></Style>',
           '<Style id="low"><IconStyle><color>ff00ffff</color><scale>0.8</scale></IconStyle>'
           '<LineStyle><color>ff00ffff</color><width>1.5</width></LineStyle><PolyStyle><fill>0</fill></PolyStyle></Style>']
    for c in cands:
        props = {key: c.get(key) for _, key, _, _ in FIELDS}
        feats.append({"type": "Feature", "properties": props,
                      "geometry": {"type": "Point", "coordinates": [c["lon"], c["lat"]]}})
        rings_ll = []
        for ring in polys.get(c["label"], []):
            xs, ys = rtransform(crs, wgs, [p[0] for p in ring], [p[1] for p in ring])
            rings_ll.append([[x, y] for x, y in zip(xs, ys)])
        if rings_ll:
            feats.append({"type": "Feature", "properties": dict(props, geom="footprint"),
                          "geometry": {"type": "Polygon", "coordinates": rings_ll}})
        style = "top" if c["rank"] <= 10 else ("mid" if c["rank"] <= 30 else "low")
        desc = (f"score {c['score']}<br/>depth {c['depth_ft']:.0f} ft ({c['depth_m']:.1f} m)<br/>"
                f"diameter {c['diam_ft']:.0f} ft<br/>depth/width {c['ratio']:.2f}<br/>"
                f"catchment {c['catch_ac']:.0f} ac<br/>chain {c['chain_n']} sinks<br/>flags {c['flags'] or 'none'}")
        kml.append(f'<Placemark><name>#{c["rank"]}</name><styleUrl>#{style}</styleUrl>'
                   f'<description><![CDATA[{desc}]]></description>'
                   f'<Point><coordinates>{c["lon"]:.7f},{c["lat"]:.7f},0</coordinates></Point></Placemark>')
        for ring in rings_ll:
            coords = " ".join(f"{x:.7f},{y:.7f},0" for x, y in ring)
            kml.append(f'<Placemark><name>#{c["rank"]} footprint</name><styleUrl>#{style}</styleUrl>'
                       f'<Polygon><outerBoundaryIs><LinearRing><coordinates>{coords}</coordinates>'
                       f'</LinearRing></outerBoundaryIs></Polygon></Placemark>')
    kml.append("</Document></kml>")
    with open(os.path.join(outdir, "candidates.geojson"), "w") as f:
        json.dump({"type": "FeatureCollection", "features": feats}, f)
    with open(os.path.join(outdir, "candidates.kml"), "w") as f:
        f.write("\n".join(kml))


def write_rasters(outdir, transform, crs, depth, acc_ha, lab2, dem=None, valid=None):
    import rasterio
    prof = dict(driver="GTiff", height=depth.shape[0], width=depth.shape[1], count=1,
                transform=transform, crs=crs, compress="deflate", tiled=True)
    if dem is not None:
        # the ground as the program saw it (after --agg / --bbox / unit conversion), in meters:
        # the one file from which any profile, hillshade or cross-section can be rebuilt
        z = np.where(valid, dem, -9999.0).astype(np.float32) if valid is not None else dem.astype(np.float32)
        with rasterio.open(os.path.join(outdir, "dem.tif"), "w", dtype="float32", nodata=-9999.0, **prof) as d:
            d.write(z, 1)
    with rasterio.open(os.path.join(outdir, "fill_depth.tif"), "w", dtype="float32", nodata=-1, **prof) as d:
        d.write(depth.astype(np.float32), 1)
    with rasterio.open(os.path.join(outdir, "flowacc_ha.tif"), "w", dtype="float32", nodata=-1, **prof) as d:
        d.write(acc_ha.astype(np.float32), 1)
    with rasterio.open(os.path.join(outdir, "labels.tif"), "w", dtype="int32", nodata=0, **prof) as d:
        d.write(lab2.astype(np.int32), 1)


def hillshade(dem, cell):
    from matplotlib.colors import LightSource
    ls = LightSource(azdeg=315, altdeg=45)
    z = np.array(dem, dtype=np.float64)
    nanmask = ~np.isfinite(z)
    if nanmask.any():
        z[nanmask] = np.nanmedian(z)
    hs = ls.hillshade(z, vert_exag=1.5, dx=cell, dy=cell)
    hs[nanmask] = np.nan
    return hs


def write_overview(outdir, dem, valid, acc_ha, lab2, cands, transform, cell, args, chains=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from rasterio.transform import array_bounds
    H, W = dem.shape
    s = max(1, int(math.ceil(max(H, W) / 2500)))
    hs = hillshade(dem[::s, ::s], cell * s)
    chan = (acc_ha >= args.channel_ha)[::s, ::s]
    foot = (lab2 > 0)[::s, ::s]
    w, so, e, n = array_bounds(H, W, transform)
    ext = (w, e, so, n)
    fig_w = 12
    fig, ax = plt.subplots(figsize=(fig_w, fig_w * H / W), dpi=130)
    ax.imshow(hs, cmap="gray", extent=ext, vmin=0, vmax=1, interpolation="nearest")
    chan_rgba = np.zeros(chan.shape + (4,), np.float32)
    chan_rgba[chan] = (0.15, 0.45, 0.95, 0.8)
    ax.imshow(chan_rgba, extent=ext, interpolation="nearest")
    foot_rgba = np.zeros(foot.shape + (4,), np.float32)
    foot_rgba[foot] = (1.0, 0.15, 0.1, 0.55)
    ax.imshow(foot_rgba, extent=ext, interpolation="nearest")
    for ch in (chains or []):
        if ch["cells"]:
            xs = [transform.c + (cc + 0.5) * transform.a for _, cc in ch["cells"]]
            ys = [transform.f + (rr + 0.5) * transform.e for rr, _ in ch["cells"]]
            ax.plot(xs, ys, color="orange", lw=1.8, alpha=0.9, zorder=4)
    if cands:
        sc = np.array([c["score"] for c in cands])
        ax.scatter([c["x"] for c in cands], [c["y"] for c in cands], s=28, c=sc, cmap="plasma",
                   edgecolors="white", linewidths=0.6, zorder=5)
        for c in cands[: min(len(cands), 40)]:
            ax.annotate(str(c["rank"]), (c["x"], c["y"]), xytext=(5, 4), textcoords="offset points",
                        fontsize=8, color="white", weight="bold", zorder=6,
                        path_effects=[pe.withStroke(linewidth=2, foreground="black")])
    ax.set_title(f"sink_scout — {len(cands)} candidates (red = footprints, blue = channels ≥ {args.channel_ha:g} ha, "
                 f"depth ≥ {args.min_depth:g} m)", fontsize=10)
    if args.no_coords:
        ax.set_xticks([])
        ax.set_yticks([])
        m_per_unit = cell / abs(transform.a)
        span_m = (e - w) * m_per_unit
        bar_m = 500.0 if span_m > 3000 else (200.0 if span_m > 1200 else 50.0)
        bar = bar_m / m_per_unit
        x0, y0 = w + 0.03 * (e - w), so + 0.04 * (n - so)
        ax.plot([x0, x0 + bar], [y0, y0], color="black", lw=5, solid_capstyle="butt", zorder=7)
        ax.plot([x0, x0 + bar], [y0, y0], color="white", lw=2.5, solid_capstyle="butt", zorder=8)
        ax.text(x0, y0 + 0.012 * (n - so), f"{bar_m:g} m   north up", fontsize=8, color="white", zorder=8,
                bbox=dict(facecolor="black", alpha=0.5, pad=1.5, lw=0))
    else:
        ax.ticklabel_format(style="plain", useOffset=False)
        ax.tick_params(labelsize=7)
        ax.set_xlabel("x (DEM CRS)", fontsize=8)
        ax.set_ylabel("y (DEM CRS)", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "overview.png"))
    plt.close(fig)


def write_thumbs(outdir, dem, valid, depth, acc_ha, lab2, cands, transform, cell, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tdir = os.path.join(outdir, "thumbs")
    os.makedirs(tdir, exist_ok=True)
    H, W = dem.shape
    half = int(round(args.thumb_size / 2 / cell))
    for c in cands[: args.thumbs]:
        r0, r1 = max(0, c["row"] - half), min(H, c["row"] + half)
        c0, c1 = max(0, c["col"] - half), min(W, c["col"] + half)
        sub = dem[r0:r1, c0:c1]
        hs = hillshade(sub, cell)
        chan = acc_ha[r0:r1, c0:c1] >= args.channel_ha
        own = lab2[r0:r1, c0:c1] == c["label"]
        others = (lab2[r0:r1, c0:c1] > 0) & ~own
        fig, ax = plt.subplots(figsize=(6, 6), dpi=110)
        ax.imshow(hs, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
        ax.imshow(np.ma.masked_where(~chan, chan), cmap="Blues", vmin=0, vmax=1.3, alpha=0.7, interpolation="nearest")
        if others.any():
            ax.contour(others.astype(float), levels=[0.5], colors="orange", linewidths=1.2)
        ax.contour(own.astype(float), levels=[0.5], colors="red", linewidths=1.6)
        ax.plot(c["col"] - c0, c["row"] - r0, "x", color="yellow", ms=10, mew=2)
        where = "" if args.no_coords else f"{c['lat']:.5f}, {c['lon']:.5f}   "
        ax.set_title(f"#{c['rank']}  score {c['score']}   depth {c['depth_ft']:.0f} ft   diam {c['diam_ft']:.0f} ft\n"
                     f"d/w {c['ratio']:.2f}   catch {c['catch_ac']:.0f} ac   chain {c['chain_n']}   "
                     f"{where}{c['flags']}", fontsize=8.5)
        ax.set_xticks([])
        ax.set_yticks([])
        bar = 50.0 / cell
        ax.plot([8, 8 + bar], [sub.shape[0] - 10, sub.shape[0] - 10], color="black", lw=4, solid_capstyle="butt")
        ax.plot([8, 8 + bar], [sub.shape[0] - 10, sub.shape[0] - 10], color="white", lw=2, solid_capstyle="butt")
        ax.text(8, sub.shape[0] - 16, f"50 m   ({args.thumb_size:.0f} m window, north up)", fontsize=7,
                color="white", va="bottom", bbox=dict(facecolor="black", alpha=0.5, pad=1.5, lw=0))
        fig.tight_layout()
        fig.savefig(os.path.join(tdir, f"rank_{c['rank']:03d}.png"))
        plt.close(fig)


def print_report(cands, n_show=15):
    if not cands:
        print("\nNo candidates passed the filters.")
        return
    print(f"\n{'rank':>4} {'score':>6} {'depth_ft':>8} {'diam_ft':>7} {'d/w':>5} {'catch_ac':>8} "
          f"{'chain':>5} {'lat':>10} {'lon':>11}  flags")
    for c in cands[:n_show]:
        print(f"{c['rank']:>4} {c['score']:>6.2f} {c['depth_ft']:>8.0f} {c['diam_ft']:>7.0f} {c['ratio']:>5.2f} "
              f"{c['catch_ac']:>8.0f} {c['chain_n']:>5} {c['lat']:>10.5f} {c['lon']:>11.5f}  {c['flags']}")
    if len(cands) > n_show:
        print(f"... {len(cands) - n_show} more in candidates.csv")


def check_point(text, cands, lab2, lab_all, depth, transform, crs, cell, ufac=1.0):
    from rasterio.crs import CRS
    from rasterio.warp import transform as rtransform
    try:
        lat, lon = [float(v) for v in text.split(",")]
    except ValueError:
        sys.exit("--check wants lat,lon in decimal degrees, e.g. --check 39.08,-108.54")
    xs, ys = rtransform(CRS.from_epsg(4326), crs, [lon], [lat])
    col, row = ~transform * (xs[0], ys[0])
    row, col = int(row), int(col)
    print(f"\n--check {lat:.5f},{lon:.5f} -> row {row}, col {col}")
    H, W = lab2.shape
    if not (0 <= row < H and 0 <= col < W):
        print("  point is outside this DEM")
        return
    win = 2
    r0, r1, c0, c1 = max(0, row - win), min(H, row + win + 1), max(0, col - win), min(W, col + win + 1)
    hit = lab2[r0:r1, c0:c1]
    hit = hit[hit > 0]
    if hit.size:
        lab = int(np.bincount(hit).argmax())
        c = next(x for x in cands if x["label"] == lab)
        print(f"  HIT: candidate rank #{c['rank']} (score {c['score']}) — depth {c['depth_ft']:.0f} ft, "
              f"diam {c['diam_ft']:.0f} ft, d/w {c['ratio']:.2f}, catch {c['catch_ac']:.0f} ac, chain of {c['chain_n']}"
              f"{', flags: ' + c['flags'] if c['flags'] else ''}")
    else:
        hit_all = lab_all[r0:r1, c0:c1]
        hit_all = hit_all[hit_all > 0]
        if hit_all.size:
            l = int(np.bincount(hit_all).argmax())
            m = lab_all == l
            print(f"  the fill found a closed basin here (max fill depth {depth[m].max():.2f} m, "
                  f"area {m.sum() * cell * cell:.0f} m2) but no pit in it survived: either nothing inside has its own "
                  f"depth >= --min-depth, or it failed --min-area/--max-area/--keep-edge")
        else:
            print(f"  no closed depression at that cell (fill depth {depth[row, col]:.2f} m)")
    if cands:
        d = [math.hypot(c["x"] - xs[0], c["y"] - ys[0]) * ufac for c in cands]
        i = int(np.argmin(d))
        print(f"  nearest candidate: rank #{cands[i]['rank']} at {d[i]:.0f} m")


# ----------------------------------------------------------------------------- self test


def make_selftest(outdir):
    """Synthetic 4-ft DEM in the Mesa CRS with planted features; returns (dem_path, roads_path, truth)."""
    import rasterio
    from rasterio.transform import from_origin
    from rasterio.crs import CRS
    from rasterio.warp import transform as rtransform

    rng = np.random.default_rng(7)
    H = W = 900
    cell_ft = 4.0
    cell_m = cell_ft * 0.3048006096
    x0, y0 = 2_270_000.0, 14_300_000.0  # arbitrary UTM12 survey-feet origin
    yy, xx = np.mgrid[0:H, 0:W]
    xm, ym = xx * cell_m, (H - yy) * cell_m
    z = 1600.0 + 0.03 * xm - 0.01 * ym                                   # regional slope
    z += 6.0 * np.sin(xm / 180.0) * np.cos(ym / 230.0)                   # rolling terrain
    # valley: axis runs from upper-right to lower-left along y = 0.45x + 200 (m), draining to lower-left
    axis_y = 0.45 * xm + 200.0
    dvalley = np.abs(ym - axis_y) / math.sqrt(1 + 0.45 ** 2)
    z -= 14.0 * np.exp(-(dvalley / 60.0) ** 2)
    z += 0.08 * rng.standard_normal((H, W))                                # lidar noise

    def pit(xc_m, yc_m, diam_m, depth_m):
        r = np.hypot(xm - xc_m, ym - yc_m)
        prof = depth_m * np.clip(1 - (r / (diam_m / 2)) ** 1.4, 0, None)
        return prof

    truth = []
    # chain of three sinks on the valley axis
    for xc, dia, dep, name in [(520, 13, 12.0, "chain sink A (40 ft)"), (590, 9, 4.5, "chain sink B"),
                               (650, 7, 2.6, "chain sink C")]:
        yc = 0.45 * xc + 200.0
        z -= pit(xc, yc, dia, dep)
        truth.append((name, xc, yc))
    # lone sink far up the valley
    xc = 870; yc = 0.45 * xc + 200
    z -= pit(xc, yc, 10, 5.0); truth.append(("lone valley sink", xc, yc))
    # hillside pit, off any channel
    z -= pit(250, 850, 9, 3.5); truth.append(("hillside pit (off-channel)", 250, 850))
    # stock pond: shallow, wide, on a flat shelf
    z -= pit(820, 180, 70, 3.2); truth.append(("stock pond (shallow/wide)", 820, 180))
    # two sinks just upstream of the road: they sit inside the dammed basin and must be split out of it
    for xc, dia, dep, name in [(380, 10, 6.0, "sink D inside road basin"), (430, 8, 3.0, "sink E inside road basin")]:
        yc = 0.45 * xc + 200.0
        z -= pit(xc, yc, dia, dep)
        truth.append((name, xc, yc))
    # road embankment crossing the valley downstream -> dammed basin
    road_x = 330.0
    road = np.abs(xm - road_x) < 4.0
    z[road] = np.maximum(z[road], (z[road] + 2.5))
    truth.append(("road dam basin", road_x + 8, 0.45 * (road_x + 8) + 200))
    # nodata corner
    z[:120, :120] = np.nan

    crs = CRS.from_user_input(CRS_SHORTCUTS["mesa"])
    transform = from_origin(x0, y0, cell_ft, cell_ft)
    dem_path = os.path.join(outdir, "selftest_dem.tif")
    zft = np.where(np.isfinite(z), z / 0.3048006096, -9999.0).astype(np.float32)
    with rasterio.open(dem_path, "w", driver="GTiff", height=H, width=W, count=1, dtype="float32",
                       crs=crs, transform=transform, nodata=-9999.0) as d:
        d.write(zft, 1)
    # road as a GeoJSON line in lon/lat (tests the reprojection path)
    col = road_x / cell_m
    xs = [x0 + col * cell_ft] * 2
    ys = [y0, y0 - H * cell_ft]
    lons, lats = rtransform(crs, CRS.from_epsg(4326), xs, ys)
    roads_path = os.path.join(outdir, "selftest_roads.geojson")
    with open(roads_path, "w") as f:
        json.dump({"type": "FeatureCollection", "features": [{"type": "Feature", "properties": {"highway": "track"},
                   "geometry": {"type": "LineString", "coordinates": [[lons[0], lats[0]], [lons[1], lats[1]]]}}]}, f)
    # truth points -> lat/lon
    txs = [x0 + (t[1] / cell_m) * cell_ft for t in truth]
    tys = [y0 - (H - t[2] / cell_m) * cell_ft for t in truth]
    tlon, tlat = rtransform(crs, CRS.from_epsg(4326), txs, tys)
    truth_ll = [(t[0], la, lo) for t, la, lo in zip(truth, tlat, tlon)]
    return dem_path, roads_path, truth_ll


# ----------------------------------------------------------------------------- main


def build_parser():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("dem", nargs="*", help="DEM file(s): GeoTIFF, .img, anything GDAL reads")
    p.add_argument("--out", default="sink_scout_out", help="output folder")
    p.add_argument("--crs", help="override DEM CRS: EPSG:26912, a proj string, or 'mesa'")
    p.add_argument("--z-units", choices=["m", "ft"], help="elevation units (default m; ft when --crs mesa)")
    p.add_argument("--nodata", type=float, help="override nodata value")
    p.add_argument("--zmin", type=float, default=-1000.0, help="elevations below this are nodata (default -1000)")
    p.add_argument("--agg", type=int, default=1, help="block-mean downsample factor (e.g. 3 for 1-ft grids)")
    p.add_argument("--mosaic", action="store_true", help="merge several DEMs before processing")
    p.add_argument("--bbox", help="clip to lat1,lon1,lat2,lon2 (two opposite corners) before processing; "
                                  "lets a 1-m tile run at full resolution over a small area")
    # detection
    p.add_argument("--floor", type=float, default=0.3, help="fill depth (m) defining a footprint (default 0.3)")
    p.add_argument("--min-depth", type=float, default=1.5, help="min max-depth (m) to keep (default 1.5)")
    p.add_argument("--max-depth", type=float, default=60.0, help="max depth (m) (default 60)")
    p.add_argument("--min-area", type=float, default=15.0, help="min footprint area m2 (default 15)")
    p.add_argument("--max-area", type=float, default=50000.0, help="max footprint area m2 (default 50000)")
    p.add_argument("--min-ratio", type=float, default=0.05, help="depth/width below this is 'shallow/wide' (default 0.05)")
    p.add_argument("--channel-ha", type=float, default=2.0, help="catchment (ha) that counts as a channel (default 2)")
    p.add_argument("--chain-dist", type=float, default=200.0,
                   help="max wash distance (m) from one pit to the next along the flow path for them to be one chain (default 200)")
    p.add_argument("--chain-ratio", type=float, default=0.1,
                   help="depth/width a candidate needs to count toward a chain (default 0.1; kept below --min-ratio "
                        "so a partly collapsed slot still links to the pit beside it)")
    p.add_argument("--max-catch-ha", type=float, default=300.0,
                   help="catchment (ha) above which a depression is flagged 'creek' and kept out of chains "
                        "(default 300 ha, about 740 acres); pipes form under small washes, dams sit on creeks")
    p.add_argument("--keep-edge", action="store_true", help="keep depressions touching the raster edge")
    p.add_argument("--rim-slope", type=float, default=0.08,
                   help="wall slope (rise/run) below which the flood-fill has left the pit and found its rim "
                        "(default 0.08; 0 disables rim trimming)")
    p.add_argument("--eps", type=float, default=1e-5, help="flat-resolution gradient (m per cell)")
    # context layers
    p.add_argument("--roads", help="road lines/polygons (.shp or GeoJSON) to mask culvert-dam false positives")
    p.add_argument("--osm-roads", action="store_true", help="download roads for the DEM bbox from OpenStreetMap")
    p.add_argument("--road-buffer", type=float, default=30.0, help="meters from a road that flags a candidate (default 30)")
    p.add_argument("--geology", help="geologic map polygons (.shp or GeoJSON)")
    p.add_argument("--geo-field", default="UNIT", help="attribute holding the unit code (default UNIT)")
    p.add_argument("--geo-units", default="", help="comma list of unit codes to favor, e.g. TRc,Jm (default: all polygons)")
    p.add_argument("--geo-buffer", type=float, default=300.0, help="meters from the geology that still counts (default 300)")
    # outputs
    p.add_argument("--thumbs", type=int, default=25, help="close-up PNGs for the top N (default 25, 0 = none)")
    p.add_argument("--thumb-size", type=float, default=300.0, help="close-up window size in meters (default 300)")
    p.add_argument("--no-overview", action="store_true")
    p.add_argument("--no-coords", action="store_true",
                   help="draw the overview and close-ups without coordinates (scale bar instead of axes, no lat/lon "
                        "in titles) so they can be published without giving locations away")
    p.add_argument("--save-rasters", action="store_true", help="write fill_depth / flowacc / labels GeoTIFFs")
    p.add_argument("--dam-width", type=float, default=4.0,
                   help="meters: ridges thinner than this (fallen logs, fences, thin berms, lidar bridges) are "
                        "shaved off and each pit measured again; pits losing over half their depth are flagged "
                        "thin-dam (default 4; 0 turns the test off; skipped on cells coarser than 2.7 m)")
    p.add_argument("--refine-margin", type=float, default=300.0,
                   help="meters of margin around each chain or pit in refine_commands.txt (default 300). "
                        "Catchment, the creek flag and chains are computed inside the box, so use 1000 or "
                        "more when a pit sits on a big wash")
    p.add_argument("--refine-top", type=int, default=15,
                   help="how many unchained top candidates get a refine command (default 15)")
    p.add_argument("--check", help="lat,lon of a known sinkhole; reports its rank or why it was dropped")
    p.add_argument("--selftest", action="store_true", help="build a synthetic DEM with planted sinks and run on it")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    truth = None
    if args.selftest:
        log("building synthetic test DEM (4-ft cells, Mesa CRS, elevations in feet)")
        dem_path, roads_path, truth = make_selftest(args.out)
        args.dem = [dem_path]
        args.crs = args.crs or "mesa"
        args.roads = args.roads or roads_path
        args.thumbs = min(args.thumbs, 8)
    if not args.dem:
        sys.exit("No DEM given. Try: python sink_scout.py --selftest   or   python sink_scout.py --help")
    if args.z_units is None:
        args.z_units = "ft" if (args.crs or "").lower() == "mesa" else "m"
        if args.z_units == "ft":
            log("--crs mesa: assuming elevations in feet (override with --z-units m)")

    t0 = time.time()
    dem, valid, transform, crs, cell, ufac = load_dem(args.dem, args.crs, args.mosaic, args.nodata, args.zmin,
                                                      args.z_units, args.agg, args.bbox)
    cell_ha = cell * cell / 10000.0

    log("filling depressions...")
    filled = fill_depressions(dem, valid)
    log("resolving flats...")
    zeps, stuck = resolve_flats(filled, valid, args.eps)
    if stuck:
        log(f"warning: {stuck} flat cells had no outlet (treated as sinks)")
    log("D8 flow directions...")
    rec, bad = d8_receivers(zeps, valid)
    if bad:
        log(f"warning: {bad} cells without a downslope neighbor after conditioning")
    log("flow accumulation...")
    acc, rounds = flow_accumulation(rec, valid)
    acc_ha = acc * cell_ha
    log(f"routing done in {rounds} rounds; longest channel {acc_ha.max():.1f} ha")
    del zeps

    road_dist = None
    road_geoms = []
    if args.roads:
        geoms, vcrs = read_vector(args.roads)
        road_geoms += [reproject_geom(g, vcrs, crs) for g in geoms]
        log(f"roads: {len(geoms)} features from {args.roads}")
    if args.osm_roads:
        road_geoms += fetch_osm_roads(crs, transform, dem.shape, os.path.join(args.out, "osm_roads.geojson"))
    if road_geoms:
        road_dist = rasterize_distance(road_geoms, dem.shape, transform, cell)
        if road_dist is None:
            log("warning: no road geometry fell inside the DEM")

    geo_dist = None
    if args.geology:
        geoms, vcrs = read_vector(args.geology)
        want = [u.strip().lower() for u in args.geo_units.split(",") if u.strip()]
        if want:
            geoms = [g for g in geoms if str(g.get("properties", {}).get(args.geo_field, "")).strip().lower() in want]
        geoms = [reproject_geom(g, vcrs, crs) for g in geoms]
        log(f"geology: {len(geoms)} polygons" + (f" matching {args.geo_units}" if want else ""))
        geo_dist = rasterize_distance(geoms, dem.shape, transform, cell)
        if geo_dist is None:
            log("warning: no geology polygon fell inside the DEM; geology ignored")

    log("measuring depressions...")
    lab_all, lab2, cands, depth = measure_depressions(dem, valid, filled, acc_ha, transform, cell, args,
                                                      road_dist, geo_dist)
    thin_dam_test(dem, valid, lab2, cands, cell, args)
    chains = trace_chains(cands, lab2, rec, dem, cell, transform, args) if cands else []
    del rec
    score_candidates(cands, args, road_dist is not None, geo_dist is not None)
    finish_fields(cands, crs)

    log("writing outputs...")
    polys = footprint_polygons(lab2, transform) if cands else {}
    write_shapefiles(cands, polys, crs, args.out)
    n_poi = write_csv(cands, args.out)
    write_geojson_kml(cands, polys, crs, args.out)
    if chains:
        write_chains(chains, cands, crs, transform, args.out)
    if args.save_rasters:
        write_rasters(args.out, transform, crs, depth, acc_ha, lab2, dem, valid)
    if not args.no_overview:
        write_overview(args.out, dem, valid, acc_ha, lab2, cands, transform, cell, args, chains)
    if args.thumbs > 0 and cands:
        write_thumbs(args.out, dem, valid, depth, acc_ha, lab2, cands, transform, cell, args)

    print_report(cands)
    print_chains(chains, cands)
    if args.agg > 1 and cands:
        n_ref = write_refine_commands(cands, chains, args)
        print(f"\n{n_ref} full-resolution follow-up command(s) written to refine_commands.txt "
              f"(one --bbox run per chain and per top candidate)")
    if n_poi:
        print(f"\n{n_poi} off-channel pit-shaped candidate(s) (cellar holes, prospects, cisterns, artifacts) "
              f"-> offchannel_pits.csv")
    if args.check:
        check_point(args.check, cands, lab2, lab_all, depth, transform, crs, cell, ufac)
    if truth:
        print("\nself-test: planted features")
        for name, la, lo in truth:
            print(f"  -- {name}")
            check_point(f"{la},{lo}", cands, lab2, lab_all, depth, transform, crs, cell, ufac)
        print("\nexpected: the five valley sinks (A-E) on top as one chain of 5 (D and E pulled out of the road-dammed\n"
              "          basin), lone valley sink next, hillside pit flagged off-channel, pond flagged shallow/wide,\n"
              "          and the dammed basin itself gone (no pit of its own).")
    log(f"done in {time.time() - t0:.1f} s -> {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
