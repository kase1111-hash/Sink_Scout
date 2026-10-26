# sink_scout, county_scan, field_kit

Finds collapse-sinkhole candidates in bare-earth lidar: pits on drainages, linked into chains down their washes, ranked by how much they look like the ground over a piping cave. Written for the book *Hollowed Ground: The Piping Caves of Colorado's Western Slope* (Kase Branham, 2026). The book's Appendix A is the beginner's path. This file is the reference.

**The outputs of these programs are maps of places where the ground is hollow. Do not publish them.** The one exception is `county_density_nocoords.png`, which shows only counts per two-kilometer square. Figures for publication are made with `--no-coords`. Report what you find to the people who keep the records (in the US, the local NSS grotto) and keep the locations among them.

## Install

Python 3.9 or newer, then:

    pip install numpy scipy scikit-image rasterio pyshp matplotlib

Put `sink_scout.py`, `county_scan.py`, and `field_kit.py` in one folder and run them from there.

## Self-test

    python sink_scout.py --selftest --out selftest

This builds a synthetic landscape with five planted valley sinks (two of them inside a road-dammed basin), a lone sink, an off-channel hillside pit, and a stock pond. It runs the detector and checks each planted feature. Expected result:
- the five valley sinks are on top as one chain of 5;
- the lone sink is next;
- the hillside pit is flagged off-channel;
- the pond is flagged shallow/wide;
- the dammed basin itself is gone.

## Input data

Any bare-earth DEM that GDAL reads (GeoTIFF, ERDAS .img), in a projected CRS with linear units.

**USGS 3DEP 1-meter tiles** (The National Map) are GeoTIFFs in 10-km squares. Each is named by the northwest corner of the square in units of 10 km: x71y433 covers eastings 710,000–720,000 and northings 4,320,000–4,330,000, UTM metres. The CRS is correct, and several tiles from one project merge with `--mosaic`. A county that crosses a UTM zone line comes as two projects with separate numbering (Mesa County's east end is in zone 13). The two projects overlap near the line, and county_scan removes the duplicates.

**Mesa County's delivered 2016 grids** (4-ft cells) are in UTM zone 12, NAD83, US survey feet, with elevations in feet, but their embedded metadata says otherwise. Use `--crs mesa`.

Clip with `--bbox LAT1,LON1,LAT2,LON2` (two opposite corners, decimal degrees, snapped to the source grid; works with `--mosaic`). Depressions touching the box edge are dropped, so draw it wide.

## Memory and cell size

About 180 bytes per processed cell. 8 GB free handles about 40 million cells: a whole 1-m tile at `--agg 2` (25 M), or a box about 6 km on a side at 1 m. `--agg N` averages N×N cells.

Results by cell size in the book's tests:
- **1 m:** the known cave at rank 1, pits and ponds cleanly separated.
- **2 m:** nearly as good; 38 of the top 40 on one hillside were pits in imagery.
- **4 m:** about half right, because small pits and ponds merge.

**Coarse then fine.** Sweep a block at `--agg 3` or `--agg 4` to find fields and long chains. The sweep writes `refine_commands.txt`, one full-resolution `--bbox` run per chain and per top candidate (`--refine-top`, default 15). Refine boxes extend `--refine-margin` (default 300 m) beyond the pits.

Catchment, the on-channel test, the creek flag, and chains are all measured inside the box. A pit on a wash that drains a square mile reports a few acres from a 600-m box (one known cave: 258 ac on the full mosaic, 6 ac in its box). When the verdict depends on the wash above, use `--refine-margin 1000`.

## Pipeline (sink_scout)

1. **Depression fill.** Flats are resolved, D8 flow directions computed, and flow accumulated.
2. **Nested split.** Pits with their own depth inside a larger basin (for example behind a road fill) are split out and measured from their own spill points.
3. **Rim trim.** Each pit is flooded upward until the wall slope falls below `--rim-slope`. That level is the rim.
4. **Measurement:**
   - depth, diameter, area, volume
   - depth/width ratio
   - elongation
   - catchment
   - on-channel test
   - road distance
   - geology distance
5. **Thin-dam test.** A grey opening with a k×k window (k ≈ `--dam-width`/cell + 1) removes raised strips thinner than `--dam-width` around each pit. The pit is then refilled in a local window with free-draining edges.
   - `dam_frac` is the fraction of depth that survives.
   - Below 0.5 the pit is flagged `thin-dam`, its score multiplied by 0.4, and it is excluded from chains.
   - The test is skipped on cells coarser than 2.7 m.
   - It catches fallen trees across ravines, fence lines and thin berms.
6. **Chains.** From each pit's spill point the flow path is traced downstream. If it reaches another pit within `--chain-dist` along the path, the two are linked. Each chain is drawn from the top pit's floor through every member's floor to the root.
7. **Score.**
   - Shape, log catchment and chain bonus add.
   - These flags multiply the score down: `off-channel`, `shallow/wide`, `creek`, `road`, `thin-dam`, `off-geology`, `nodata-edge`.
   - `nested(n)` marks a pit split out of a basin holding n of them; n ≥ 5 in one hollow is the talus-pit signature.

## Options that matter

| Option | Default | Meaning |
|---|---|---|
| `--min-depth` | 1.5 m | Smallest pit that counts |
| `--min-area` | 15 m² | Drops single-cell noise |
| `--min-ratio` | 0.05 | Depth/width below this is flagged `shallow/wide` |
| `--channel-ha` | 2 ha | Drainage needed for a pit to count as on a channel |
| `--chain-dist` | 200 m | Longest along-wash gap between linked pits |
| `--chain-ratio` | 0.1 | Depth/width a pit needs to join a chain |
| `--max-catch-ha` | 300 ha | Catchment above which a depression is a `creek` and kept out of chains |
| `--rim-slope` | 0.08 | Wall slope that ends the rim trim; 0 turns rim trimming off |
| `--dam-width` | 4 m | Thickest ridge treated as a possible false dam; 0 turns the test off |
| `--osm-roads` | off | Downloads OpenStreetMap roads and tracks for the footprint and flags culvert dams |
| `--roads file.shp` | none | Uses a local road layer instead of OpenStreetMap |
| `--geology geol.shp --geo-field UNIT --geo-units TRc,Jm` | none | Favors candidates near outcrops of the named units |
| `--check LAT,LON` | none | Reports the rank given to the pit at those coordinates; if it was dropped, reports its fill depth and area and which threshold to loosen. Plain decimals, no spaces: `39.08,-108.54` |
| `--save-rasters` | off | Writes dem.tif (the ground as processed), fill depth, flow accumulation and labels |
| `--no-coords` | off | Scale bar instead of coordinate axes; no lat/lon on close-ups |

Notes on the options:
- `--chain-ratio` is set below `--min-ratio` so that a partly collapsed slot still links to the pit beside it.
- `--chain-dist`: a test at 300 m on a known pit field linked nothing new, so keep 200.

## Outputs (sink_scout)

- `candidates.csv`: one row per pit, with every measurement in metres and feet, plus `flags`, `chain_id`, `depth_open_ft` and `dam_frac`.
- `candidates.kml` and `candidates.geojson`: the same data for Google Earth and for GIS.
- `candidates.shp` and `depressions.shp`: points and footprints, for GIS.
- `chains.csv`, `chains.kml`, `chains.geojson`: one line per chain, drawn down the wash. Each records:
  - the along-wash length;
  - the floor-to-floor and rim-to-rim drop;
  - the grade.
  
  These are a preliminary figure for the pipe underneath, never a surveyed one.
- `offchannel_pits.csv`: pit-shaped holes with no drainage into them, such as cellar holes, cisterns, prospects, mine shafts and lidar artifacts.
- `overview.png`: hillshade, channels in blue, ranked candidates numbered.
- `thumbs/`: a close-up of each top candidate.
- `refine_commands.txt`: written after any `--agg` run.

## county_scan.py

    python county_scan.py tiles --out county --dry-run
    python county_scan.py tiles --out county --agg 1 --split 3 --osm-roads

`--dry-run` prints the plan, the memory per run and the expected time, and runs nothing.

**How tiles are run:**
- Each tile (or each of its `--split` N×N pieces) is run with a 500-m buffer of same-project ground.
- A candidate is kept only if its centre falls in the tile's core. Cores follow the nominal 10-km grid, because USGS tiles overlap their neighbours by about 6 m.
- At 1 m, `--split 3` uses about 4 GB per run. A 100-tile county takes most of a day.
- Runs are resumable. A finished piece is rerun automatically if its source tiles change (a neighbour added) or if `--agg` changes.
- Across runs from different pieces, candidates less than 6 m apart are merged into one. This handles UTM-zone project overlap.

**What is kept per tile:** the top 50 pits on washes, the top 10 off-channel holes, and the top 10 on creeks. Ponds and culvert dams are dropped. Road data is fetched once per tile.

**Outputs:**
- `county_candidates.kml` (three groups plus chains), with the same data as CSV and GeoJSON
- `county_chains.geojson`
- an overview map
- `county_density_nocoords.png`

Keep `county/runs/`: field_kit reads chain members from it.

## field_kit.py

    python field_kit.py county/county_candidates.geojson --tiles tiles --chains county/county_chains.geojson --sheets top25 --out field

It writes:
- **A KMZ with a two-panel lidar close-up in every balloon.** Both panels are shaded relief with 0.5-m contours:
  - 300 m around the pit;
  - a 60-m zoom at full resolution, showing the hollow and its lip: whether it closes against a broad bank or a thin ridge, and how clean the rim is.
  
  Bare-earth models smooth boulders away; a talus pit shows as one hollow with many low spots.
- **A GPX file** of waypoints and chain tracks.
- **Chain profiles:** the ground down the wash, the pit floors, and the ground above the line joining them. These show relative heights only, so they are publishable. Members are read from `county/runs/<piece>/candidates.geojson`; pits the county list dropped are labelled `r12` and drawn in grey.
- **A field-sheets PDF:** one page per pit with the close-up, the numbers, the profile and a field checklist. It carries a PRIVATE banner.

Options:
- `--tag land=file.shp:FIELD` (repeatable) labels candidates by polygon, for land manager or geology. Even-odd point-in-polygon.
- `--truth verdicts.csv` (lat, lon, verdict, name, notes, date) writes:
  - `truth_report.txt`: the pipe rate by layer and rank band, and recall of independently known caves;
  - `ground_truth_worklist.csv`.
  
  Verdicts follow the pit, not its rank, across reruns.
- The z unit is auto-detected.

## Troubleshooting

- **The log reports elevations of 0.0 to 0.0.** The files were not read. Some tiles declare their nodata value as the float32 minimum. Current versions merge into NaN to handle this; update and rerun.
- **Out of memory.** Use a larger `--agg`, a smaller `--bbox`, or `--split` with county_scan.
- **Nothing found in a box.** Depressions touching the edge are dropped. Draw the box wider.
- **OpenStreetMap fetch fails.** The program retries a mirror. Otherwise pass a local `--roads` file.

## Repository rule

Never commit run outputs: no `county/`, `field/` or run folders, and no `ground_truth*.csv`. They carry locations.
