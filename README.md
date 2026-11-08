# Sink Scout

`sink_scout`, `county_scan` and `field_kit` find collapse-sinkhole candidates in bare-earth lidar: pits on drainages, linked into chains down their washes, ranked by how much they look like the ground over a piping cave.

Written for the book *Hollowed Ground: The Piping Caves of Colorado's Western Slope* (Kase Branham, 2026). The book's Appendix A, [included here](Appendix_A_Running_Sink_Scout.md), is the beginner's path from nothing to a first result. This file is the reference: every setting and every file the programs write. When the programs change, it changes with them; [CHANGELOG.md](CHANGELOG.md) says what changed.

**The outputs of these programs are maps of places where the ground is hollow. Do not publish them.** The one exception is `county_density_nocoords.png`, which shows only counts per two-kilometer square. Figures for publication are made with `--no-coords`. Report what you find to the people who keep the records (in the US, the local NSS grotto) and keep the locations among them. More under [The rule](#the-rule).

## The three programs

| Program | Takes | Does |
|---|---|---|
| `sink_scout.py` | one elevation file, several merged, or a box clipped from them | Finds the pits, links them into chains, scores and ranks them, and writes the list, the maps and the Google Earth files. |
| `county_scan.py` | a folder of USGS 1-meter tiles | Runs `sink_scout.py` over every tile, a piece at a time, and merges the results into one county list and two county maps. |
| `field_kit.py` | a finished list from either, and the same tiles | Makes what goes in the truck: a Google Earth file with a lidar close-up of every pit, GPS waypoints, chain profiles, printed field sheets, and a running tally of what the pits turned out to be. |

The usual order:

```
python sink_scout.py --selftest --out selftest                       # once, to prove the install
python sink_scout.py tile.tif --bbox LAT1,LON1,LAT2,LON2 --out box1  # one box, to learn the numbers
python county_scan.py tiles --out county --dry-run                   # a county: the plan
python county_scan.py tiles --out county --osm-roads                 # a county: the run (hours)
python field_kit.py county/county_candidates.geojson --tiles tiles --chains county/county_chains.geojson --sheets top25 --out field
python field_kit.py ... --truth verdicts.csv                         # again as verdicts come in
```

## Contents

- [The rule](#the-rule)
- [What it finds](#what-it-finds)
- [Install](#install)
- [Self-test](#self-test)
- [A first run](#a-first-run)
- [Input data](#input-data)
- [Memory, time and cell size](#memory-time-and-cell-size)
- [How sink_scout works](#how-sink_scout-works)
- [sink_scout.py](#sink_scoutpy): options, outputs
- [county_scan.py](#county_scanpy): how tiles are run, layers, options, outputs
- [field_kit.py](#field_kitpy): options, outputs, ground truth, tags
- [Output columns](#output-columns)
- [Troubleshooting](#troubleshooting)
- [Known issues](#known-issues)
- [Files in this repository](#files-in-this-repository)
- [License and citation](#license-and-citation)

## The rule

Everything these programs write is a map of places where the ground is hollow, and that makes everything they write private. The reasons are in the book's chapter on the people who map mud.

**Private:** every file the programs write except the three kinds below, and their console output too: the ranked table prints latitude and longitude, and the log prints the `--bbox`.

**Publishable:**
- `county_density_nocoords.png`: counts of pits per two-kilometer square, with no points, labels or coordinates.
- Chain profiles from `field_kit.py` (`profiles/*.png`): heights along the wash relative to the top pit, with no coordinates.
- `overview.png` and close-ups made with `--no-coords`. That option takes the coordinates off a figure, not the ground. A wide hillshade of distinctive terrain can be matched to the public lidar, so crop to what the figure needs.

**Less obvious ways a location gets out:**
- A tile name with a rank (`x71y433, rank 3`): the lidar is public and the program gives the same list for the same tile, so anyone can run it and read off the pit.
- A `--check` or `--bbox` value in a question or a bug report. [CONTRIBUTING.md](CONTRIBUTING.md) says how to report a problem without one.
- A file or folder named after a place.

**Never commit run outputs.** The [.gitignore](.gitignore) keeps out the usual output folders (`county/`, `field/`, `field_*/`, `check_*/`, `runs/`, `selftest/`, `*_out/`, `*_1m/`), anything with `PRIVATE` in its name, elevation data, `ground_truth*.csv`, and every `.kml`, `.kmz` and `.gpx`. Run folders named any other way are not covered: keep them outside the repository.

Walk the candidates, report what you find to the people who keep the records, and do not publish locations.

## What it finds

A piping cave, a pipe in colluvium or bentonitic fill, shows at the surface as a chain of steep-sided collapse pits on the axis of a wash, with real catchment draining into them. A normal wash has no closed depressions at all, so anything that survives the filters is worth a look.

It does not find caves. It finds places where the ground has been hollow. Most of the top of a good list will be real pits; plenty below that will be ponds, culverts, cellar holes, and things nobody has a name for. The only way to know which is which is to look at each one in imagery and then walk the ones that survive.

It mistakes for a pit anything that holds water or shows a floor below the ground around it. In clay country that means stock ponds and culvert ponds, which the flags catch, and the occasional cellar hole or mine shaft, which lands on the off-channel list. In the high timber, a fallen tree across a ravine dams it in the model, and a row of them reads as a chain until the thin-dam test drains them. In sandstone there is more: potholes, plunge pools under pour-offs, narrow slots, and arches, where the laser sees the floor through the opening and calls it a hole. The program does better on clay than on rock.

Checked in October 2026 against one surveyor's list of 229 Mesa County caves, a county scan put a candidate within 100 m of the entrance of 124 of them. Many of the small caves in Wasatch claystone enter through pits on open slopes with no wash above them, so they land in the off-channel layer, not the main list.

## Install

Python 3.9 or newer, then, in the folder with the files:

    pip install -r requirements.txt

or, without the file:

    pip install numpy scipy scikit-image rasterio pyshp matplotlib

The rasterio wheels from pip include GDAL, so nothing else needs installing. Put `sink_scout.py`, `county_scan.py`, and `field_kit.py` in one folder and run them from there: `county_scan.py` runs the `sink_scout.py` beside it, and `field_kit.py` imports it. If `pip` is not recognized on Windows, Python was installed without "Add Python to PATH"; Appendix A covers the fix.

## Self-test

    python sink_scout.py --selftest --out selftest

This builds a synthetic landscape (4-ft cells in the Mesa County CRS, elevations in feet) with five planted valley sinks (two of them inside a road-dammed basin), a lone sink, an off-channel hillside pit, and a stock pond. It runs the detector, checks each planted feature with `--check`, and takes a few seconds. Expected result:
- the five valley sinks are on top as one chain of 5;
- the lone sink is next;
- the hillside pit is flagged off-channel;
- the pond is flagged shallow/wide;
- the dammed basin itself is gone.

The two lines Appendix A quotes are printed exactly so:

    HIT: candidate rank #6 (score 3.127) — depth 15 ft, diam 30 ft, d/w 0.49, catch 10 ac, chain of 1
    HIT: candidate rank #8 (score 0.418) — depth 8 ft, diam 164 ft, d/w 0.05, catch 5 ac, chain of 1, flags: shallow/wide

Nothing else will work until this does. The current version gives exactly these lines with Python 3.13, numpy 2.5, scipy 1.18, scikit-image 0.26, rasterio 1.5 and matplotlib 3.11. Besides the usual outputs, the self-test leaves its synthetic ground and road in the output folder as `selftest_dem.tif` and `selftest_roads.geojson`.

## A first run

The example from Appendix A is a city park in Grand Junction, picked because there is no cave under it:

    python sink_scout.py USGS_one_meter_x71y433_CO_MesaCo_QL2_UTM12_2016.tif --bbox 39.075,-108.545,39.085,-108.530 --out park_1m

Use your own tile and your own box. A box half a mile to a mile on a side runs in a minute or two at 1 m. When it finishes, the program says where it put the results. Open `overview.png` first, then `candidates.kml` in Google Earth, then read the flags in `candidates.csv` ([below](#flags)).

Before walking to anything, find it in satellite imagery. Winter imagery with a low sun shows pits best. A dark hole with a lit crescent on one wall is what a pit looks like from above. If the imagery shows a pond, a tank, or a building, cross it off.

The best way to learn what the numbers mean is to ask about a pit you already know:

    python sink_scout.py tile.tif --bbox LAT1,LON1,LAT2,LON2 --check LAT,LON --out known_pit

`--check` prints the rank the pit got. If it was dropped, it says so, gives the depth and area of the closed basin the fill found there, and names the limits it may have failed.

## Input data

Any bare-earth DEM that GDAL reads (GeoTIFF, ERDAS .img), in a projected CRS with linear units.

**USGS 3DEP 1-meter tiles** (The National Map) are GeoTIFFs in 10-km squares. Each is named by the northwest corner of the square in units of 10 km: x71y433 covers eastings 710,000–720,000 and northings 4,320,000–4,330,000, UTM meters. The CRS is correct, and several tiles from one project merge with `--mosaic`. A county that crosses a UTM zone line comes as two projects with separate numbering (Mesa County's east end is in zone 13). The two projects overlap near the line, and county_scan removes the duplicates.

**Mesa County's delivered 2016 grids** (4-ft cells) are in UTM zone 12, NAD83, US survey feet, with elevations in feet, but their embedded metadata says otherwise. Use `--crs mesa`, which sets that CRS and assumes feet for elevations (override with `--z-units m`).

**Other DEMs.** A DEM in latitude and longitude must be reprojected first. rasterio installs a command that does it:

    rio warp dem_lonlat.tif dem_utm.tif --dst-crs EPSG:26912 --resampling bilinear

If a file's CRS tag is wrong, override it with `--crs` (an EPSG code or a PROJ string). If its elevations are in feet, add `--z-units ft`. Nodata comes from the file; `--nodata` overrides it, and anything below `--zmin` (default -1000) is treated as nodata as well.

Clip with `--bbox LAT1,LON1,LAT2,LON2` (two opposite corners, decimal degrees, snapped to the source grid; works with `--mosaic`). Depressions touching the box edge are dropped, so draw it wide.

## Memory, time and cell size

About 180 bytes per processed cell, and about 4 seconds per million cells on the author's desktop (25 million cells in about 95 seconds). 8 GB free handles about 40 million cells: a whole 1-m tile at `--agg 2` (25 M), or a box about 6 km on a side at 1 m. `--agg N` averages N×N cells; a block needs more than half its cells to hold data, or it becomes nodata.

Results by cell size in the book's tests:
- **1 m:** the known cave at rank 1, pits and ponds cleanly separated.
- **2 m:** nearly as good; 38 of the top 40 on one hillside were pits in imagery.
- **4 m:** about half right, because small pits and ponds merge.

**Coarse then fine.** Sweep a block at `--agg 3` or `--agg 4` to find fields and long chains. The sweep writes `refine_commands.txt`, one full-resolution `--bbox` run per chain and per top candidate (`--refine-top`, default 15). Refine boxes extend `--refine-margin` (default 300 m) beyond the pits. Each line carries the sweep's DEM files, `--crs`, `--min-ratio` and road options; add any other options by hand.

Catchment, the on-channel test, the creek flag, and chains are all measured inside the box. A pit on a wash that drains a square mile reports a few acres from a 600-m box (one known cave: 258 ac on the full mosaic, 6 ac in its box). When the verdict depends on the wash above, use `--refine-margin 1000`.

## How sink_scout works

1. **Depression fill.** The ground is filled until water would run off (grayscale reconstruction); filled surface minus ground is the depth of water each cell would hold. Flats left by the fill are given a slight slope toward their outlets, D8 flow directions are computed, and flow is accumulated, so drainage runs through each basin and the flow into a basin is its catchment.
2. **Nested split.** Pits with their own depth inside a larger basin (for example behind a road fill) are split out and measured from their own spill points. A pit needs `--min-depth` of its own to be split out.
3. **Rim trim.** Each pit is flooded upward until the wall slope falls below `--rim-slope`. That level is the rim, and depth is measured from it.
4. **Measurement:**
   - depth, diameter, area, volume
   - depth/width ratio
   - elongation
   - catchment
   - on-channel test
   - road distance
   - geology distance

   Pits outside the depth and area limits are dropped here, and so are depressions touching the edge of the area.
5. **Thin-dam test.** A grey opening with a k×k window (k ≈ `--dam-width`/cell + 1) removes raised strips thinner than `--dam-width` around each pit. The pit is then refilled in a local window with free-draining edges.
   - `dam_frac` is the fraction of depth that survives.
   - Below 0.5 the pit is flagged `thin-dam`, its score multiplied by 0.4, and it is excluded from chains.
   - The test is skipped on cells coarser than 2.7 m.
   - It catches fallen trees across ravines, fence lines and thin berms.
6. **Chains.** From each pit's spill point the flow path is traced downstream. If it reaches another pit within `--chain-dist` along the path, the two are linked. Each pit links to at most one pit below it, so the links form trees; the longest path from a top pit to the bottom pit is the chain's stem. Each chain is drawn from the top pit's floor through every member's floor to the root. A pit joins a chain only if its depth/width is at least `--chain-ratio`, its catchment is at most `--max-catch-ha`, and it is not flagged `road` or `thin-dam`.
7. **Score.** Shape, depth, log catchment and the chain bonus add; elongation subtracts; each flag multiplies the score down. Candidates are ranked by score.

### The score

    points = depth + shape + catchment + chain − elongation      (never below 0.05)
    score  = points × the factor of every flag the pit carries

| Term | Adds | At most |
|---|---|---|
| depth | depth in meters ÷ 4 (a 4-m pit adds 1) | 1.5 |
| shape | depth/width ÷ 0.3 (a pit 0.3 as deep as it is wide adds 1) | 1.5 |
| catchment | log₁₀(catchment ÷ `--channel-ha` + 1) (18 ha adds 1 at the default) | 1.5 |
| chain | 0.75 for each other pit in the chain | 1.5 |
| elongation | subtracts 0.25 for each unit of elongation above 2 (length to width) | 0.5 |

The highest possible score is 6.

### Flags

| Flag | Set when | Score × | Usually means |
|---|---|---|---|
| `off-channel` | no channel of `--channel-ha` or more touches the footprint | 0.5 | No water reaches it, so it is probably not a swallet. |
| `shallow/wide` | depth/width below `--min-ratio` | 0.4 | Saucer rather than well: usually a pond. |
| `creek` | catchment above `--max-catch-ha` | 0.5 | Too much water drains to it for a pipe: usually a dam. |
| `road` | within `--road-buffer` of a road (only with `--roads` or `--osm-roads`) | 0.15 | It sits against a road fill: usually a culvert. |
| `thin-dam` | keeps less than half its depth when ridges thinner than `--dam-width` are removed | 0.4 | Most of its depth is held up by something narrow, a log or a berm. |
| `off-geology` | farther than `--geo-buffer` from the chosen units (only with `--geology`) | 0.6 | Not on the rock you asked for. |
| `nodata-edge` | touches a gap in the laser data | 0.7 | Part of it may be missing. |
| `nested(n)` | split out of a basin holding n pits | 1 | n ≥ 5 in one hollow is the talus-pit signature. |

A candidate with no flags, a depth near half its width, and a few to a few dozen acres draining into it, in a chain down a wash, is the kind worth a walk.

## sink_scout.py

    python sink_scout.py DEM [DEM ...] [options]
    python sink_scout.py --selftest --out selftest

### Options

**Input**

| Option | Default | Meaning |
|---|---|---|
| `DEM` | | One or more elevation files. More than one needs `--mosaic`. |
| `--out DIR` | `sink_scout_out` | Output folder; created if missing. |
| `--mosaic` | off | Merge several DEMs into one before processing. |
| `--bbox LAT1,LON1,LAT2,LON2` | whole file | Clip to a box: two opposite corners in decimal degrees, no spaces, snapped to the source grid. |
| `--agg N` | 1 | Average N×N cells into one before processing (`--agg 3` for 1-ft grids). |
| `--crs CRS` | from the file | Override the file's CRS: `EPSG:26912`, a PROJ string, or `mesa`. |
| `--z-units m` or `ft` | `m` (`ft` with `--crs mesa`) | Units of the elevations. |
| `--nodata VALUE` | from the file | Override the nodata value. |
| `--zmin Z` | -1000 | Elevations below this are nodata. |

**Detection**

| Option | Default | Meaning |
|---|---|---|
| `--min-depth` | 1.5 m | Smallest pit that counts, and the depth of its own a pit inside a larger basin needs to be split out |
| `--max-depth` | 60 m | Deepest pit kept |
| `--min-area` | 15 m² | Drops single-cell noise |
| `--max-area` | 50,000 m² | Largest footprint kept |
| `--floor` | 0.3 m | Fill depth that puts a cell inside a depression's footprint |
| `--min-ratio` | 0.05 | Depth/width below this is flagged `shallow/wide` (county_scan passes 0.2) |
| `--channel-ha` | 2 ha | Drainage needed for a pit to count as on a channel; also the unit of the catchment score |
| `--chain-dist` | 200 m | Longest along-wash gap between linked pits |
| `--chain-ratio` | 0.1 | Depth/width a pit needs to join a chain |
| `--max-catch-ha` | 300 ha | Catchment (about 740 acres) above which a depression is a `creek` and kept out of chains |
| `--rim-slope` | 0.08 | Wall slope (rise over run) that ends the rim trim; 0 turns rim trimming off |
| `--dam-width` | 4 m | Thickest ridge treated as a possible false dam; 0 turns the test off. Skipped on cells coarser than `--dam-width` ÷ 1.5 |
| `--keep-edge` | off | Keep depressions touching the edge of the area. Their true size is unknown, so they are dropped by default |
| `--eps` | 0.00001 m | Slope given to flats so water can cross them; capped below the smallest real step in the DEM |

**Context layers**

| Option | Default | Meaning |
|---|---|---|
| `--osm-roads` | off | Downloads OpenStreetMap roads and tracks for the footprint and flags culvert dams. Cached as `osm_roads.geojson` in the output folder and reused |
| `--roads file.shp` | none | Uses a local road layer (lines or polygons, shapefile or GeoJSON) instead of, or as well as, OpenStreetMap |
| `--road-buffer` | 30 m | Distance from a road that flags `road` and keeps a pit out of chains |
| `--geology geol.shp` | none | Geologic map polygons (shapefile or GeoJSON) |
| `--geo-field` | `UNIT` | Attribute holding the unit code |
| `--geo-units TRc,Jm` | every polygon | Units to favor, comma-separated; case is ignored |
| `--geo-buffer` | 300 m | Distance from those units that still counts as on them; farther is flagged `off-geology` |

**Outputs**

| Option | Default | Meaning |
|---|---|---|
| `--thumbs N` | 25 | Close-ups of the top N candidates; 0 for none (8 at most in the self-test) |
| `--thumb-size` | 300 m | Width of each close-up |
| `--no-overview` | off | Skip `overview.png` |
| `--no-coords` | off | Scale bar instead of coordinate axes; no lat/lon on close-ups |
| `--save-rasters` | off | Writes dem.tif (the ground as processed), fill depth, flow accumulation and labels |

**Follow-up**

| Option | Default | Meaning |
|---|---|---|
| `--check LAT,LON` | none | Reports the rank given to the pit at those coordinates; if it was dropped, reports the depth and area of the basin there and which limits to loosen. Plain decimals, no spaces: `39.08,-108.54` |
| `--refine-top` | 15 | How many unchained top candidates get a line in `refine_commands.txt` |
| `--refine-margin` | 300 m | Margin around each chain or pit in `refine_commands.txt` |
| `--selftest` | off | Builds the synthetic landscape described under [Self-test](#self-test) and runs on it |

Notes on the options:
- `--chain-ratio` is set below `--min-ratio` so that a partly collapsed slot still links to the pit beside it.
- `--chain-dist`: a test at 300 m on a known pit field linked nothing new, so keep 200.
- A coordinate that begins with a minus sign (a southern latitude) must be joined to its option with `=`: `--bbox=-33.91,151.20,-33.90,151.22`.

### Outputs

All in `--out`.

| File | Contents |
|---|---|
| `candidates.csv` | One row per pit, best first, with every measurement in meters and feet, plus `flags`, `chain_id`, `depth_open_ft` and `dam_frac`. [Columns](#candidatescsv). |
| `candidates.geojson` | The same data in longitude/latitude: a point at each pit's deepest cell and a polygon for its footprint. For GIS. |
| `candidates.kml` | The same for Google Earth: a pin per pit named by rank, with its measurements, and its outline. Red for ranks 1–10, orange to 30, yellow below. |
| `candidates.shp`, `depressions.shp` | Points at the deepest cell and footprint polygons, in the DEM's CRS, with `.shx`, `.dbf` and `.prj`. For GIS. Field names are cut to 10 characters. |
| `offchannel_pits.csv` | Pit-shaped holes (depth/width 0.25 or more) with no drainage into them, deepest first: cellar holes, cisterns, prospects, mine shafts, lidar artifacts. Same columns. |
| `chains.csv`, `chains.kml`, `chains.geojson` | One line per chain, drawn down the wash. Written when there is at least one chain. [Columns](#chainscsv). |
| `overview.png` | Hillshade, channels in blue, footprints in red, chains in orange, candidates colored by score; the top 40 numbered. |
| `thumbs/rank_001.png` … | A close-up of each top candidate: hillshade, channels, the pit's outline in red and its neighbors' in orange, its deepest cell marked, its measurements in the title. |
| `refine_commands.txt` | Written after any `--agg` run that found something: one full-resolution command per chain and per top candidate. |
| `osm_roads.geojson` | With `--osm-roads`: the roads downloaded for this area. |
| `dem.tif`, `fill_depth.tif`, `flowacc_ha.tif`, `labels.tif` | With `--save-rasters`: the ground as processed (meters, after `--agg`, `--bbox` and unit conversion), the fill depth, the drainage area in hectares, and each pit's footprint numbered by its `label`. `dem.tif` is the file to send, privately, when a single pit needs a second opinion. |
| `selftest_dem.tif`, `selftest_roads.geojson` | With `--selftest`: the synthetic ground and its road. |

The chain lines record the along-wash length, the floor-to-floor and rim-to-rim drop, and the grade. They measure wash between pits, not cave. Neither chain length nor sink size says much about how much cave is there; the walk decides.

### What it prints

Timestamped progress (the DEM's size, cell size, share of valid cells, elevation range and CRS, and the count left after each step), then a table of the top 15 candidates (rank, score, depth and diameter in feet, depth/width, catchment in acres, chain size, latitude, longitude, flags), the top 12 chains, the number of refine commands and off-channel pits written, and the `--check` report if one was asked for. Treat it as private: the table carries coordinates and the log carries the box.

## county_scan.py

    python county_scan.py tiles --out county --dry-run
    python county_scan.py tiles --out county --agg 1 --split 3 --osm-roads
    python county_scan.py tiles --out county --agg 1 --split 3 --merge-only        # rebuild the lists only

`--dry-run` prints the plan, the memory per run and the expected time, and runs nothing.

The folder holds USGS 1-m tiles under the names The National Map gives them (`USGS_one_meter_x71y433_<project>.tif`); other files are ignored. Mesa County's own grids are not named that way: run them with `sink_scout.py --mosaic --crs mesa`.

**How tiles are run:**
- Each tile (or each of its `--split` N×N pieces) is run with a 500-m buffer of same-project ground.
- A candidate is kept only if its center falls in the tile's core. Cores follow the nominal 10-km grid, because USGS tiles overlap their neighbors by about 6 m.
- At 1 m, `--split 3` uses about 4 GB per run. A 100-tile county takes most of a day.
- Runs are resumable. A finished piece is rerun automatically if its source tiles change (a neighbor added) or if `--agg` changes.
- Across runs from different pieces, candidates less than 6 m apart are merged into one. This handles UTM-zone project overlap.

Each run is one `sink_scout.py` command: the piece and its buffer as a `--bbox`, `--agg`, `--min-ratio`, `--max-catch-ha`, `--thumbs 0 --no-overview`, the tile's roads, and anything in `--extra`. Its command line and everything it prints go to `runs/<piece>/run.log`; a `DONE` file marks it finished. Only the tile set and `--agg` are compared when deciding whether to rerun. After changing `--extra`, `--min-ratio`, `--max-catch-ha` or the roads, delete `runs/` (or the `DONE` files) to redo the runs. Settings that only affect the merge (the `--top` limits, `--keep-shallow`) need only `--merge-only`.

Stop it at any time (Ctrl+C) and run the same command again; finished pieces are skipped. A run that fails is reported at the end; fix the cause and run the same command again.

**What is kept:** every pit on a wash, every off-channel hole, and every pit on a creek, with no per-tile limit (set one with `--top`, `--top-offchannel`, `--top-bigwash` if a county is too crowded). Road-flagged depressions are dropped. Shallow, wide hollows are dropped unless `--keep-shallow`, which writes them as their own layer. Road data is fetched once per tile. To rebuild the county files from finished runs with different settings, add `--merge-only`.

*Why no limits (Oct 2026):* checked against one surveyor's list of 229 Mesa County caves, the old limits (50 / 10 / 10 per tile) cut the caves with a candidate within 100 m of an entrance from 124 to 91, mostly small caves in Wasatch claystone that enter through pits on open slopes with no wash above them and so land in the off-channel layer. Keeping shallow hollows added 2,437 candidates for 8 more caves within 25 m, which is why it is optional. The score itself held up: a refit on the same caves did no better out of sample (AUC 0.63 against 0.72 for the current score).

county_scan passes `--min-ratio 0.2` to each run, not sink_scout's own 0.05, so a pit shallower than a fifth of its width is `shallow/wide` here.

### Layers

Each candidate goes to the first layer that fits, in this order:

| Layer | Ids | Holds | Limit per tile |
|---|---|---|---|
| (dropped) | | flagged `road` | |
| `shallow` | S1, S2 … | flagged `shallow/wide`, only with `--keep-shallow` (dropped otherwise): mostly ponds, a few caves | `--top-shallow` |
| `offchannel` | O1, O2 … | flagged `off-channel`: shafts, cellars, arches, and in claystone hillsides many real cave entrances | `--top-offchannel` |
| `bigwash` | B1, B2 … | flagged `creek`: on a channel draining more than `--max-catch-ha` | `--top-bigwash` |
| `pits` | P1, P2 … | everything else: pits on washes | `--top` |

County ranks are by score within each layer. `field_kit.py` uses the same ids.

### Options

| Option | Default | Meaning |
|---|---|---|
| `folder` | | Folder of USGS 1-m tiles |
| `--out DIR` | `county` | Output folder |
| `--agg N` | 2 | Cell size of the runs in meters (with 1-m tiles) |
| `--split N` | 1 | Split each tile into N×N pieces; use 3 with `--agg 1` so each run fits in memory |
| `--buffer` | 500 m | Ground borrowed from neighboring tiles around each piece |
| `--osm-roads` | off | Fetch OpenStreetMap roads once per tile and use them in every run |
| `--min-ratio` | 0.2 | Passed to each run |
| `--max-catch-ha` | 300 ha | Passed to each run |
| `--extra "..."` | none | More sink_scout options, quoted: `--extra "--min-depth 2"`. Split on spaces, so paths in it cannot contain spaces |
| `--top N` | 0 (all) | Pits on washes kept per tile |
| `--top-offchannel N` | 0 (all) | Off-channel pits kept per tile |
| `--top-bigwash N` | 0 (all) | Big-wash pits kept per tile |
| `--keep-shallow` | off | Keep shallow, wide hollows as their own layer instead of dropping them |
| `--top-shallow N` | 0 (all) | Shallow hollows kept per tile |
| `--only IDS` | all tiles | Run only these tiles, by id, comma-separated (`x70y433,x71y433`): for a test. A tile number used by two UTM-zone projects gets `_z12` or `_z13` added |
| `--dry-run` | off | Print the plan and estimates; run nothing |
| `--merge-only` | off | Skip the runs; rebuild the county files and maps from finished runs |
| `--no-maps` | off | Skip the two county maps |
| `--map-res` | 30 m | Hillshade cell for the county maps |
| `--map-dpi` | 200 | Resolution of the county maps |
| `--density-km` | 2 km | Square size of the publishable density map |
| `--sink-scout PATH` | beside this script | sink_scout.py to run |
| `--python PATH` | this Python | Python used for the runs |

The plan printed at the start of every run, dry or not, warns when the largest run is over 30 million cells (try a larger `--split` or `--agg`). It also notes when the tiles come from more than one project: projects (UTM zones) are not mosaicked together, so pits on the zone line are measured from one side only.

### Outputs

| File | Contents |
|---|---|
| `runs/<tile>/` or `runs/<tile>_r<i>c<j>/` | One sink_scout run per tile or piece, as listed under [sink_scout outputs](#outputs) but without the overview and close-ups, plus `run.log` and `DONE` |
| `roads/<tile>.geojson` | With `--osm-roads`: roads for the tile and its buffer |
| `county_candidates.geojson` | Kept candidates as points and footprints, in longitude/latitude, with every sink_scout field (except `label`, `parent`, `x`, `y`) plus `tile`, `piece`, `layer`, `run_rank` and `county_rank` |
| `county_candidates.csv` | The same table without shapes. [Columns](#county_candidatescsv) |
| `county_candidates.kml` | Google Earth: a folder per layer, pins named by id and colored orange to yellow by score, and a folder of chains |
| `county_chains.geojson` | Chains whose top pit lies in a run's core, with `tile` and `piece` |
| `county_overview.png` | Hillshade of the whole area with every kept pit, off-channel and big-wash candidate, and chain; the top 100 pits numbered. Private |
| `county_density_nocoords.png` | Pits on washes and on big washes per two-kilometer square over the hillshade, with a 10-km scale bar and no coordinates. Publishable |
| `county_summary.txt` | Per tile: runs finished, candidates found, how many fell inside the core, and how many were kept of those available in each layer |

Keep `county/runs/`: field_kit reads chain members from it.

## field_kit.py

    python field_kit.py county/county_candidates.geojson --tiles tiles --chains county/county_chains.geojson --sheets top25 --out field
    python field_kit.py county/county_candidates.geojson --tiles tiles --chains county/county_chains.geojson \
        --truth verdicts.csv --tag land=sma.shp:ADMIN_AGEN --tag geology=geology.shp:UNIT_NAME --sheets top25 --out field

It reads a candidate list and the elevation tiles the scan used. Nothing here reruns the scan, so it works on any finished run and can be run again whenever a verdict or a map layer is added.

The list can be a county list or one sink_scout run's `candidates.geojson` (with that run's `chains.geojson`). From a single run every candidate is included, ids are `#1`, `#2` …, and the layers are only `pits` (on a channel) and `offchannel`.

### Options

| Option | Default | Meaning |
|---|---|---|
| `candidates` | | `county_candidates.geojson`, or a run's `candidates.geojson` |
| `--tiles DIR` | required | Folder with the elevation tiles the scan used (`.tif` and `.img`) |
| `--chains FILE` | none | `county_chains.geojson` or a run's `chains.geojson`: profiles and tracks |
| `--out DIR` | `field` | Output folder |
| `--crs CRS` | from the tiles | CRS override for tiles with bad metadata (`mesa`) |
| `--layers` | `pits,offchannel,bigwash` | Layers to include. Add `shallow` for a county run with `--keep-shallow`. Verdicts are matched only against the layers included |
| `--sheets` | none | Field sheets: `top25` (the first 25 pits on washes), a chain (`C3`: its members in the list), or ids (`P31,P36,O4`), comma-separated. Quote single-run ids, since a shell reads `#` as a comment: `--sheets "#3,#7"` |
| `--truth FILE` | none | CSV of verdicts: lat, lon, verdict, and optionally name, notes, date |
| `--match-m` | 40 m | How close a verdict must be to a candidate |
| `--tag NAME=FILE:FIELD` | none | Label each candidate by the polygon it falls in; repeatable |
| `--closeup-m` | 300 m | Width of the close-up |
| `--zoom-m` | 60 m | Width of the zoom panel |
| `--no-closeups` | off | Skip the close-ups: fast, and a KMZ without pictures |

The z unit is auto-detected: if the tiles' elevations are about 3.28 times the list's, they are feet and are converted.

### Outputs

| File | Contents |
|---|---|
| `field_candidates.kmz` | For Google Earth Pro: a pin for every candidate, one folder per layer, with the close-up and the measurements in its balloon; the chains, each with its profile |
| `field_candidates_web.kml` | The same pins and lines without pictures |
| `field_candidates.gpx` | A waypoint for every candidate, named by its id, with the rim elevation and the depth, width, catchment, tags and verdict in its description; each chain as a track |
| `field_candidates.csv` | The list with close-up names, tags and verdicts. [Columns](#field_candidatescsv) |
| `closeups/<id>.jpg` | The two-panel close-up (PNG if JPEG cannot be written). A single run's `#12` is saved as `r12.jpg` |
| `profiles/<chain>.png` | The chain profile. Publishable |
| `profiles/chains.csv` | Per chain: the pits on the line, wash length, floor drop, grade, and the least and median ground above the line joining the floors |
| `sheets.pdf` | With `--sheets`: one page per pit for the clipboard |
| `truth_report.txt` | With `--truth`: how often each kind of candidate turned out to be a pipe, and which known caves the scan found or missed |
| `ground_truth_worklist.csv` | Every candidate, with the verdicts already given filled in. [Columns](#ground_truth_worklistcsv) |

**The close-up.** Both panels are shaded relief with contours:
- 300 m around the pit: contours every 1 m (2 m where the window spans more than 40 m of relief, 5 m beyond 100 m), every fifth one bold, with the pit's footprint in red, the other candidates in orange with their ids, and a yellow box marking the zoom.
- A 60-m zoom at full resolution with 0.5-m contours, showing the hollow and its lip: whether it closes against a broad bank or a thin ridge, and how clean the rim is.

Bare-earth models smooth boulders away; a talus pit shows as one hollow with many low spots. The title gives the id, depth, width and catchment, the flags, the depth without thin ridges when that is 2 ft or more shallower, and any verdict.

**The balloon** lists depth, width, depth/width, catchment, floor elevation, depth without thin ridges, chain, distance to the nearest road or track, the boulder-pit note when five or more low spots share one hollow, the tags, the flags and the verdict.

**Pin colors** in the KMZ and web KML:

| Color | Meaning |
|---|---|
| green | a pipe: verdict `pipe`, `known`, `walked`, `surveyed` or `cave` |
| blue-green | another real hole: `talus`, `mine`, `shaft`, `cellar`, `arch`, `alcove`, `slot`, `plunge`, `pothole` or `sinkhole` |
| red | not a hole: `pond`, `reservoir`, `culvert`, `road`, `logs`, `log`, `dam`, `artifact`, `building`, `none` or `nothing` |
| orange to yellow | no verdict yet; yellow is the highest score in its layer |

**Google Earth on the web** cannot show pictures packed inside a KMZ ("Only external images are currently supported"), and hosting the close-ups online would publish the locations. Use the .kmz in Google Earth Pro (desktop) for the pictures; `field_candidates_web.kml` gives the web version the pins without them.

**Chain profiles** show the ground down the wash, the pit floors, the dashed line joining them, and, shaded between pits, the ground above that line, labeled with its least and median height. Distance is feet along the wash from the top pit; height is feet below the top pit's rim. These show relative heights only, so they are publishable. The shading is not roof thickness, and the profile is not a survey. Members are read from `county/runs/<piece>/candidates.geojson`; pits the county list dropped are labeled `r12` and drawn in grey.

**Field sheets** carry a PRIVATE banner and, for each pit: id, layer and rank; the 300-m close-up and a 60-m zoom; latitude and longitude, and grid coordinates in the tile's CRS; floor and rim elevations, depth, width, depth/width, catchment, depth without thin ridges, chain, distance to the nearest road or track (straight line), tags, flags, the boulder-pit note, and any verdict; the chain profile when the pit is on one; and the field checks from the book's safety chapter with boxes to tick.

### Ground truth

Keep the verdicts in a CSV with a header row. `lat`, `lon` and `verdict` are required; `name`, `notes` and `date` are optional. Column names can be in any order and any case, other columns are ignored, and the UTF-8 file Excel saves is fine:

    lat,lon,verdict,name,notes,date
    LAT,LON,pipe,,walked 120 ft to a choke,2026-10-01
    LAT,LON,pond,,stock tank,2026-10-01
    LAT,LON,known,Lone Sink Cave,on the grotto list,2026-09-12

The first word of the verdict counts; case is ignored:

| Class | Words |
|---|---|
| pipe | `pipe`, `known`, `walked`, `surveyed`, `cave` |
| other real hole | `talus`, `mine`, `shaft`, `cellar`, `arch`, `alcove`, `slot`, `plunge`, `pothole`, `sinkhole` |
| not a hole | `pond`, `reservoir`, `culvert`, `road`, `logs`, `log`, `dam`, `artifact`, `building`, `none`, `nothing` |
| not counted | anything else, such as `unknown` |

Use `known` only for caves known before or apart from the program (on the list, surveyed by someone, found on foot): those are the ones that measure what it misses.

Each verdict goes to the nearest candidate within `--match-m`. Verdicts follow the pit, not its rank, across reruns. When two rows land on one candidate, the first row with a counted verdict wins, so change a verdict in place instead of adding a row below it.

`truth_report.txt`, also printed, gives:
- how many verdicts matched a candidate;
- for each layer, the candidates checked and how many were pipes, other real holes, and not holes (with percentages for pipes and not-holes), and a tally of the verdict words;
- for pits on washes, the pipe rate and the real-hole rate by rank band (1–25, 26–50, 51–100, 101–200, 201 on);
- for every cave marked `known`, whether the scan found it (with the candidate and its distance from the recorded point) or missed it (with the nearest candidate and the `--check` command that shows why).

`ground_truth_worklist.csv` lists every candidate with the verdicts given so far, followed by any verdicts that matched no candidate. Add verdicts as you check pits and pass it as the next `--truth` file.

### Tags

`--tag NAME=FILE:FIELD` labels each candidate with the `FIELD` value of the polygon it falls in, for land manager, geology, or anything else drawn as polygons. Repeat it for several layers:

    --tag land=sma.shp:ADMIN_AGEN --tag geology=geology.shp:UNIT_NAME

`FILE` is a shapefile (its `.prj` is read; without one, longitude/latitude is assumed) or GeoJSON in longitude/latitude. Point-in-polygon uses the even-odd rule, so holes are honored; a candidate in no polygon is tagged `(none)`. If `FIELD` is missing or misspelled, the program stops and lists the fields the file has. Tags appear in the balloons, the waypoints, both CSVs and the field sheets, and the log prints how many pits fall under each value.

## Output columns

### candidates.csv

The same fields are in `candidates.geojson`, `candidates.kml` (in part), `offchannel_pits.csv`, and the two shapefiles (names in brackets where the shapefile shortens them). In the CSV and shapefiles a value that was not measured is -1; in GeoJSON it is empty (null).

| Column | Meaning |
|---|---|
| `rank` | Position by score, 1 = best |
| `score` | See [The score](#the-score) |
| `depth_m`, `depth_ft` | Rim to the deepest cell |
| `diam_m`, `diam_ft` | Diameter of a circle with the footprint's area |
| `area_m2` | Footprint area |
| `volume_m3` [`vol_m3`] | Water the pit would hold to its rim |
| `ratio` | Depth/width |
| `elong` | Length to width of the footprint; 1 is round |
| `catch_ha`, `catch_ac` | Area draining into the pit, inside the processed area |
| `on_channel` [`on_chan`] | 1 if a channel touches the footprint |
| `road_m` | Distance to the nearest road; -1 without a road layer |
| `geo_m` | Distance to the chosen geology; -1 without `--geology` |
| `chain_id` | Chain number, 1 = longest stem; 0 = not in a chain |
| `chain_n` | Pits in the chain; 1 = alone |
| `chain_len_m` [`chain_m`] | Along-wash length of the chain's stem |
| `min_z_m` | Elevation of the floor |
| `spill_z_m` | Elevation of the rim the depth is measured from |
| `depth_open_ft` [`dep_opn_ft`] | Depth left when ridges thinner than `--dam-width` are removed; -1 if the test did not run |
| `dam_frac` | `depth_open` ÷ depth; below 0.5 is `thin-dam`; -1 if the test did not run |
| `lat`, `lon` | The deepest cell, WGS84 |
| `x`, `y` | The deepest cell in the DEM's CRS and units |
| `flags` | See [Flags](#flags) |
| `label` | The pit's number in `labels.tif` |
| `parent` | The closed depression the pit came from; pits sharing one were split out of one basin |
| `nsib` | How many pits were split out of that depression; 1 = the pit is the whole depression |

### chains.csv

| Column | Meaning |
|---|---|
| `chain_id` | 1 = longest stem |
| `n_pits` | Pits in the chain |
| `ranks` | Their ranks, space-separated |
| `top_rank`, `bottom_rank` | The pit at the top of the stem and the one at the bottom |
| `stem_ft`, `stem_m` | Along-wash length from the top pit to the bottom pit |
| `straight_ft` | Straight-line distance between them |
| `floor_drop_ft` | Top pit's floor to bottom pit's floor |
| `rim_drop_ft` | Top pit's rim to bottom pit's rim |
| `grade_pct` | Rim drop ÷ stem length |
| `top_lat`, `top_lon`, `bottom_lat`, `bottom_lon` | The two ends |

`chains.geojson` carries `chain_id`, `n_pits`, `stem_ft`, `floor_drop_ft`, `rim_drop_ft`, `grade_pct` and `ranks`; `county_chains.geojson` adds `tile` and `piece`.

### county_candidates.csv

`layer`, `county_rank`, `tile`, `piece` (the run it was kept from), `run_rank` (its rank in that run), then `score`, `depth_ft`, `diam_ft`, `ratio`, `elong`, `catch_ac`, `road_m`, `chain_n`, `chain_len_m`, `min_z_m`, `flags`, `lat`, `lon` as in `candidates.csv`. `road_m` is empty without roads.

### field_candidates.csv

| Column | Meaning |
|---|---|
| `id` | P12, O3, B1, S7, or `#12` from a single run |
| `layer` | `pits`, `offchannel`, `bigwash` or `shallow` |
| `rank` | County rank in its layer, or the run's rank |
| `score`, `depth_ft`, `diam_ft`, `ratio`, `catch_ac`, `flags` | As in `candidates.csv` |
| `chain`, `chain_n` | The chain's id in this kit (C1 = longest) and its size |
| `floor_ft` | Elevation of the floor |
| `road_ft` | Distance to the nearest road or track; empty without roads |
| `depth_open_ft` | Depth without thin ridges |
| `tag_<name>` | One column per `--tag` |
| `verdict` | From `--truth` |
| `lat`, `lon` | The deepest cell |
| `closeup` | The close-up's file name in `closeups/` |

### ground_truth_worklist.csv

`id`, `layer`, `rank`, `score`, `depth_ft`, `diam_ft`, `catch_ac`, `chain_n`, `flags`, one `tag_<name>` column per `--tag`, `lat`, `lon`, `verdict`, `name`, `notes`, `date`. Verdicts that matched no candidate follow, with the layer `(no candidate)`.

### profiles/chains.csv

`chain`, `n_pits`, `pits_on_line` (ids along the line from the top), `stem_ft`, `floor_drop_ft`, `grade_pct`, `least_cover_ft` and `median_cover_ft` (the least, and the median, ground above the line joining the pit floors between pits: not roof thickness).

## Troubleshooting

- **The log reports elevations of 0.0 to 0.0.** The files were not read. Some tiles declare their nodata value as the float32 minimum. Current versions merge into NaN to handle this; update and rerun.
- **Out of memory.** Use a larger `--agg`, a smaller `--bbox`, or `--split` with county_scan.
- **Nothing found in a box.** Depressions touching the edge are dropped. Draw the box wider. If a pit you know is missing, run with `--check` on it.
- **OpenStreetMap fetch fails.** The program tries three Overpass servers in turn and, if all fail, carries on without roads and says so. Otherwise pass a local `--roads` file.
- **"Several DEMs given: add --mosaic".** Merge them with `--mosaic`, or run them one at a time.
- **"DEM is in geographic (lat/lon) coordinates".** Reproject it (see [Input data](#input-data)), or pass `--crs` if the tag is simply wrong.
- **"DEM carries no CRS".** Pass `--crs` (`EPSG:26912`, a PROJ string, or `mesa`).
- **"--bbox does not overlap this DEM".** The corners are latitude first, then longitude, in decimal degrees.
- **"argument --bbox: expected one argument".** A value that begins with a minus sign needs `=`: `--bbox=-33.91,151.20,-33.90,151.22`.
- **county_scan finds no tiles.** Tile names must contain `x<col>y<row>_` and end in `.tif`, as USGS names them.
- **"field_kit.py needs sink_scout.py in the same folder".** Keep the three programs together.
- **A close-up says "no elevation data here".** The `--tiles` folder does not hold the ground under that candidate.
- **No pictures in Google Earth on the web.** Expected: open the .kmz in Google Earth Pro, or use `field_candidates_web.kml` without them.
- **A thin-dam test "skipped" message.** The cells are coarser than `--dam-width` ÷ 1.5 (2.7 m at the default); the test needs finer ground.

## Known issues

Found while checking this version; none changes which pits are found or how they score.

- **Chains that cross a tile or piece seam.** field_kit matches chain members to the county list by the run they came from. A member the county list kept from the neighboring run is labeled by its run rank (`r1`, in grey) as though it had been dropped, is not marked as part of the chain, and is left out of that chain's `--sheets C<n>` pages.
- **`county_summary.txt`** has no per-tile column for the shallow layer; its first lines give the total.
- **Field sheets** use a 60-m zoom whatever `--zoom-m` says.
- **`--help` text** is out of date in places: county_scan's lists a `county_chains.kml` that is not written (the chains are a folder in `county_candidates.kml`) and says the summary includes failures and run times; field_kit's describes single-panel PNG close-ups with 1-m contours. This file describes what the programs do now.

## Files in this repository

| File | What it is |
|---|---|
| `sink_scout.py` | The detector |
| `county_scan.py` | County-wide runs over USGS tiles |
| `field_kit.py` | Field files, sheets and the ground-truth tally |
| `requirements.txt` | The six libraries |
| `Appendix_A_Running_Sink_Scout.md` | The book's Appendix A: the beginner's path |
| `CHANGELOG.md` | What changed, and when |
| `CONTRIBUTING.md` | Reporting problems and proposing changes without giving away locations |
| `CITATION.cff` | How to cite the programs and the book |
| `LICENSE` | CC0 1.0 |
| `.gitignore` | Keeps run outputs, tiles and ground-truth files out of the repository |
| `.github/ISSUE_TEMPLATE/problem_report.md` | The issue form, which asks for no locations |

## License and citation

CC0 1.0 Universal: dedicated to the public domain. Use, change, and share it freely; see `LICENSE`. The rule about locations above is a request, not a license term.

If you use the programs in published work, cite the book: Kase Branham, *Hollowed Ground: The Piping Caves of Colorado's Western Slope* (2026). GitHub's "Cite this repository" button gives the same from [CITATION.cff](CITATION.cff).
