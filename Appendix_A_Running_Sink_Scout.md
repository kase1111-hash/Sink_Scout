# Appendix A

## Running sink_scout

**If you do not write software, read this page.**

sink_scout is the program from the chapter on finding caves from the sky. You give it a lidar elevation file, the free kind the government serves for most of the country, and it does the flood test on every square yard of it. It lists the closed pits that sit on drainages, links the ones strung down the same wash into chains, and ranks them by how much they look like the ground over a pipe. It writes the list as a spreadsheet, a picture of the ground with the pits numbered, and a file that opens in Google Earth with a pin on every pit.

It does not find caves. It finds places where the ground has been hollow. Most of the top of a good list will be real pits; plenty below that will be ponds, culverts, cellar holes, and things nobody has a name for, and the only way to know which is which is to look at each one in imagery and then walk the ones that survive.

Everything it writes is a map of places where the ground is hollow, and that means everything it writes is private. Look at it, walk it, report what you find to the grotto, and do not post it. The reasons are in the chapter on the people who map mud.

The rest of this appendix is one path from nothing to a first result. The program's full reference, every setting and every file it writes, is in the README at the repository address on the copyright page, where it can be kept current as the program changes.

---

**What you need**

A computer running Windows, Mac, or Linux, with 8 GB of memory or more. The runs in this book were made on a desktop with 8 GB free.

Python, which is free. Download it from python.org and install it. On Windows, tick the box that says Add Python to PATH on the first screen of the installer.

The three program files, sink_scout.py, county_scan.py, and field_kit.py, from the repository. Put them in one folder, say a folder called sinks in your Documents.

**Opening a terminal**

The program has no buttons. It runs from a terminal, a window where you type a line and press Enter. On Windows it is called Command Prompt; on a Mac it is Terminal, in Applications, Utilities. Open it and move into your folder by typing

    cd Documents/sinks

and pressing Enter. Then install the six libraries the program uses, with one line:

    pip install numpy scipy scikit-image rasterio pyshp matplotlib

It takes a few minutes and prints a great deal. If it says pip is not recognized, Python was installed without the PATH box ticked; run the installer again and tick it.

**The self-test**

    python sink_scout.py --selftest --out selftest

This builds a made-up hillside with sinkholes planted in it, a road dam, a stock pond, and a pit on a slope with no water running to it, runs the program on it, and checks each planted feature. It takes a few seconds. The lines to look for read like this:

    HIT: candidate rank #6 (score 3.127) — depth 15 ft, diam 30 ft, d/w 0.49, catch 10 ac, chain of 1
    HIT: candidate rank #8 (score 0.418) — depth 8 ft, diam 164 ft, d/w 0.05, catch 5 ac, chain of 1, flags: shallow/wide

The first is a lone sinkhole in a valley, found and ranked. The second is the stock pond, found and flagged shallow and wide, which pushes it to the bottom. At the end the program prints what it expected, and if what it found matches, the installation works. Nothing else in this appendix will work until this does.

**One real run**

Get a lidar file. In the United States the place is the USGS National Map Downloader, on the web at apps.nationalmap.gov/downloader. Find your area on the map, choose Elevation Products (3DEP) and then the 1 meter DEM, search, and download one of the tiles it offers. Each is a GeoTIFF about half a gigabyte, covering a square about six miles on a side. Put it in your sinks folder.

A whole tile is too much for most computers at full detail, so clip a box out of it. The box is two opposite corners in latitude and longitude, typed as decimals with a comma between and no spaces. This example is a city park in Grand Junction, picked because there is no cave under it:

    python sink_scout.py USGS_one_meter_x71y433_CO_MesaCo_QL2_UTM12_2016.tif --bbox 39.075,-108.545,39.085,-108.530 --out park_1m

Use your own tile's name and your own box. Anything touching the edge of the box is dropped, so draw it wider than the ground you care about; a box half a mile to a mile on a side runs in a minute or two. When it finishes, the program says where it put the results.

**What it wrote**

Open the folder it names. The files that matter at first:

overview.png is a picture of the ground, lit as if by a low sun, with the washes in blue and the candidates numbered by rank.

candidates.kml opens in Google Earth, with a pin on each candidate. Click a pin to see its measurements.

candidates.csv opens in any spreadsheet. One row per candidate, best first.

chains.kml draws each chain as a line down its wash.

offchannel_pits.csv is the side list of holes with no water running into them, where cellar holes, prospects, and mine shafts land.

The thumbs folder has a close-up of each of the top candidates.

**Reading a candidate**

Each row gives the rank and the score; the depth and the width, in feet and meters; how many acres drain into it; whether it sits on a channel; the chain it belongs to and how many pits are in that chain; and how far it is from a road. The flags column says what the program doubts. Off-channel means no water reaches it, so it is probably not a swallet. Shallow/wide means saucer rather than well, which is usually a pond. Road means it sits against a road fill, which is usually a culvert. Creek means too much water drains to it for a pipe, which is usually a dam. Thin-dam means most of its depth is held up by something narrow, a log or a berm. Nodata-edge means it touches a gap in the laser data. A candidate with no flags, a depth near half its width, and a few to a few dozen acres draining into it, in a chain down a wash, is the kind worth a walk.

Before walking to anything, find it in satellite imagery. Winter imagery with a low sun shows pits best. A dark hole with a lit crescent on one wall is what a pit looks like from above. If the imagery shows a pond, a tank, or a building, cross it off.

**When it goes wrong**

If the program reports elevations from 0.0 to 0.0, it did not read the file, whatever else it says; get the current version from the repository and run it again. If it runs out of memory, make the box smaller, or add --agg 2 to the line, which averages the ground into cells twice as wide and uses a quarter of the memory. If it finds nothing, the box is probably too tight around the ground you care about. Mesa County's own lidar files, the four-foot grids the county delivers, need --crs mesa added to the line, because their labels do not match their contents; the USGS one-meter tiles do not.

**Going further**

Three more things the README covers in full. county_scan.py runs a whole folder of tiles, a county if you have the disk and the patience, and merges the results into one list; the run in the chapter on finding caves from the sky took most of a day. field_kit.py turns a finished list into what goes in the truck: a Google Earth file with a lidar close-up in every pin, waypoints for a GPS or a phone, a profile of each chain, and printed field sheets with the checks from the chapter on walking them. And sink_scout.py --check followed by a latitude and longitude tells you what the program made of a pit you already know, which is the best way to learn what its numbers mean.

**What it mistakes for a pit**

Anything that holds water or shows a floor below the ground around it. In clay country that is stock ponds and culvert ponds, which the flags catch, and the occasional cellar hole or mine shaft, which lands on the side list. In the high timber, a fallen tree across a ravine dams it in the model, and a row of them reads as a chain until the thin-dam test drains them. In sandstone there is more: potholes, plunge pools under pour-offs, narrow slots, and arches, where the laser sees the floor through the opening and calls it a hole. The program does better on clay than on rock.

**The rule**

Walk the candidates, report what you find to the people who keep the records, and do not publish locations.
