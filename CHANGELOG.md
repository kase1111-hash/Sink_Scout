# Changelog

The programs carry no version numbers. Entries are dated by when the change reached this repository.

## 2026-11-08: documentation

- `README.md` rewritten as the full reference the book points to: every option of the three programs with its default, every file they write, the columns of each table, how the score is built, troubleshooting, and known issues.
- Added this changelog, `CONTRIBUTING.md` (reporting problems and proposing changes without giving away locations), `CITATION.cff`, and an issue template that asks for no locations.
- `.gitignore` also keeps out the `*_refine_*` folders that the commands in `refine_commands.txt` write to.

## 2026-11-08: latest version from the working session

The programs, `README.md`, `requirements.txt`, `LICENSE` and `.gitignore` from the working session, replacing the copies uploaded on 2026-10-25.

**county_scan.py**
- No per-tile limit on what is kept. `--top`, `--top-offchannel` and `--top-bigwash` now default to 0, meaning all; they were 50, 10 and 10. Checked against one surveyor's list of 229 Mesa County caves, the old limits cut the caves with a candidate within 100 m of an entrance from 124 to 91.
- `--keep-shallow` keeps shallow, wide hollows as their own layer (ids S1, S2 …), limited per tile by `--top-shallow`. Without it they are dropped, as before.
- The off-channel layer is described as holding, in claystone hillsides, many real cave entrances.

**field_kit.py**
- Writes `field_candidates_web.kml`: the same pins and lines without pictures, for Google Earth on the web, which cannot show pictures packed in a KMZ.
- Reads the shallow layer: its own folder in the KMZ, and its own line in the truth report.
- The zoom panel is titled "the hollow and its lip". It no longer claims to show boulders, which bare-earth models smooth away.

**sink_scout.py**
- No change to detection. The self-test prints what it printed before, including the two lines Appendix A quotes.
- Chain wording: chain length and drop describe the wash between the pits, not the cave under it.

**README.md**
- Install from `requirements.txt`.
- `--merge-only` example, and "What is kept" rewritten for the new defaults with the reasons for them.
- `field_candidates_web.kml` described.
- License section.

**New files**
- `requirements.txt`: the six libraries.
- `LICENSE`: CC0 1.0 Universal.
- `.gitignore`: keeps run outputs, elevation tiles and ground-truth files out of the repository.

## 2026-10-25: first upload

`sink_scout.py`, `county_scan.py`, `field_kit.py`, `README.md` and `Appendix_A_Running_Sink_Scout.md`, uploaded through the GitHub web page.
