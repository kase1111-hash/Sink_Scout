# Contributing

Problem reports and changes are welcome. One rule comes before everything else: this repository is public, so nothing posted to it may give away where a pit is. That covers issues, pull requests, commits, comments, test data and file names.

## What gives a location away

- Coordinates of any kind: latitude and longitude, UTM, a `--bbox` or `--check` value.
- Any file the programs write, except the three publishable kinds the README lists under [The rule](README.md#the-rule). That includes KML, KMZ, GPX, GeoJSON and CSV files, shapefiles, `dem.tif` and the other rasters, close-ups, field sheets, and `overview.png` made without `--no-coords`.
- Console output. The ranked table prints latitude and longitude, and the log prints the `--bbox`.
- A tile name with a rank or a run piece (`x71y433_r1c2, rank 3`). The lidar is public and the program gives the same list for the same tile, so anyone can run it and read off the pit.
- A place name tied to a pit, in text or in a file name.

## Reporting a problem

1. Run the self-test and say whether it passes:

       python sink_scout.py --selftest --out selftest

   Its output is safe to post in full: the landscape is synthetic.
2. If you can, show the problem on the self-test, or on the public box Appendix A uses, a city park in Grand Junction with no cave under it: `--bbox 39.075,-108.545,39.085,-108.530`.
3. Otherwise describe it without the location: the command with paths and coordinates replaced by `...`, the error, and the lines of output around it with the `--bbox` line and the latitude and longitude columns deleted.
4. If it shows only on real ground and the files are needed, say so in the issue and ask for a private way to send them. `--save-rasters` writes `dem.tif`, the ground as the program processed it; that is the file to send. Never attach it to an issue.

Give the Python version, the operating system, the versions of the six libraries (`pip freeze` lists them), and where and when you got the programs. The issue template asks for all of this.

If you find a location already in this repository, do not repeat it in an issue. Open one that says only that there is something to remove, and ask for a private channel.

## Proposing a change

- **Run the self-test before and after.** The expected result is in the README. Appendix A of the book prints two of its lines:

      HIT: candidate rank #6 (score 3.127) — depth 15 ft, diam 30 ft, d/w 0.49, catch 10 ac, chain of 1
      HIT: candidate rank #8 (score 0.418) — depth 8 ft, diam 164 ft, d/w 0.05, catch 5 ac, chain of 1, flags: shallow/wide

  A change that alters them breaks the book's instructions for every reader who runs it. If that is the point of the change, say so plainly in the pull request.
- **Keep the README in step.** The book sends readers to the README for every setting and every file the programs write. A new or changed option, default, output file or column belongs in the README in the same pull request, and the change belongs in `CHANGELOG.md`.
- **Keep the install to one line.** The programs need only the six libraries in `requirements.txt`, and Appendix A installs them with one pip command. Ask before adding a library.
- **Keep Python 3.9 working.** It is the oldest version the README promises.
- **Commit no run output.** The `.gitignore` covers the usual folders, not every name you might give one. Check `git status` before committing from a folder you have run in. Test data must be synthetic, like the self-test's.

## License

Contributions are released under CC0 1.0, like the rest of the repository.
