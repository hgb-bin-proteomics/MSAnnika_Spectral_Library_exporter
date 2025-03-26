# Post Processing

After running a search with [Spectronaut](https://biognosys.com/software/spectronaut/)
using the created spectral library you will want to validate your results. For
this purpose we have created a post processing script that uses the Spectronaut
result file and the spectral library and creates a result file that can then be
used with [xiFDR](https://www.rappsilberlab.org/software/xifdr/) for validation.

## Usage of the post processing script

- First you need to set a few parameters in the `post_process.py` script:
  ```python
  CROSSLINKER = "PhoX" # name of the crosslinker
  CROSSLINKER_MASS = 209.97181 # delta mass of the crosslinker
  SPECTRONAUT_DELIM = "," # delimiter in Spectronaut output file, e.g. "," for comma delimited files, "\t" for tab delimited files
  SPECTRONAUT_MATCH_TOLERANCE = 0.05 # match tolerance in Da
  SPECTRONAUT_FRAGMENT_MZ_COLUMN_NAME = "F.CalibratedMz" # which F Mz to use for matching
  SPECTRONAUT_CSCORE_COLUMN_NAME = "EG.Cscore" # which Cscore to use for re-soring
  ```
- Make sure that `post_process.py`, the Spectronaut result file, and the spectral
  library are in the same directory.
- Launch a shell/terminal in that directory.
- Execute the command `python post_process.py SPECTRONAUT_RESULT.csv`.
  - You should replace `SPECTRONAUT_RESULT.csv` with the filename of your
    Spectronaut result file.
- You will be asked to make sure that you set the correct parameters in the
  script. If you did you can accept with `yes`.
- Then the script should run through and give you information about the annotation.
- When the script is finished it will create three new files:
  - `*_annotated.csv`: The Spectronaut result file with additionally calculated scores.
  - `*_grouped_by_residue_pair.csv`: Additionally to the previous file, the results
    are grouped by residue pair to create "pseudo CSMs". This file can be used in
    xiFDR, but you will most likely need to manually map a few columns in the xiFDR
    GUI.
  - `*_xiFDR.csv`: Like the last file but this time only the necessary columns for
    xiFDR are included. There should not be any manual mapping in the xiFDR GUI
    necessary, except to choose which score to use.

> [!Note]
> The spectral library file name is directly parsed from the Spectronaut result
> file from the column `EG.Library`. Make sure that this matches up with the
> filename of the spectral library in the current directory!
