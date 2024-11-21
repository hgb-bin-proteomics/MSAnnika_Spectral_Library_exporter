##### PARAMETERS #####

# name of the mgf file containing the MS2 spectra
SPECTRA_FILE = ["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mgf"]
# you can process multiple files like this:
# SPECTRA_FILE = ["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mgf", "20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_002.mgf"]
# name of the CSM file exported from Proteome Discoverer
CSMS_FILE = "20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.xlsx"
# name of the experiment / run (any descriptive text is allowed)
RUN_NAME = "20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001-(1)"
# name of the sample organism that should be reported in the spectral library
ORGANISM = "Homo sapiens"
# name of the crosslink modification
CROSSLINKER = "DSSO"
# possible modifications and their monoisotopic masses
MODIFICATIONS = \
    {"Oxidation": [15.994915],
     "Carbamidomethyl": [57.021464],
     "DSSO": [54.01056, 85.98264, 103.99320]}
# modifications mapping for xiFDR sequences
MODIFICATIONS_XI = \
    {"Ccm": ["C", "Carbamidomethyl"],
     "Mox": ["M", "Oxidation"]}
# expected ion types (any of a, b, c, x, y, z)
ION_TYPES = ("b", "y")
# maximum expected charge of fragment ions
MAX_CHARGE = 4
# tolerance for matching peaks
MATCH_TOLERANCE = 0.02
# parameters for calculating iRT
iRT_PARAMS = {"iRT_m": 1.3066, "iRT_t": 29.502}
# regex pattern used for parsing scan number from the spectrum title
PARSER_PATTERN = "\\.\\d+\\."
