##### PARAMETERS #####

# name of the mgf file containing the MS2 spectra
SPECTRA_FILE = ["XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"]
# name of the CSM file exported from Proteome Discoverer
CSMS_FILE = "example_CSM_xiFDR2.2.1.csv"
# name of the experiment / run (any descriptive text is allowed)
RUN_NAME = "XLpeplib_Beveridge_DSS_rep1_xiSearch"
# name of the sample organism that should be reported in the spectral library
ORGANISM = "S. pyogenes"
# name of the crosslink modification
CROSSLINKER = "BS3"
# possible modifications and their monoisotopic masses
MODIFICATIONS = \
    {"Oxidation": [15.994915],
     "Carbamidomethyl": [57.021464],
     "BS3": [138.06808]}
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
# only take the best CSM per unique peptidoform and charge (True) or not (False)
GROUP_PRECURSORS = True
# raise an error if spectra do not contain FAIMS compensation voltage information
ERROR_ON_NO_FAIMS = True
