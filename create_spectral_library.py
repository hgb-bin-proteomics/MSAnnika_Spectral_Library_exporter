#!/usr/bin/env python3

# MS ANNIKA SPECTRAL LIBRARY EXPORTER
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# version tracking
__version = "1.0.0"
__date = "2023-10-10"

############################### WIP ###########################################

# REQUIREMENTS
# pip install pandas
# pip install openpyxl
# pip install tqdm
# pip install pyteomics

##### PARAMETERS #####

SPECTRA_FILE = "20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mgf"
CSMS_FILE = "20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001-(1).xlsx"
RUN_NAME = "20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001-(1)"
CROSSLINKER = "DSSO"
MODIFICATIONS = \
    {"Oxidation": [15.994915],
     "Carbamidomethyl": [57.021464],
     "DSSO": [54.01056, 85.98264, 103.99320]}
ION_TYPES = ("b", "y")
MAX_CHARGE = 4
MATCH_TOLERANCE = 0.02
iRT_PARAMS = {"iRT_m": 1.3066, "iRT_t": 29.502}

######################

"""
NOTES
Column1 -> iterator
linkId -> proteinName1_proteinName2-prot1LinkPos_prot2LinkPos
ProteinID -> proteinName1_proteinName2
StrippedPeptide -> Both peptides without modifications concatenated
FragementGroupId -> peptide1_peptide2-pep1LinkPos_pep2LinkPos:charge
PrecursorCharge -> precursorCharge
PrecursorMz -> precursorMz
ModifiedPeptide -> peptide1(Mod[noXL])_peptide2(Mod[noXL])
IsotopeLabel -> 0
scanID -> scanNr
run -> run
searchID -> MS Annika
crosslinkedResidues -> prot1LinkPos_prot2LinkPos
LabeledSequence -> ModifiedPeptide
iRT -> can maybe calculate
RT -> rt
CCS -> 0
IonMobility -> 0
FragmentCharge -> int
FragmentType -> b/y
FragmentNumber -> int
FragmentPepId -> 0 for first pep, 1 for second pep
FragmentMz -> fragmentMz
RelativeIntensity -> [0..1]
FragmentLossType -> empty
CLContainingFragment -> TRUE / FALSE (XL in frag)
LossyFragment -> FALSE

Every row is a fragment ion?
"""

######################

# import packages
import os
import pandas as pd
from tqdm import tqdm
from pyteomics import mgf, mass

from typing import Dict
from typing import List
from typing import Tuple
from typing import Set
import warnings

# reading spectra
def read_spectra(filename: str) -> Dict[int, Dict]:
    """
    Returns a dictionary that maps scan numbers to spectra:
    Dict[int -> Dict["precursor"        -> float
                     "charge"           -> int
                     "max_intensity"    -> float
                     "peaks"            -> Dict[m/z -> intensity]]
    """

    result_dict = dict()

    with mgf.read(filename) as reader:
        for spectrum in reader:
            scan_nr = int(spectrum["params"]["title"].split("scan=")[1].strip("\""))
            spectrum_dict = dict()
            spectrum_dict["precursor"] = spectrum["params"]["pepmass"]
            spectrum_dict["charge"] = spectrum["params"]["charge"]
            spectrum_dict["max_intensity"] = float(max(spectrum["intensity array"]))
            peaks = dict()
            for i, mz in enumerate(spectrum["m/z array"]):
                peaks[mz] = spectrum["intensity array"][i]
            spectrum_dict["peaks"] = peaks
            result_dict[scan_nr] = spectrum_dict
        reader.close()

    return result_dict

# generate a position to modification mass mapping
def generate_modifications_dict(peptide: str, modification_str: str) -> Dict[int, List[float]]:
    """
    Returns a mapping of peptide positions (0 based) to possible modification masses.
    modification_str is the modification string as returned by MS Annika e.g. K5(DSSO);M7(Oxidation)
    the modification in braces has to be defined in MODIFICATIONS
    """

    modifications_dict = dict()

    modifications = modification_str.split(";")
    for modification in modifications:
        # remove possible white spaces
        modification = modification.strip()
        # get modified amino acid and modification position
        aa_and_pos = modification.split("(")[0]
        # get modification type
        mod = modification.split("(")[1].rstrip(")")

        if aa_and_pos == "Nterm":
            pos = -1
        elif aa_and_pos == "Cterm":
            pos = len(peptide)
        else:
            pos = int(aa_and_pos[1:]) - 1
        if mod in MODIFICATIONS:
            modifications_dict[pos] = MODIFICATIONS[mod]
        else:
            warnings.warn("Modification '" + mod + "' not found!")

    return modifications_dict

# generate all theoretical fragments
# adapted from https://pyteomics.readthedocs.io/en/latest/examples/example_msms.html
def generate_theoretical_fragments(peptide: str, modifications: Dict[int, List[float]], ion_types: Tuple[str] = ("b", "y"), max_charge: int = 1) -> Dict[float, str]:
    """
    Generates a set of theoretical fragment ion masses of the specified peptide with the modifications.
    """

    fragments = dict()

    for i in range(1, len(peptide)):
        for ion_type in ion_types:
            for charge in range(1, max_charge + 1):
                if ion_type[0] in "abc":
                    frag_mass = mass.fast_mass(peptide[:i], ion_type = ion_type, charge = charge)
                    mass_possibilites = set()
                    for mod_pos in modifications.keys():
                        # if the modification is within the fragment:
                        if mod_pos < i:
                            # add the modification mass / charge if its a normal modification
                            if len(modifications[mod_pos]) == 1:
                                frag_mass += modifications[mod_pos][0] / charge
                            else:
                                # if it's a crosslinking modification we add the crosslinker fragment masses
                                # to a set of possible modification mass additions to generate a fragment mass
                                # for every crosslinker fragment
                                for modification in modifications[mod_pos]:
                                    mass_possibilites.add(modification / charge)
                    # we add all possible fragment masses including all crosslinker fragment possibilites
                    if len(mass_possibilites) == 0:
                        if frag_mass not in fragments:
                            fragments[frag_mass] = ion_type + str(i) + "+" + str(charge) + ": " + peptide[:i]
                    else:
                        for mass_possibility in mass_possibilites:
                            if frag_mass + mass_possibility not in fragments:
                                fragments[frag_mass + mass_possibility] = ion_type + str(i) + "+" + str(charge) + ": " + peptide[:i]
                else:
                    frag_mass = mass.fast_mass(peptide[i:], ion_type = ion_type, charge = charge)
                    mass_possibilites = set()
                    for mod_pos in modifications.keys():
                        if mod_pos >= i:
                            if len(modifications[mod_pos]) == 1:
                                frag_mass += modifications[mod_pos][0] / charge
                            else:
                                for modification in modifications[mod_pos]:
                                    mass_possibilites.add(modification / charge)
                    if len(mass_possibilites) == 0:
                        if frag_mass not in fragments:
                            fragments[frag_mass] = ion_type + str(len(peptide) - i) + "+" + str(charge) + ": " + peptide[i:]
                    else:
                        for mass_possibility in mass_possibilites:
                            if frag_mass + mass_possibility not in fragments:
                                fragments[frag_mass + mass_possibility] = ion_type + str(len(peptide) - i) + "+" + str(charge) + ": " + peptide[i:]

    return fragments

def get_matches(row: pd.Series, alpha: bool, spectra: Dict[int, Dict]) -> Tuple[float, Dict[float, str], Dict[float, str]]:

    #todo

    scan_nr = row["First Scan"]

    if alpha:
        sequence = row["Sequence A"]
        modifications = row["Modifications A"]
    else:
        sequence = row["Sequence B"]
        modifications = row["Modifications B"]

    spectrum = spectra[scan_nr]

    modifications_processed = generate_modifications_dict(sequence, modifications)
    theoretical_fragments = generate_theoretical_fragments(sequence, modifications_processed, ion_types = ION_TYPES, max_charge = MAX_CHARGE)

    total_intensity = 0
    matched_fragments = dict()

    for peak_mz in spectrum["peaks"].keys():
        for fragment in theoretical_fragments.keys():
            if round(peak_mz, 4) < round(fragment + MATCH_TOLERANCE, 4) and round(peak_mz, 4) > round(fragment - MATCH_TOLERANCE, 4):
                total_intensity += spectrum["peaks"][peak_mz]
                matched_fragments[peak_mz] = theoretical_fragments[fragment]
                break

    return total_intensity, matched_fragments, theoretical_fragments

def get_positions_in_protein(row: pd.Series) -> Dict[str, int]:

    pep_pos_A = int(row["A in protein"])
    pep_pos_B = int(row["B in protein"])
    xl_pos_A = int(row["Crosslinker Position A"])
    xl_pos_B = int(row["Crosslinker Position B"])

    return {"A": pep_pos_A + xl_pos_A, "B": pep_pos_B + xl_pos_B}

def get_linkId(row: pd.Series) -> str:

    positions = get_positions_in_protein(row)

    return str(row["Accession A"]) + "_" + str(row["Accession B"]) + "-" + str(positions["A"]) + "_" + str(positions["B"])

def get_ProteinID(row: pd.Series) -> str:

    return str(row["Accession A"]) + "_" + str(row["Accession B"])

def get_StrippedPeptide(row: pd.Series) -> str:

    return str(row["Sequence A"]) + str(row["Sequence B"])

def get_FragmentGroupId(row: pd.Series) -> str:

    return str(row["Sequence A"]) + "_" + str(row["Sequence B"]) + "-" + str(row["Crosslinker Position A"]) + "_" + str(row["Crosslinker Position B"]) + ":" + str(row["Charge"])

def get_PrecursorCharge(row: pd.Series) -> int:

    return int(row["Charge"])

def get_PrecursorMz(row: pd.Series) -> float:

    return float(row["m/z [Da]"])

def get_ModifiedPeptide(row: pd.Series) -> str:

    def parse_mod_str(mod_str):
        modifications_dict = dict()
        modifications = mod_str.split(";")
        for modification in modifications:
            aa_and_pos = modification.strip().split("(")[0]
            mod = modification.strip().split("(")[1].rstrip(")")

            if mod == CROSSLINKER:
                continue

            if aa_and_pos == "Nterm":
                pos = 0
            elif aa_and_pos == "Cterm":
                pos = len(peptide)
            else:
                pos = int(aa_and_pos[1:])

            if pos in modifications_dict:
                modifications_dict[pos].append(mod)
            else:
                modifications_dict[pos] = [mod]

        return modifications_dict

    def str_insert(string, index, character):
        return string[:index] + character + string[:index]

    mods_A = parse_mod_str(str(row["Modifications A"]))
    mods_B = parse_mod_str(str(row["Modifications B"]))

    shift = 0
    mod_A_template_str = str(row["Sequence A"])
    for pos in mods_A.keys():
        current_mods = "(" + ", ".join(mods_A[pos]) + ")"
        mod_A_template_str = str_insert(mod_A_template_str, pos + shift, current_mods)
        shift += len(current_mods)

    shift = 0
    mod_B_template_str = str(row["Sequence B"])
    for pos in mods_B.keys():
        current_mods = "(" + ", ".join(mods_B[pos]) + ")"
        mod_B_template_str = str_insert(mod_B_template_str, pos + shift, current_mods)
        shift += len(current_mods)

    return mod_A_template_str + "_" + mod_B_template_str

def get_IsotopeLabel() -> int:

    return 0

def get_scanID(row: pd.Series) -> int:

    return int(row["First Scan"])

def get_run() -> str:

    return RUN_NAME

def get_searchID() -> str:

    return str(row["Crosslink Strategy"])

def get_crosslinkedResidues(row: pd.Series) -> str:

    positions = get_positions_in_protein(row)

    return str(positions["A"]) + "_" + str(positions["B"])

def get_LabeledSequence(row: pd.Series) -> str:

    return get_ModifiedPeptide(row)

def get_iRT(row: pd.Series) -> float:

    return (float(row["RT [min]"]) - iRT_PARAMS["iRT_t"]) / iRT_PARAMS["iRT_m"]

def get_RT(row: pd.Series) -> float:

    return float(row["RT [min]"])

def get_CCS() -> float:

    return 0.0

def get_IonMobility() -> float:

    return 0.0

def get_fragment_values(csm: pd.Series, spectra: Dict) -> Dict:
    # todo
    pass

def main() -> pd.DataFrame:

    print("INFO: Running CSM annotation with input files:\nSpectra: " + SPECTRA_FILE + "\nCSMs: " + CSMS_FILE + "\nDoublets: " + DOUBLETS_FILE)
    print("INFO: Using the following modifications:")
    print(MODIFICATIONS)
    print("INFO: Using the following ion types:")
    print(ION_TYPES)
    print("INFO: Using the following charge states:")
    print([i for i in range(1, MAX_CHARGE + 1)])
    print("INFO: Using a match tolerance of: " + str(MATCH_TOLERANCE) + " Da")
    print("INFO: Starting annotation process...")

    print("INFO: Reading spectra...")
    spectra = read_spectra(SPECTRA_FILE)
    print("INFO: Done reading spectra!")

    print("INFO: Reading CSMs...")
    csms = pd.read_excel(CSMS_FILE)
    print("INFO: Done reading CSMs! Starting fragment ion annotation...")

    tqdm.pandas(desc = "INFO: Progress bar - Annotating alpha peptide fragments")
    csms[["Fragment Intensities A (Sum)", "Matched Ions A", "Theoretical Ions A"]] = \
        csms.progress_apply(lambda row: get_intensities(row, True, spectra),
                            axis = 1,
                            result_type = "expand")
    print("INFO: Done processing alpha peptides!")

    tqdm.pandas(desc = "INFO: Progress bar - Annotating beta peptide fragments")
    csms[["Fragment Intensities B (Sum)", "Matched Ions B", "Theoretical Ions B"]] = \
        csms.progress_apply(lambda row: get_intensities(row, False, spectra),
                            axis = 1,
                            result_type = "expand")
    print("INFO: Done processing beta peptides!")

    csms["Fragment Intensities Total"] = csms.apply(lambda row: get_total_fragment_intensity(row), axis = 1)

    print("INFO: Done annotating fragment ions!")

    if DOUBLETS_FILE is not None and os.path.isfile(DOUBLETS_FILE):
        print("INFO: Doublet file was provided! Reading doublet file...")
        doublets = read_doublets(DOUBLETS_FILE)
        print("INFO: Done reading doublet file! Starting doublet annotation...")
        tqdm.pandas(desc = "INFO: Progress bar - Annotating doublets")
        csms[["Alpha Light Intensities (Sum)", "Alpha Heavy Intensities (Sum)", "Beta Light Intensities (Sum)", "Beta Heavy Intensities (Sum)",
              "Alpha Light Peaks", "Alpha Heavy Peaks", "Beta Light Peaks", "Beta Heavy Peaks"]] = \
              csms.progress_apply(lambda row: get_doublets(row, spectra, doublets),
                                  axis = 1,
                                  result_type = "expand")
        print("INFO: Done annotating doublets! Calculating intensities...")
        csms["Alpha Doublet Intensities Total"] = csms.apply(lambda row: row["Alpha Light Intensities (Sum)"] + row["Alpha Heavy Intensities (Sum)"], axis = 1)
        csms["Beta Doublet Intensities Total"] = csms.apply(lambda row: row["Beta Light Intensities (Sum)"] + row["Beta Heavy Intensities (Sum)"], axis = 1)
        csms["Doublet Intensities Total"] = csms.apply(lambda row: get_total_doublet_intensity(row), axis = 1)
        print("INFO: Done calculating doublet intensities!")
    else:
        print("INFO: Doublet file was not provided or not found! Skipping doublet annotation.")

    print("INFO: Calculating total spectrum intensities...")
    tqdm.pandas(desc = "INFO: Progress bar - Intensity per spectrum")
    csms["Total Intensity in Spectrum"] = csms.progress_apply(lambda row: get_spectrum_intensity(row, spectra), axis = 1)
    print("INFO: Done calculating total spectrum intensities!")

    csms.to_excel(".".join(CSMS_FILE.split(".")[:-1]) + "_with_intensities.xlsx")
    print("SUCCESS: Output file generated as '" + ".".join(CSMS_FILE.split(".")[:-1]) + "_with_intensities.xlsx" + "'!")

    return csms

if __name__ == "__main__":

        main()
