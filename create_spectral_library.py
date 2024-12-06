#!/usr/bin/env python3

# MS ANNIKA SPECTRAL LIBRARY EXPORTER
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# version tracking
__version = "1.4.3"
__date = "2024-12-06"

# REQUIREMENTS
# pip install pandas
# pip install openpyxl
# pip install pyteomics

##### PARAMETERS #####

from config import SPECTRA_FILE
from config import CSMS_FILE
from config import RUN_NAME
from config import CROSSLINKER
from config import MODIFICATIONS
from config import MODIFICATIONS_XI
from config import ION_TYPES
from config import MAX_CHARGE
from config import MATCH_TOLERANCE
from config import iRT_PARAMS
from config import ORGANISM
from config import PARSER_PATTERN

######################

# import packages
import re
import pandas as pd
from pyteomics import mgf, mass

from typing import Dict
from typing import List
from typing import Tuple
from typing import Set
from typing import Union
from typing import BinaryIO
from typing import Any
import warnings

##################### FILE READERS #####################

def parse_xi(result_file: str, spectra: Dict[str, Any]) -> pd.DataFrame:
    """Parses the xiFDR CSM result file and returns it in MS Annika format for
    spectral library creation.
    """
    xi = pd.read_csv(result_file)
    ## needed cols
    # Sequence A
    # Modifications A
    # Sequence B
    # Modifications B
    # First Scan
    # Spectrum File
    # A in protein
    # B in protein
    # Crosslinker Position A
    # Crosslinker Position B
    # Accession A
    # Accession B
    # Charge
    # m/z [Da]
    # Crosslink Strategy
    # RT [min]
    # Compensation Voltage
    ms_annika_struc = {"Sequence A": [],
                       "Modifications A": [],
                       "Sequence B": [],
                       "Modifications B": [],
                       "First Scan": [],
                       "Spectrum File": [],
                       "A in protein": [],
                       "B in protein": [],
                       "Crosslinker Position A": [],
                       "Crosslinker Position B": [],
                       "Accession A": [],
                       "Accession B": [],
                       "Charge": [],
                       "m/z [Da]": [],
                       "Crosslink Strategy": [],
                       "RT [min]": [],
                       "Compensation Voltage": []}

    # parsing functions
    def xi_get_sequence(row: pd.Series, alpha: bool = True) -> str:
        seq = str(row["PepSeq1"]).strip() if alpha else str(row["PepSeq2"]).strip()
        seq_a = ""
        for aa in seq:
            if aa.isupper():
                seq_a += aa
        return seq_a

    def xi_get_modifications(row: pd.Series, alpha: bool = True) -> str:
        seq = str(row["PepSeq1"]).strip() if alpha else str(row["PepSeq2"]).strip()
        clean_seq = xi_get_sequence(row, alpha)
        xl_pos = int(row["LinkPos1"]) if alpha else int(row["LinkPos2"])

        if len(MODIFICATIONS_XI) > 10:
            msg = "Found more than 10 possible modifications for xi. " + \
                  "Maximum number of modifications supported is 10. " + \
                  "Please update MODIFICATIONS_XI in the config file!"
            raise RuntimeError(msg)

        mod_map = dict()
        mod_map_rev = dict()
        for i, key in enumerate(MODIFICATIONS_XI.keys()):
            mod_map[str(i)] = key
            mod_map_rev[key] = str(i)

        for mod in MODIFICATIONS_XI.keys():
            seq = seq.replace(mod, mod_map_rev[mod])

        mod_str = ""
        for i, aa in enumerate(seq):
            if aa in mod_map:
                mod_str += f"{MODIFICATIONS_XI[mod_map[aa]][0]}{i+1}({MODIFICATIONS_XI[mod_map[aa]][1]});"

        mod_str += f"{clean_seq[xl_pos-1]}{xl_pos}({str(row['Crosslinker']).strip()})"

        return mod_str

    def xi_get_rt(row: pd.Series, spectra: Dict[str, Any]) -> float:
        spec_file_name = ".".join(str(row["PeakListFileName"]).split(".")[:-1]).strip()
        rt = spectra[spec_file_name][int(row["scan"])]["rt"]
        return rt / 60.0

    def xi_get_cv(row: pd.Series, spectra: Dict[str, Any]) -> float:
        # I don't think we get this from the MGF file?
        return 0.0

    for i, row in xi.iterrows():
        if row["isDecoy"]:
            continue
        ms_annika_struc["Sequence A"].append(xi_get_sequence(row, True))
        ms_annika_struc["Sequence B"].append(xi_get_sequence(row, False))
        ms_annika_struc["Modifications A"].append(xi_get_modifications(row, True))
        ms_annika_struc["Modifications B"].append(xi_get_modifications(row, False))
        ms_annika_struc["First Scan"].append(int(row["scan"]))
        ms_annika_struc["Spectrum File"].append(str(row["PeakListFileName"]).strip())
        ms_annika_struc["A in protein"].append(int(row["PepPos1"])-1)
        ms_annika_struc["B in protein"].append(int(row["PepPos2"])-1)
        ms_annika_struc["Crosslinker Position A"].append(int(row["LinkPos1"]))
        ms_annika_struc["Crosslinker Position B"].append(int(row["LinkPos2"]))
        ms_annika_struc["Accession A"].append(str(row["Protein1"]).strip())
        ms_annika_struc["Accession B"].append(str(row["Protein2"]).strip())
        ms_annika_struc["Charge"].append(int(row["exp charge"]))
        ms_annika_struc["m/z [Da]"].append(float(row["exp m/z"]))
        ms_annika_struc["Crosslink Strategy"].append("xi")
        ms_annika_struc["RT [min]"].append(xi_get_rt(row, spectra))
        ms_annika_struc["Compensation Voltage"].append(xi_get_cv(row, spectra))

    return pd.DataFrame(ms_annika_struc)

############### SPECTRAL LIBRARY CREATOR ###############

# parse scan number from pyteomics mgf params
def parse_scannr(params: Dict, pattern: str = PARSER_PATTERN, i: int = 0) -> Tuple[int, int]:
    """Parses the scan number from the params dictionary of the pyteomics mgf
    spectrum.

    Parameters
    ----------
    params : Dict
        The "params" dictionary of the pyteomics mgf spectrum.

    pattern : str
        Regex pattern to use for parsing the scan number from the title if it
        can't be infered otherwise.

    i : int
        The scan number to be returned in case of failure.

    Returns
    -------
    (exit_code, scan_nr) : Tuple
        A tuple with the exit code (0 if successful, 1 if parsing failed) at the
        first position [0] and the scan number at the second position [1].
    """

    # prefer scans attr over title attr
    if "scans" in params:
        try:
            return (0, int(params["scans"]))
        except:
            pass

    # try parse title
    if "title" in params:

        # if there is a scan token in the title, try parse scan_nr
        if "scan" in params["title"]:
            try:
                return (0, int(params["title"].split("scan=")[1].strip("\"")))
            except:
                pass

        # else try to parse by pattern
        try:
            scan_nr = re.findall(pattern, params["title"])[0]
            scan_nr = re.sub(r"[^0-9]", "", scan_nr)
            if len(scan_nr) > 0:
                return (0, int(scan_nr))
        except:
            pass

        # else try parse whole title
        try:
            return (0, int(params["title"]))
        except:
            pass

    # return unsuccessful parse
    return (1, i)

# reading spectra
# def read_spectra(filename: str | BinaryIO) -> Dict[int, Dict]:
# for backward compatibility >>
def read_spectra(filename: Union[str, BinaryIO]) -> Dict[int, Dict]:
    """
    Returns a dictionary that maps scan numbers to spectra:
    Dict[int -> Dict["precursor"        -> float
                     "charge"           -> int
                     "rt"               -> float
                     "max_intensity"    -> float
                     "peaks"            -> Dict[m/z -> intensity]]
    """

    result_dict = dict()

    with mgf.read(filename) as reader:
        for spectrum in reader:
            parser_result = parse_scannr(spectrum["params"])
            if parser_result[0] != 0:
                raise RuntimeError(f"Could not parse scan number for spectrum {spectrum}. Please adjust PARSER_PATTERN in the config file!")
            scan_nr = parser_result[1]
            spectrum_dict = dict()
            spectrum_dict["precursor"] = spectrum["params"]["pepmass"]
            spectrum_dict["charge"] = spectrum["params"]["charge"]
            spectrum_dict["rt"] = spectrum["params"]["rtinseconds"] if "rtinseconds" in spectrum["params"] else 0.0
            spectrum_dict["max_intensity"] = float(max(spectrum["intensity array"])) if len(spectrum["intensity array"]) > 0 else 0.0
            peaks = dict()
            for i, mz in enumerate(spectrum["m/z array"]):
                peaks[mz] = spectrum["intensity array"][i]
            spectrum_dict["peaks"] = peaks
            result_dict[scan_nr] = spectrum_dict
        reader.close()

    return result_dict

# read multiple spectra files
def read_multiple_spectra(filenames: List[str]) -> Dict[str, Dict[int, Dict]]:
    """
    Returns a dictionary that maps filenames to scan numbers to spectra:
    Dict[str -> Dict[int -> Dict["precursor"        -> float
                                 "charge"           -> int
                                 "max_intensity"    -> float
                                 "peaks"            -> Dict[m/z -> intensity]]
    """

    result_dict = dict()

    for filename in filenames:
        current_spectra_file = ".".join(filename.split(".")[:-1]).strip()
        result_dict[current_spectra_file] = read_spectra(filename)

    return result_dict

# read multiple spectra files - streamlit version
def read_multiple_spectra_streamlit(st_files) -> Dict[str, Dict[int, Dict]]:
    """
    Returns a dictionary that maps filenames to scan numbers to spectra:
    Dict[str -> Dict[int -> Dict["precursor"        -> float
                                 "charge"           -> int
                                 "max_intensity"    -> float
                                 "peaks"            -> Dict[m/z -> intensity]]
    """

    result_dict = dict()

    for st_file in st_files:
        current_spectra_file = ".".join(st_file.name.split(".")[:-1]).strip()
        result_dict[current_spectra_file] = read_spectra(st_file)

    return result_dict

# generate a position to modification mass mapping
def generate_modifications_dict(peptide: str,
                                modification_str: str,
                                possible_modifications: Dict[str, List[float]] = MODIFICATIONS) -> Dict[int, List[float]]:
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
        if mod in possible_modifications:
            modifications_dict[pos] = possible_modifications[mod]
        else:
            warnings.warn("Modification '" + mod + "' not found!")

    return modifications_dict

# generate all theoretical fragments
# adapted from https://pyteomics.readthedocs.io/en/latest/examples/example_msms.html
def generate_theoretical_fragments(peptide: str,
                                   modifications: Dict[int, List[float]],
                                   ion_types: Tuple[str] = ("b", "y"),
                                   max_charge: int = 1) -> Dict[float, str]:
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

# get all fragments and their annotations
def get_fragments(row: pd.Series,
                  alpha: bool,
                  spectra: Dict[str, Dict[int, Dict]],
                  crosslinker: str = CROSSLINKER,
                  possible_modifications: Dict[str, List[float]] = MODIFICATIONS,
                  ion_types: Tuple[str] = ION_TYPES,
                  max_charge: int = MAX_CHARGE,
                  match_tolerance: float = MATCH_TOLERANCE) -> List[Dict]:
    """
    Generates all fragments with the necessary spectral library annotations for a given CSM peptide.
    """

    # function to check if the fragment contains the crosslinker
    def check_if_xl_in_frag(row, alpha, ion_type, fragment, crosslinker):

        if alpha:
            peptide = row["Sequence A"]
            mods = row["Modifications A"]
        else:
            peptide = row["Sequence B"]
            mods = row["Modifications B"]

        pos = 0
        mods_list = mods.split(";")
        for mod_in_list in mods_list:
            aa_and_pos = mod_in_list.strip().split("(")[0]
            mod = mod_in_list.strip().split("(")[1].rstrip(")")

            if mod == crosslinker:
                if aa_and_pos == "Nterm":
                    pos = -1
                elif aa_and_pos == "Cterm":
                    pos = len(peptide)
                else:
                    pos = int(aa_and_pos[1:]) - 1
                break

        if ion_type in "abc":
            if len(fragment) > pos:
                return True
        else:
            if len(peptide) - len(fragment) <= pos:
                return True

        return False
    # end function

    fragments = list()

    scan_nr = row["First Scan"]

    if alpha:
        sequence = row["Sequence A"]
        modifications = row["Modifications A"]
    else:
        sequence = row["Sequence B"]
        modifications = row["Modifications B"]

    current_spectra_file = ".".join(row["Spectrum File"].split(".")[:-1]).strip()
    spectrum = spectra[current_spectra_file][scan_nr]

    modifications_processed = generate_modifications_dict(sequence, modifications, possible_modifications)
    theoretical_fragments = generate_theoretical_fragments(sequence, modifications_processed, ion_types, max_charge)

    matched_fragments = dict()

    # match fragments
    for peak_mz in spectrum["peaks"].keys():
        for fragment in theoretical_fragments.keys():
            if round(peak_mz, 4) < round(fragment + match_tolerance, 4) and round(peak_mz, 4) > round(fragment - match_tolerance, 4):
                matched_fragments[peak_mz] = theoretical_fragments[fragment]
                break

    # get annotations
    for match in matched_fragments.keys():
        fragment_charge = int(matched_fragments[match].split("+")[1].split(":")[0])
        fragment_type = str(matched_fragments[match][0])
        fragment_number = int(matched_fragments[match].split("+")[0][1:])
        fragment_pep_id = 0 if alpha else 1
        fragment_mz = match
        fragment_rel_intensity = float(spectrum["peaks"][match] / spectrum["max_intensity"])
        fragment_loss_type = ""
        fragment_contains_xl = check_if_xl_in_frag(row, alpha, fragment_type, matched_fragments[match].split(":")[1].strip(), crosslinker)
        fragment_lossy = False
        fragments.append({"FragmentCharge": fragment_charge,
                          "FragmentType": fragment_type,
                          "FragmentNumber": fragment_number,
                          "FragmentPepId": fragment_pep_id,
                          "FragmentMz": fragment_mz,
                          "RelativeIntensity": fragment_rel_intensity,
                          "FragmentLossType": fragment_loss_type,
                          "CLContainingFragment": fragment_contains_xl,
                          "LossyFragment": fragment_lossy
                          })

    return fragments

# get crosslink position in proteins
def get_positions_in_protein(row: pd.Series) -> Dict[str, int]:
    """
    Returns the crosslink position of the first protein of peptide alpha and the first protein of peptide beta.
    """

    pep_pos_A = int(row["A in protein"]) if ";" not in str(row["A in protein"]) else int(row["A in protein"].split(";")[0])
    pep_pos_B = int(row["B in protein"]) if ";" not in str(row["B in protein"]) else int(row["B in protein"].split(";")[0])
    xl_pos_A = int(row["Crosslinker Position A"])
    xl_pos_B = int(row["Crosslinker Position B"])

    return {"A": pep_pos_A + xl_pos_A, "B": pep_pos_B + xl_pos_B}

##### DECOY GENERATION #####
# decoy generation implemented as described by Zhang et al. here: https://doi.org/10.1021/acs.jproteome.7b00614

def generate_decoy_csm_dd(row: pd.Series, crosslinker: str = CROSSLINKER) -> pd.Series:
    """
    """

    decoy_csm = row.copy(deep = True)

    # decoy seq
    seq_a = str(decoy_csm["Sequence A"]).strip()
    decoy_seq_a = seq_a[:-1][::-1] + seq_a[-1]
    seq_b = str(decoy_csm["Sequence B"]).strip()
    decoy_seq_b = seq_b[:-1][::-1] + seq_b[-1]
    decoy_csm["Sequence A"] = decoy_seq_a
    decoy_csm["Sequence B"] = decoy_seq_b

    # decoy mods
    ## <calculate_new_position>
    def calculate_new_position(csm, alpha, modification):
        sequence = csm["Sequence B"]
        if alpha:
            sequence = csm["Sequence A"]
        if "Nterm" in modification:
            return modification
        if "Cterm" in modification:
            return modification
        pos = int(modification.split("(")[0][1:]) - 1
        aa = modification.split("(")[0][0].strip()
        ptm = modification.split("(")[1].split(")")[0].strip()
        if pos == (len(sequence) - 1):
            return modification
        new_pos = len(sequence) - 2 - pos
        if aa != sequence[new_pos]:
            warnings.warn(f"Target and decoy modification positions may to match (decoy position may be incorrect)!", UserWarning)
        #assert aa == sequence[new_pos]
        if ptm == crosslinker:
            if alpha:
                csm["Crosslinker Position A"] = new_pos + 1
            else:
                csm["Crosslinker Position B"] = new_pos + 1

        return f"{aa}{new_pos + 1}({ptm})"
    ## </calculate_new_position>

    decoy_mods_a = [calculate_new_position(decoy_csm, True, mod.strip()) for mod in str(decoy_csm["Modifications A"]).split(";")]
    decoy_mods_b = [calculate_new_position(decoy_csm, False, mod.strip()) for mod in str(decoy_csm["Modifications B"]).split(";")]
    decoy_csm["Modifications A"] = ";".join(decoy_mods_a)
    decoy_csm["Modifications B"] = ";".join(decoy_mods_b)

    return decoy_csm

def generate_decoy_csm_td(row: pd.Series, crosslinker: str = CROSSLINKER) -> pd.Series:
    """
    """

    decoy_csm = row.copy(deep = True)

    # decoy seq
    seq_b = str(decoy_csm["Sequence B"]).strip()
    decoy_seq_b = seq_b[:-1][::-1] + seq_b[-1]
    decoy_csm["Sequence B"] = decoy_seq_b

    # decoy mods
    ## <calculate_new_position>
    def calculate_new_position(csm, alpha, modification):
        sequence = csm["Sequence B"]
        if alpha:
            sequence = csm["Sequence A"]
        if "Nterm" in modification:
            return modification
        if "Cterm" in modification:
            return modification
        pos = int(modification.split("(")[0][1:]) - 1
        aa = modification.split("(")[0][0].strip()
        ptm = modification.split("(")[1].split(")")[0].strip()
        if pos == (len(sequence) - 1):
            return modification
        new_pos = len(sequence) - 2 - pos
        if aa != sequence[new_pos]:
            warnings.warn(f"Target and decoy modification positions may to match (decoy position may be incorrect)!", UserWarning)
        #assert aa == sequence[new_pos]
        if ptm == crosslinker:
            if alpha:
                csm["Crosslinker Position A"] = new_pos + 1
            else:
                csm["Crosslinker Position B"] = new_pos + 1

        return f"{aa}{new_pos + 1}({ptm})"
    ## </calculate_new_position>

    decoy_mods_b = [calculate_new_position(decoy_csm, False, mod.strip()) for mod in str(decoy_csm["Modifications B"]).split(";")]
    decoy_csm["Modifications B"] = ";".join(decoy_mods_b)

    return decoy_csm

def generate_decoy_csm_dt(row: pd.Series, crosslinker: str = CROSSLINKER) -> pd.Series:
    """
    """

    decoy_csm = row.copy(deep = True)

    # decoy seq
    seq_a = str(decoy_csm["Sequence A"]).strip()
    decoy_seq_a = seq_a[:-1][::-1] + seq_a[-1]
    decoy_csm["Sequence A"] = decoy_seq_a

    # decoy mods
    ## <calculate_new_position>
    def calculate_new_position(csm, alpha, modification):
        sequence = csm["Sequence B"]
        if alpha:
            sequence = csm["Sequence A"]
        if "Nterm" in modification:
            return modification
        if "Cterm" in modification:
            return modification
        pos = int(modification.split("(")[0][1:]) - 1
        aa = modification.split("(")[0][0].strip()
        ptm = modification.split("(")[1].split(")")[0].strip()
        if pos == (len(sequence) - 1):
            return modification
        new_pos = len(sequence) - 2 - pos
        if aa != sequence[new_pos]:
            warnings.warn(f"Target and decoy modification positions may to match (decoy position may be incorrect)!", UserWarning)
        #assert aa == sequence[new_pos]
        if ptm == crosslinker:
            if alpha:
                csm["Crosslinker Position A"] = new_pos + 1
            else:
                csm["Crosslinker Position B"] = new_pos + 1

        return f"{aa}{new_pos + 1}({ptm})"
    ## </calculate_new_position>

    decoy_mods_a = [calculate_new_position(decoy_csm, True, mod.strip()) for mod in str(decoy_csm["Modifications A"]).split(";")]
    decoy_csm["Modifications A"] = ";".join(decoy_mods_a)

    return decoy_csm

def get_decoy_fragments(decoy_csm: pd.Series,
                        target_fragments: List[Dict],
                        possible_modifications: Dict[str, List[float]] = MODIFICATIONS,
                        crosslinker: str = CROSSLINKER) -> List[Dict]:
    """
    """

    decoy_fragments = list()

    ## <get_decoy_mzs>
    def get_decoy_mzs(decoy_csm, pep_id, ion_type, ion_number, charge, possible_modifications) -> List[float]:
        decoy_mzs = list()
        seq = decoy_csm["Sequence A"]
        mods_str = decoy_csm["Modifications A"]
        if pep_id == 1:
            seq = decoy_csm["Sequence B"]
            mods_str = decoy_csm["Modifications B"]
        mods = generate_modifications_dict(seq, mods_str, possible_modifications)
        if ion_type in "abc":
            end_pos = ion_number
            frag_mass = mass.fast_mass(seq[:end_pos], ion_type = ion_type, charge = charge)
            mz_possibilites = set()
            for mod_pos in mods.keys():
                # if the modification is within the fragment:
                if mod_pos < end_pos:
                    # add the modification mass / charge if its a normal modification
                    if len(mods[mod_pos]) == 1:
                        frag_mass += mods[mod_pos][0] / charge
                    else:
                        # if it's a crosslinking modification we add the crosslinker fragment masses
                        # to a set of possible modification mass additions to generate a fragment mass
                        # for every crosslinker fragment
                        for mod in mods[mod_pos]:
                            mz_possibilites.add(mod / charge)
            # we add all possible fragment masses including all crosslinker fragment possibilites
            if len(mz_possibilites) == 0:
                if frag_mass not in decoy_mzs:
                    decoy_mzs.append(frag_mass)
            else:
                for mz_possibility in mz_possibilites:
                    if frag_mass + mz_possibility not in decoy_mzs:
                        decoy_mzs.append(frag_mass + mz_possibility)
            return decoy_mzs
        else:
            start_pos = len(seq) - ion_number
            frag_mass = mass.fast_mass(seq[start_pos:], ion_type = ion_type, charge = charge)
            mz_possibilites = set()
            for mod_pos in mods.keys():
                if mod_pos >= start_pos:
                    if len(mods[mod_pos]) == 1:
                        frag_mass += mods[mod_pos][0] / charge
                    else:
                        for mod in mods[mod_pos]:
                            mz_possibilites.add(mod / charge)
            if len(mz_possibilites) == 0:
                if frag_mass not in decoy_mzs:
                    decoy_mzs.append(frag_mass)
            else:
                for mz_possibility in mz_possibilites:
                    if frag_mass + mz_possibility not in decoy_mzs:
                        decoy_mzs.append(frag_mass + mz_possibility)
            return decoy_mzs
    ## </get_decoy_mzs>

    ## <check_if_xl_in_frag>
    def check_if_xl_in_frag(decoy_csm, pep_id, ion_type, ion_number, crosslinker) -> bool:
        if pep_id == 0:
            peptide = decoy_csm["Sequence A"]
            mods = decoy_csm["Modifications A"]
        else:
            peptide = decoy_csm["Sequence B"]
            mods = decoy_csm["Modifications B"]

        pos = 0
        mods_list = mods.split(";")
        for mod_in_list in mods_list:
            aa_and_pos = mod_in_list.strip().split("(")[0]
            mod = mod_in_list.strip().split("(")[1].rstrip(")")

            if mod == crosslinker:
                if aa_and_pos == "Nterm":
                    pos = -1
                elif aa_and_pos == "Cterm":
                    pos = len(peptide)
                else:
                    pos = int(aa_and_pos[1:]) - 1
                break

        if ion_type in "abc":
            if ion_number > pos:
                return True
        else:
            if len(peptide) - ion_number <= pos:
                return True

        return False
    ## </check_if_xl_in_frag>

    for fragment in target_fragments:
        fragment_charge = fragment["FragmentCharge"]
        fragment_type = fragment["FragmentType"]
        fragment_number = fragment["FragmentNumber"]
        fragment_pep_id = fragment["FragmentPepId"]
        fragment_mzs = get_decoy_mzs(decoy_csm, fragment_pep_id, fragment_type, fragment_number, fragment_charge, possible_modifications)
        fragment_rel_intensity = fragment["RelativeIntensity"]
        fragment_loss_type = ""
        fragment_contains_xl = check_if_xl_in_frag(decoy_csm, fragment_pep_id, fragment_type, fragment_number, crosslinker)
        fragment_lossy = False
        for fragment_mz in fragment_mzs:
            decoy_fragments.append({"FragmentCharge": fragment_charge,
                                    "FragmentType": fragment_type,
                                    "FragmentNumber": fragment_number,
                                    "FragmentPepId": fragment_pep_id,
                                    "FragmentMz": fragment_mz,
                                    "RelativeIntensity": fragment_rel_intensity,
                                    "FragmentLossType": fragment_loss_type,
                                    "CLContainingFragment": fragment_contains_xl,
                                    "LossyFragment": fragment_lossy
                                  })

    return decoy_fragments

##### SPECTRAL LIBRARY COLUMNS #####

# get the linkId value
def get_linkId(row: pd.Series) -> str:
    """
    Returns the first accession of the alpha peptide and the first accession of the beta peptide + the corresponding crosslink positions.
    """

    positions = get_positions_in_protein(row)
    accession_a = row["Accession A"] if ";" not in row["Accession A"] else row["Accession A"].split(";")[0]
    accession_b = row["Accession B"] if ";" not in row["Accession B"] else row["Accession B"].split(";")[0]

    return str(accession_a) + "_" + str(accession_b) + "-" + str(positions["A"]) + "_" + str(positions["B"])

# get the ProteinID value
def get_ProteinID(row: pd.Series) -> str:
    """
    Returns the first accession of the alpha peptide and the first accession of the beta peptide.
    """

    accession_a = row["Accession A"] if ";" not in row["Accession A"] else row["Accession A"].split(";")[0]
    accession_b = row["Accession B"] if ";" not in row["Accession B"] else row["Accession B"].split(";")[0]

    return str(accession_a) + "_" + str(accession_b)

def get_Organism() -> str:
    """
    Returns the organism of the sample.
    """

    return str(ORGANISM)

# get the StrippedPeptide value
def get_StrippedPeptide(row: pd.Series) -> str:
    """
    Returns the sequences of the cross-linked peptides concatenated.
    """

    return str(row["Sequence A"]) + str(row["Sequence B"])

# get the FragmentGroupId value
def get_FragmentGroupId(row: pd.Series) -> str:
    """
    Returns 'SequenceA_SequenceB-CrosslinkerPositionA_CrosslinkerPositionB:Charge'.
    """

    return str(row["Sequence A"]) + "_" + str(row["Sequence B"]) + "-" + str(row["Crosslinker Position A"]) + "_" + str(row["Crosslinker Position B"]) + ":" + str(row["Charge"])

# get the PrecursorCharge value
def get_PrecursorCharge(row: pd.Series) -> int:
    """
    Returns the precursor charge.
    """

    return int(row["Charge"])

# get the PrecursorMz value
def get_PrecursorMz(row: pd.Series) -> float:
    """
    Returns the precursor m/z.
    """

    return float(row["m/z [Da]"])

# get the ModifiedPeptide value
def get_ModifiedPeptide(row: pd.Series,
                        crosslinker: str = CROSSLINKER) -> str:
    """
    Returns SequenceA with modification annotations (without crosslinker) _ SequenceB with modification annotations (without crosslinker).
    """

    # helper function to parse MS Annika modification string
    def parse_mod_str(peptide, mod_str, crosslinker):
        modifications_dict = dict()
        modifications = mod_str.split(";")
        for modification in modifications:
            aa_and_pos = modification.strip().split("(")[0]
            mod = modification.strip().split("(")[1].rstrip(")")

            if mod == crosslinker:
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
    # end function

    # helper function to insert string into string
    def str_insert(string, index, character):
        return string[:index] + character + string[index:]
    # end function

    mods_A = parse_mod_str(str(row["Sequence A"]), str(row["Modifications A"]), crosslinker)
    mods_B = parse_mod_str(str(row["Sequence B"]), str(row["Modifications B"]), crosslinker)

    # generate annotation for sequence A
    shift = 0
    mod_A_template_str = str(row["Sequence A"])
    for pos in mods_A.keys():
        current_mods = "[" + ", ".join(mods_A[pos]) + "]"
        mod_A_template_str = str_insert(mod_A_template_str, pos + shift, current_mods)
        shift += len(current_mods)

    # generate annotation for sequence B
    shift = 0
    mod_B_template_str = str(row["Sequence B"])
    for pos in mods_B.keys():
        current_mods = "[" + ", ".join(mods_B[pos]) + "]"
        mod_B_template_str = str_insert(mod_B_template_str, pos + shift, current_mods)
        shift += len(current_mods)

    return mod_A_template_str + "_" + mod_B_template_str

# get the IsotopeLabel value
def get_IsotopeLabel() -> int:
    """
    Dummy function.
    """
    return 0

# get the file value
def get_filename(row: pd.Series) -> str:
    """
    Returns the filename of the corresponding RAW/MGF of the CSM.
    """

    return str(row["Spectrum File"])

# get the scanID value
def get_scanID(row: pd.Series) -> int:
    """
    Returns the scan nr. of the CSM.
    """

    return int(row["First Scan"])

# get the run value
def get_run(run_name: str = RUN_NAME) -> str:
    """
    Returns the run name specified in config.py
    """

    return run_name

# get the searchID value
def get_searchID(row: pd.Series) -> str:
    """
    Returns the identifying search engine name.
    """

    return str(row["Crosslink Strategy"])

# get the crosslinkedResidues value
def get_crosslinkedResidues(row: pd.Series) -> str:
    """
    Returns the positions of the cross-linked residues of the first proteins of the cross-linked peptides respectively, seperated by '_'.
    """

    positions = get_positions_in_protein(row)

    return str(positions["A"]) + "_" + str(positions["B"])

# get the LabeledSequence value
def get_LabeledSequence(row: pd.Series,
                        crosslinker: str = CROSSLINKER) -> str:
    """
    Returns SequenceA with modification annotations (without crosslinker) _ SequenceB with modification annotations (without crosslinker).
    """

    return get_ModifiedPeptide(row, crosslinker)

# get the iRT value
def get_iRT(row: pd.Series,
            iRT_t: float = iRT_PARAMS["iRT_t"],
            iRT_m: float = iRT_PARAMS["iRT_m"]) -> float:
    """
    Returns the calculated iRT using the values specified in config.py.
    """

    return (float(row["RT [min]"]) - iRT_t) / iRT_m

# get the RT value
def get_RT(row: pd.Series) -> float:
    """
    Returns the RT of a CSM.
    """

    return float(row["RT [min]"])

# get the CCS value
def get_CCS() -> float:
    """
    Dummy function.
    """
    return 0.0

# get the IonMobility value
def get_IonMobility(csm: pd.Series) -> float:
    """
    Dummy function.
    """
    return float(csm["Compensation Voltage"])

# get the values for all fragments of a CSM
def get_fragment_values(csm: pd.Series,
                        spectra: Dict,
                        crosslinker: str = CROSSLINKER,
                        possible_modifications: Dict[str, List[float]] = MODIFICATIONS,
                        ion_types: Tuple[str] = ION_TYPES,
                        max_charge: int = MAX_CHARGE,
                        match_tolerance: float = MATCH_TOLERANCE) -> Dict[str, List]:
    """
    Returns the annotated fragments of both cross-linked peptides.
    """

    fragments_A = get_fragments(csm, True, spectra, crosslinker, possible_modifications, ion_types, max_charge, match_tolerance)
    fragments_B = get_fragments(csm, False, spectra, crosslinker, possible_modifications, ion_types, max_charge, match_tolerance)

    return {"Fragments_A": fragments_A, "Fragments_B": fragments_B}

##### MAIN FUNCTION #####

# generates the spectral library
# def main(spectra_file: List[str] | List[BinaryIO] = SPECTRA_FILE,
#          csms_file: str | BinaryIO = CSMS_FILE,
# for backward compatibility >>
def main(spectra_file: Union[List[str], List[BinaryIO]] = SPECTRA_FILE,
         csms_file: Union[str, BinaryIO] = CSMS_FILE,
         run_name: str = RUN_NAME,
         crosslinker: str = CROSSLINKER,
         modifications: Dict[str, List[float]] = MODIFICATIONS,
         ion_types: Tuple[str] = ION_TYPES,
         max_charge: int = MAX_CHARGE,
         match_tolerance: float = MATCH_TOLERANCE,
         iRT_m: float = iRT_PARAMS["iRT_m"],
         iRT_t: float = iRT_PARAMS["iRT_t"],
         is_streamlit: bool = False,
         save_output: bool = True) -> Dict[str, pd.DataFrame]:

    if is_streamlit:
        print("INFO: Creating spectral library with input files:\nSpectra: " +
              "\n".join([spectrum_file.name for spectrum_file in spectra_file]) +
              "\nCSMs: " + str(csms_file.name))
    else:
        print("INFO: Creating spectral library with input files:\nSpectra: " + "\n".join(spectra_file) + "\nCSMs: " + csms_file)
    print("INFO: Using the following modifications:")
    print(modifications)
    print("INFO: Using the following ion types:")
    print(ion_types)
    print("INFO: Using the following charge states:")
    print([i for i in range(1, max_charge + 1)])
    print("INFO: Using a match tolerance of: " + str(match_tolerance) + " Da")
    print("INFO: Starting annotation process...")

    print("INFO: Reading spectra...")
    spectra = read_multiple_spectra(spectra_file) if not is_streamlit else read_multiple_spectra_streamlit(spectra_file)
    print("INFO: Done reading spectra!")

    print("INFO: Reading CSMs...")
    if "xlsx" in csms_file.split(".")[-1]:
        csms = pd.read_excel(csms_file)
    else:
        csms = parse_xi(csms_file, spectra)
        csms.to_csv(csms_file + ".converted.csv", index = False)
        if csms.shape[0] < 1000000:
            csms.to_excel(csms_file + ".converted.csv", index = False)
    print("INFO: Done reading CSMs! Starting spectral library creation...")

    # columns
    linkId_s = list()
    ProteinID_s = list()
    Organism_s = list()
    StrippedPeptide_s = list()
    FragmentGroupId_s = list()
    PrecursorCharge_s = list()
    PrecursorMz_s = list()
    ModifiedPeptide_s = list()
    IsotopeLabel_s = list()
    file_s = list()
    scanID_s = list()
    run_s = list()
    searchID_s = list()
    crosslinkedResidues_s = list()
    LabeledSequence_s = list()
    iRT_s = list()
    RT_s = list()
    CCS_s = list()
    IonMobility_s = list()
    FragmentCharge_s = list()
    FragmentType_s = list()
    FragmentNumber_s = list()
    FragmentPepId_s = list()
    FragmentMz_s = list()
    RelativeIntensity_s = list()
    FragmentLossType_s = list()
    CLContainingFragment_s = list()
    LossyFragment_s = list()
    Is_Decoy_s = list()
    Decoy_Type_s = list()

    # decoy dd columns
    linkId_s_decoy = list()
    ProteinID_s_decoy = list()
    Organism_s_decoy = list()
    StrippedPeptide_s_decoy = list()
    FragmentGroupId_s_decoy = list()
    PrecursorCharge_s_decoy = list()
    PrecursorMz_s_decoy = list()
    ModifiedPeptide_s_decoy = list()
    IsotopeLabel_s_decoy = list()
    file_s_decoy = list()
    scanID_s_decoy = list()
    run_s_decoy = list()
    searchID_s_decoy = list()
    crosslinkedResidues_s_decoy = list()
    LabeledSequence_s_decoy = list()
    iRT_s_decoy = list()
    RT_s_decoy = list()
    CCS_s_decoy = list()
    IonMobility_s_decoy = list()
    FragmentCharge_s_decoy = list()
    FragmentType_s_decoy = list()
    FragmentNumber_s_decoy = list()
    FragmentPepId_s_decoy = list()
    FragmentMz_s_decoy = list()
    RelativeIntensity_s_decoy = list()
    FragmentLossType_s_decoy = list()
    CLContainingFragment_s_decoy = list()
    LossyFragment_s_decoy = list()
    Is_Decoy_s_decoy = list()
    Decoy_Type_s_decoy = list()

    # decoy dt columns
    linkId_s_decoy_dt = list()
    ProteinID_s_decoy_dt = list()
    Organism_s_decoy_dt = list()
    StrippedPeptide_s_decoy_dt = list()
    FragmentGroupId_s_decoy_dt = list()
    PrecursorCharge_s_decoy_dt = list()
    PrecursorMz_s_decoy_dt = list()
    ModifiedPeptide_s_decoy_dt = list()
    IsotopeLabel_s_decoy_dt = list()
    file_s_decoy_dt = list()
    scanID_s_decoy_dt = list()
    run_s_decoy_dt = list()
    searchID_s_decoy_dt = list()
    crosslinkedResidues_s_decoy_dt = list()
    LabeledSequence_s_decoy_dt = list()
    iRT_s_decoy_dt = list()
    RT_s_decoy_dt = list()
    CCS_s_decoy_dt = list()
    IonMobility_s_decoy_dt = list()
    FragmentCharge_s_decoy_dt = list()
    FragmentType_s_decoy_dt = list()
    FragmentNumber_s_decoy_dt = list()
    FragmentPepId_s_decoy_dt = list()
    FragmentMz_s_decoy_dt = list()
    RelativeIntensity_s_decoy_dt = list()
    FragmentLossType_s_decoy_dt = list()
    CLContainingFragment_s_decoy_dt = list()
    LossyFragment_s_decoy_dt = list()
    Is_Decoy_s_decoy_dt = list()
    Decoy_Type_s_decoy_dt = list()

    # decoy td columns
    linkId_s_decoy_td = list()
    ProteinID_s_decoy_td = list()
    Organism_s_decoy_td = list()
    StrippedPeptide_s_decoy_td = list()
    FragmentGroupId_s_decoy_td = list()
    PrecursorCharge_s_decoy_td = list()
    PrecursorMz_s_decoy_td = list()
    ModifiedPeptide_s_decoy_td = list()
    IsotopeLabel_s_decoy_td = list()
    file_s_decoy_td = list()
    scanID_s_decoy_td = list()
    run_s_decoy_td = list()
    searchID_s_decoy_td = list()
    crosslinkedResidues_s_decoy_td = list()
    LabeledSequence_s_decoy_td = list()
    iRT_s_decoy_td = list()
    RT_s_decoy_td = list()
    CCS_s_decoy_td = list()
    IonMobility_s_decoy_td = list()
    FragmentCharge_s_decoy_td = list()
    FragmentType_s_decoy_td = list()
    FragmentNumber_s_decoy_td = list()
    FragmentPepId_s_decoy_td = list()
    FragmentMz_s_decoy_td = list()
    RelativeIntensity_s_decoy_td = list()
    FragmentLossType_s_decoy_td = list()
    CLContainingFragment_s_decoy_td = list()
    LossyFragment_s_decoy_td = list()
    Is_Decoy_s_decoy_td = list()
    Decoy_Type_s_decoy_td = list()

    # process CSMs
    for i, row in csms.iterrows():
        # target
        link_Id = get_linkId(row)
        ProteinID = get_ProteinID(row)
        Organism = get_Organism()
        StrippedPeptide = get_StrippedPeptide(row)
        FragmentGroupId = get_FragmentGroupId(row)
        PrecursorCharge = get_PrecursorCharge(row)
        PrecursorMz = get_PrecursorMz(row)
        ModifiedPeptide = get_ModifiedPeptide(row, crosslinker)
        IsotopeLabel = get_IsotopeLabel()
        cfile = get_filename(row)
        scanID = get_scanID(row)
        run = get_run(run_name)
        searchID = get_searchID(row)
        crosslinkedResidues = get_crosslinkedResidues(row)
        LabeledSequence = get_LabeledSequence(row)
        iRT = get_iRT(row, iRT_t, iRT_m)
        RT = get_RT(row)
        CCS = get_CCS()
        IonMobility = get_IonMobility(row)
        fragments = get_fragment_values(row, spectra, crosslinker, modifications, ion_types, max_charge, match_tolerance)

        for k in fragments.keys():
            pep = fragments[k]
            for frag in pep:
                linkId_s.append(link_Id)
                ProteinID_s.append(ProteinID)
                Organism_s.append(Organism)
                StrippedPeptide_s.append(StrippedPeptide)
                FragmentGroupId_s.append(FragmentGroupId)
                PrecursorCharge_s.append(PrecursorCharge)
                PrecursorMz_s.append(PrecursorMz)
                ModifiedPeptide_s.append(ModifiedPeptide)
                IsotopeLabel_s.append(IsotopeLabel)
                file_s.append(cfile)
                scanID_s.append(scanID)
                run_s.append(run)
                searchID_s.append(searchID)
                crosslinkedResidues_s.append(crosslinkedResidues)
                LabeledSequence_s.append(LabeledSequence)
                iRT_s.append(iRT)
                RT_s.append(RT)
                CCS_s.append(CCS)
                IonMobility_s.append(IonMobility)
                FragmentCharge_s.append(frag["FragmentCharge"])
                FragmentType_s.append(frag["FragmentType"])
                FragmentNumber_s.append(frag["FragmentNumber"])
                FragmentPepId_s.append(frag["FragmentPepId"])
                FragmentMz_s.append(frag["FragmentMz"])
                RelativeIntensity_s.append(frag["RelativeIntensity"])
                FragmentLossType_s.append(frag["FragmentLossType"])
                CLContainingFragment_s.append(frag["CLContainingFragment"])
                LossyFragment_s.append(frag["LossyFragment"])
                Is_Decoy_s.append(False)
                Decoy_Type_s.append("TT")

        # decoy dd
        decoy_csm = generate_decoy_csm_dd(row, crosslinker)
        decoy_link_Id = get_linkId(decoy_csm)
        decoy_ProteinID = get_ProteinID(decoy_csm)
        decoy_Organism = get_Organism()
        decoy_StrippedPeptide = get_StrippedPeptide(decoy_csm)
        decoy_FragmentGroupId = get_FragmentGroupId(decoy_csm)
        decoy_PrecursorCharge = get_PrecursorCharge(decoy_csm)
        decoy_PrecursorMz = get_PrecursorMz(decoy_csm)
        decoy_ModifiedPeptide = get_ModifiedPeptide(decoy_csm, crosslinker)
        decoy_IsotopeLabel = get_IsotopeLabel()
        decoy_cfile = get_filename(decoy_csm)
        decoy_scanID = get_scanID(decoy_csm)
        decoy_run = get_run(run_name)
        decoy_searchID = get_searchID(decoy_csm)
        decoy_crosslinkedResidues = get_crosslinkedResidues(decoy_csm)
        decoy_LabeledSequence = get_LabeledSequence(decoy_csm)
        decoy_iRT = get_iRT(decoy_csm, iRT_t, iRT_m)
        decoy_RT = get_RT(decoy_csm)
        decoy_CCS = get_CCS()
        decoy_IonMobility = get_IonMobility(decoy_csm)
        decoy_fragments = {"Fragments_A": get_decoy_fragments(decoy_csm, fragments["Fragments_A"], modifications, crosslinker),
                           "Fragments_B": get_decoy_fragments(decoy_csm, fragments["Fragments_B"], modifications, crosslinker)}

        for k in decoy_fragments.keys():
            decoy_pep = decoy_fragments[k]
            decoy_frag_mzs = list()
            for decoy_frag in decoy_pep:
                if decoy_frag["FragmentMz"] in decoy_frag_mzs:
                    continue
                linkId_s_decoy.append(decoy_link_Id)
                ProteinID_s_decoy.append(decoy_ProteinID)
                Organism_s_decoy.append(decoy_Organism)
                StrippedPeptide_s_decoy.append(decoy_StrippedPeptide)
                FragmentGroupId_s_decoy.append(decoy_FragmentGroupId)
                PrecursorCharge_s_decoy.append(decoy_PrecursorCharge)
                PrecursorMz_s_decoy.append(decoy_PrecursorMz)
                ModifiedPeptide_s_decoy.append(decoy_ModifiedPeptide)
                IsotopeLabel_s_decoy.append(decoy_IsotopeLabel)
                file_s_decoy.append(decoy_cfile)
                scanID_s_decoy.append(decoy_scanID)
                run_s_decoy.append(decoy_run)
                searchID_s_decoy.append(decoy_searchID)
                crosslinkedResidues_s_decoy.append(decoy_crosslinkedResidues)
                LabeledSequence_s_decoy.append(decoy_LabeledSequence)
                iRT_s_decoy.append(decoy_iRT)
                RT_s_decoy.append(decoy_RT)
                CCS_s_decoy.append(decoy_CCS)
                IonMobility_s_decoy.append(decoy_IonMobility)
                FragmentCharge_s_decoy.append(decoy_frag["FragmentCharge"])
                FragmentType_s_decoy.append(decoy_frag["FragmentType"])
                FragmentNumber_s_decoy.append(decoy_frag["FragmentNumber"])
                FragmentPepId_s_decoy.append(decoy_frag["FragmentPepId"])
                FragmentMz_s_decoy.append(decoy_frag["FragmentMz"])
                RelativeIntensity_s_decoy.append(decoy_frag["RelativeIntensity"])
                FragmentLossType_s_decoy.append(decoy_frag["FragmentLossType"])
                CLContainingFragment_s_decoy.append(decoy_frag["CLContainingFragment"])
                LossyFragment_s_decoy.append(decoy_frag["LossyFragment"])
                Is_Decoy_s_decoy.append(True)
                Decoy_Type_s_decoy.append("DD")
                decoy_frag_mzs.append(decoy_frag["FragmentMz"])

        # decoy dt
        decoy_csm_dt = generate_decoy_csm_dt(row, crosslinker)
        decoy_link_Id_dt = get_linkId(decoy_csm_dt)
        decoy_ProteinID_dt = get_ProteinID(decoy_csm_dt)
        decoy_Organism_dt = get_Organism()
        decoy_StrippedPeptide_dt = get_StrippedPeptide(decoy_csm_dt)
        decoy_FragmentGroupId_dt = get_FragmentGroupId(decoy_csm_dt)
        decoy_PrecursorCharge_dt = get_PrecursorCharge(decoy_csm_dt)
        decoy_PrecursorMz_dt = get_PrecursorMz(decoy_csm_dt)
        decoy_ModifiedPeptide_dt = get_ModifiedPeptide(decoy_csm_dt, crosslinker)
        decoy_IsotopeLabel_dt = get_IsotopeLabel()
        decoy_cfile_dt = get_filename(decoy_csm_dt)
        decoy_scanID_dt = get_scanID(decoy_csm_dt)
        decoy_run_dt = get_run(run_name)
        decoy_searchID_dt = get_searchID(decoy_csm_dt)
        decoy_crosslinkedResidues_dt = get_crosslinkedResidues(decoy_csm_dt)
        decoy_LabeledSequence_dt = get_LabeledSequence(decoy_csm_dt)
        decoy_iRT_dt = get_iRT(decoy_csm_dt, iRT_t, iRT_m)
        decoy_RT_dt = get_RT(decoy_csm_dt)
        decoy_CCS_dt = get_CCS()
        decoy_IonMobility_dt = get_IonMobility(decoy_csm_dt)
        decoy_fragments_dt = {"Fragments_A": get_decoy_fragments(decoy_csm, fragments["Fragments_A"], modifications, crosslinker),
                              "Fragments_B": fragments["Fragments_B"]}

        for k in decoy_fragments_dt.keys():
            decoy_pep_dt = decoy_fragments_dt[k]
            decoy_frag_mzs_dt = list()
            for decoy_frag_dt in decoy_pep_dt:
                if decoy_frag_dt["FragmentMz"] in decoy_frag_mzs_dt:
                    continue
                linkId_s_decoy_dt.append(decoy_link_Id_dt)
                ProteinID_s_decoy_dt.append(decoy_ProteinID_dt)
                Organism_s_decoy_dt.append(decoy_Organism_dt)
                StrippedPeptide_s_decoy_dt.append(decoy_StrippedPeptide_dt)
                FragmentGroupId_s_decoy_dt.append(decoy_FragmentGroupId_dt)
                PrecursorCharge_s_decoy_dt.append(decoy_PrecursorCharge_dt)
                PrecursorMz_s_decoy_dt.append(decoy_PrecursorMz_dt)
                ModifiedPeptide_s_decoy_dt.append(decoy_ModifiedPeptide_dt)
                IsotopeLabel_s_decoy_dt.append(decoy_IsotopeLabel_dt)
                file_s_decoy_dt.append(decoy_cfile_dt)
                scanID_s_decoy_dt.append(decoy_scanID_dt)
                run_s_decoy_dt.append(decoy_run_dt)
                searchID_s_decoy_dt.append(decoy_searchID_dt)
                crosslinkedResidues_s_decoy_dt.append(decoy_crosslinkedResidues_dt)
                LabeledSequence_s_decoy_dt.append(decoy_LabeledSequence_dt)
                iRT_s_decoy_dt.append(decoy_iRT_dt)
                RT_s_decoy_dt.append(decoy_RT_dt)
                CCS_s_decoy_dt.append(decoy_CCS_dt)
                IonMobility_s_decoy_dt.append(decoy_IonMobility_dt)
                FragmentCharge_s_decoy_dt.append(decoy_frag_dt["FragmentCharge"])
                FragmentType_s_decoy_dt.append(decoy_frag_dt["FragmentType"])
                FragmentNumber_s_decoy_dt.append(decoy_frag_dt["FragmentNumber"])
                FragmentPepId_s_decoy_dt.append(decoy_frag_dt["FragmentPepId"])
                FragmentMz_s_decoy_dt.append(decoy_frag_dt["FragmentMz"])
                RelativeIntensity_s_decoy_dt.append(decoy_frag_dt["RelativeIntensity"])
                FragmentLossType_s_decoy_dt.append(decoy_frag_dt["FragmentLossType"])
                CLContainingFragment_s_decoy_dt.append(decoy_frag_dt["CLContainingFragment"])
                LossyFragment_s_decoy_dt.append(decoy_frag_dt["LossyFragment"])
                Is_Decoy_s_decoy_dt.append(True)
                Decoy_Type_s_decoy_dt.append("DT")
                decoy_frag_mzs_dt.append(decoy_frag_dt["FragmentMz"])

        # decoy td
        decoy_csm_td = generate_decoy_csm_td(row, crosslinker)
        decoy_link_Id_td = get_linkId(decoy_csm_td)
        decoy_ProteinID_td = get_ProteinID(decoy_csm_td)
        decoy_Organism_td = get_Organism()
        decoy_StrippedPeptide_td = get_StrippedPeptide(decoy_csm_td)
        decoy_FragmentGroupId_td = get_FragmentGroupId(decoy_csm_td)
        decoy_PrecursorCharge_td = get_PrecursorCharge(decoy_csm_td)
        decoy_PrecursorMz_td = get_PrecursorMz(decoy_csm_td)
        decoy_ModifiedPeptide_td = get_ModifiedPeptide(decoy_csm_td, crosslinker)
        decoy_IsotopeLabel_td = get_IsotopeLabel()
        decoy_cfile_td = get_filename(decoy_csm_td)
        decoy_scanID_td = get_scanID(decoy_csm_td)
        decoy_run_td = get_run(run_name)
        decoy_searchID_td = get_searchID(decoy_csm_td)
        decoy_crosslinkedResidues_td = get_crosslinkedResidues(decoy_csm_td)
        decoy_LabeledSequence_td = get_LabeledSequence(decoy_csm_td)
        decoy_iRT_td = get_iRT(decoy_csm_td, iRT_t, iRT_m)
        decoy_RT_td = get_RT(decoy_csm_td)
        decoy_CCS_td = get_CCS()
        decoy_IonMobility_td = get_IonMobility(decoy_csm_td)
        decoy_fragments_td = {"Fragments_A": fragments["Fragments_A"],
                              "Fragments_B": get_decoy_fragments(decoy_csm, fragments["Fragments_B"], modifications, crosslinker)}

        for k in decoy_fragments_td.keys():
            decoy_pep_td = decoy_fragments_td[k]
            decoy_frag_mzs_td = list()
            for decoy_frag_td in decoy_pep_td:
                if decoy_frag_td["FragmentMz"] in decoy_frag_mzs_td:
                    continue
                linkId_s_decoy_td.append(decoy_link_Id_td)
                ProteinID_s_decoy_td.append(decoy_ProteinID_td)
                Organism_s_decoy_td.append(decoy_Organism_td)
                StrippedPeptide_s_decoy_td.append(decoy_StrippedPeptide_td)
                FragmentGroupId_s_decoy_td.append(decoy_FragmentGroupId_td)
                PrecursorCharge_s_decoy_td.append(decoy_PrecursorCharge_td)
                PrecursorMz_s_decoy_td.append(decoy_PrecursorMz_td)
                ModifiedPeptide_s_decoy_td.append(decoy_ModifiedPeptide_td)
                IsotopeLabel_s_decoy_td.append(decoy_IsotopeLabel_td)
                file_s_decoy_td.append(decoy_cfile_td)
                scanID_s_decoy_td.append(decoy_scanID_td)
                run_s_decoy_td.append(decoy_run_td)
                searchID_s_decoy_td.append(decoy_searchID_td)
                crosslinkedResidues_s_decoy_td.append(decoy_crosslinkedResidues_td)
                LabeledSequence_s_decoy_td.append(decoy_LabeledSequence_td)
                iRT_s_decoy_td.append(decoy_iRT_td)
                RT_s_decoy_td.append(decoy_RT_td)
                CCS_s_decoy_td.append(decoy_CCS_td)
                IonMobility_s_decoy_td.append(decoy_IonMobility_td)
                FragmentCharge_s_decoy_td.append(decoy_frag_td["FragmentCharge"])
                FragmentType_s_decoy_td.append(decoy_frag_td["FragmentType"])
                FragmentNumber_s_decoy_td.append(decoy_frag_td["FragmentNumber"])
                FragmentPepId_s_decoy_td.append(decoy_frag_td["FragmentPepId"])
                FragmentMz_s_decoy_td.append(decoy_frag_td["FragmentMz"])
                RelativeIntensity_s_decoy_td.append(decoy_frag_td["RelativeIntensity"])
                FragmentLossType_s_decoy_td.append(decoy_frag_td["FragmentLossType"])
                CLContainingFragment_s_decoy_td.append(decoy_frag_td["CLContainingFragment"])
                LossyFragment_s_decoy_td.append(decoy_frag_td["LossyFragment"])
                Is_Decoy_s_decoy_td.append(True)
                Decoy_Type_s_decoy_td.append("TD")
                decoy_frag_mzs_td.append(decoy_frag_td["FragmentMz"])

        if (i + 1) % 100 == 0:
            print("INFO: Processed " + str(i + 1) + " CSMs in total...")

    # generate dataframe
    tt_dict = {"linkId": linkId_s,
               "ProteinID": ProteinID_s,
               "Organism": Organism_s,
               "StrippedPeptide": StrippedPeptide_s,
               "FragmentGroupId": FragmentGroupId_s,
               "PrecursorCharge": PrecursorCharge_s,
               "PrecursorMz": PrecursorMz_s,
               "ModifiedPeptide": ModifiedPeptide_s,
               "IsotopeLabel": IsotopeLabel_s,
               "File": file_s,
               "scanID": scanID_s,
               "run": run_s,
               "searchID": searchID_s,
               "crosslinkedResidues": crosslinkedResidues_s,
               "LabeledSequence": LabeledSequence_s,
               "iRT": iRT_s,
               "RT": RT_s,
               "CCS": CCS_s,
               "IonMobility": IonMobility_s,
               "FragmentCharge": FragmentCharge_s,
               "FragmentType": FragmentType_s,
               "FragmentNumber": FragmentNumber_s,
               "FragmentPepId": FragmentPepId_s,
               "FragmentMz": FragmentMz_s,
               "RelativeIntensity": RelativeIntensity_s,
               "FragmentLossType": FragmentLossType_s,
               "CLContainingFragment": CLContainingFragment_s,
               "LossyFragment": LossyFragment_s,
               "IsDecoy": Is_Decoy_s,
               "DecoyType": Decoy_Type_s}

    spectral_library = pd.DataFrame(tt_dict)

    dd_dict = {"linkId": linkId_s_decoy,
               "ProteinID": ProteinID_s_decoy,
               "Organism": Organism_s_decoy,
               "StrippedPeptide": StrippedPeptide_s_decoy,
               "FragmentGroupId": FragmentGroupId_s_decoy,
               "PrecursorCharge": PrecursorCharge_s_decoy,
               "PrecursorMz": PrecursorMz_s_decoy,
               "ModifiedPeptide": ModifiedPeptide_s_decoy,
               "IsotopeLabel": IsotopeLabel_s_decoy,
               "File": file_s_decoy,
               "scanID": scanID_s_decoy,
               "run": run_s_decoy,
               "searchID": searchID_s_decoy,
               "crosslinkedResidues": crosslinkedResidues_s_decoy,
               "LabeledSequence": LabeledSequence_s_decoy,
               "iRT": iRT_s_decoy,
               "RT": RT_s_decoy,
               "CCS": CCS_s_decoy,
               "IonMobility": IonMobility_s_decoy,
               "FragmentCharge": FragmentCharge_s_decoy,
               "FragmentType": FragmentType_s_decoy,
               "FragmentNumber": FragmentNumber_s_decoy,
               "FragmentPepId": FragmentPepId_s_decoy,
               "FragmentMz": FragmentMz_s_decoy,
               "RelativeIntensity": RelativeIntensity_s_decoy,
               "FragmentLossType": FragmentLossType_s_decoy,
               "CLContainingFragment": CLContainingFragment_s_decoy,
               "LossyFragment": LossyFragment_s_decoy,
               "IsDecoy": Is_Decoy_s_decoy,
               "DecoyType": Decoy_Type_s_decoy}

    spectral_library_decoy_dd = pd.DataFrame(dd_dict)

    dt_dict = {"linkId": linkId_s_decoy_dt,
               "ProteinID": ProteinID_s_decoy_dt,
               "Organism": Organism_s_decoy_dt,
               "StrippedPeptide": StrippedPeptide_s_decoy_dt,
               "FragmentGroupId": FragmentGroupId_s_decoy_dt,
               "PrecursorCharge": PrecursorCharge_s_decoy_dt,
               "PrecursorMz": PrecursorMz_s_decoy_dt,
               "ModifiedPeptide": ModifiedPeptide_s_decoy_dt,
               "IsotopeLabel": IsotopeLabel_s_decoy_dt,
               "File": file_s_decoy_dt,
               "scanID": scanID_s_decoy_dt,
               "run": run_s_decoy_dt,
               "searchID": searchID_s_decoy_dt,
               "crosslinkedResidues": crosslinkedResidues_s_decoy_dt,
               "LabeledSequence": LabeledSequence_s_decoy_dt,
               "iRT": iRT_s_decoy_dt,
               "RT": RT_s_decoy_dt,
               "CCS": CCS_s_decoy_dt,
               "IonMobility": IonMobility_s_decoy_dt,
               "FragmentCharge": FragmentCharge_s_decoy_dt,
               "FragmentType": FragmentType_s_decoy_dt,
               "FragmentNumber": FragmentNumber_s_decoy_dt,
               "FragmentPepId": FragmentPepId_s_decoy_dt,
               "FragmentMz": FragmentMz_s_decoy_dt,
               "RelativeIntensity": RelativeIntensity_s_decoy_dt,
               "FragmentLossType": FragmentLossType_s_decoy_dt,
               "CLContainingFragment": CLContainingFragment_s_decoy_dt,
               "LossyFragment": LossyFragment_s_decoy_dt,
               "IsDecoy": Is_Decoy_s_decoy_dt,
               "DecoyType": Decoy_Type_s_decoy_dt}

    spectral_library_decoy_dt = pd.DataFrame(dt_dict)

    td_dict = {"linkId": linkId_s_decoy_td,
               "ProteinID": ProteinID_s_decoy_td,
               "Organism": Organism_s_decoy_td,
               "StrippedPeptide": StrippedPeptide_s_decoy_td,
               "FragmentGroupId": FragmentGroupId_s_decoy_td,
               "PrecursorCharge": PrecursorCharge_s_decoy_td,
               "PrecursorMz": PrecursorMz_s_decoy_td,
               "ModifiedPeptide": ModifiedPeptide_s_decoy_td,
               "IsotopeLabel": IsotopeLabel_s_decoy_td,
               "File": file_s_decoy_td,
               "scanID": scanID_s_decoy_td,
               "run": run_s_decoy_td,
               "searchID": searchID_s_decoy_td,
               "crosslinkedResidues": crosslinkedResidues_s_decoy_td,
               "LabeledSequence": LabeledSequence_s_decoy_td,
               "iRT": iRT_s_decoy_td,
               "RT": RT_s_decoy_td,
               "CCS": CCS_s_decoy_td,
               "IonMobility": IonMobility_s_decoy_td,
               "FragmentCharge": FragmentCharge_s_decoy_td,
               "FragmentType": FragmentType_s_decoy_td,
               "FragmentNumber": FragmentNumber_s_decoy_td,
               "FragmentPepId": FragmentPepId_s_decoy_td,
               "FragmentMz": FragmentMz_s_decoy_td,
               "RelativeIntensity": RelativeIntensity_s_decoy_td,
               "FragmentLossType": FragmentLossType_s_decoy_td,
               "CLContainingFragment": CLContainingFragment_s_decoy_td,
               "LossyFragment": LossyFragment_s_decoy_td,
               "IsDecoy": Is_Decoy_s_decoy_td,
               "DecoyType": Decoy_Type_s_decoy_td}

    spectral_library_decoy_td = pd.DataFrame(td_dict)

    if save_output:
        # save spectral library
        spectral_library.to_csv(".".join(csms_file.split(".")[:-1]) + "_spectralLibrary.csv", index = True)
        spectral_library_decoy_dd.to_csv(".".join(csms_file.split(".")[:-1]) + "_spectralLibraryDECOY_DD.csv", index = True)
        spectral_library_decoy_dt.to_csv(".".join(csms_file.split(".")[:-1]) + "_spectralLibraryDECOY_DT.csv", index = True)
        spectral_library_decoy_td.to_csv(".".join(csms_file.split(".")[:-1]) + "_spectralLibraryDECOY_TD.csv", index = True)

        print("SUCCESS: Spectral library created with filename:")
        print(".".join(csms_file.split(".")[:-1]) + "_spectralLibrary.csv")
        print("SUCCESS: Decoy Spectral libraries created with filenames:")
        print(".".join(csms_file.split(".")[:-1]) + "_spectralLibraryDECOY_DD.csv")
        print(".".join(csms_file.split(".")[:-1]) + "_spectralLibraryDECOY_DT.csv")
        print(".".join(csms_file.split(".")[:-1]) + "_spectralLibraryDECOY_TD.csv")

        print("Creating merged library...")
        merged_spec_lib = pd.concat([spectral_library, spectral_library_decoy_dd, spectral_library_decoy_dt, spectral_library_decoy_td], ignore_index = True)
        should_shape = spectral_library.shape[0] + spectral_library_decoy_dd.shape[0] + spectral_library_decoy_dt.shape[0] + spectral_library_decoy_td.shape[0]
        if merged_spec_lib.shape[0] != should_shape:
            warnings.warn(f"Merged spectral library has {merged_spec_lib.shape[0]} rows, should be {should_shape} rows! Potential loss of data, please consider merging manually!")
        merged_spec_lib.to_csv(".".join(csms_file.split(".")[:-1]) + "_spectralLibraryFULL.csv", index = True)
        print("SUCCESS: Merged spectral library created with filename:")
        print(".".join(csms_file.split(".")[:-1]) + "_spectralLibraryFULL.csv")

    return {"TargetLib": spectral_library, "DecoyLib": spectral_library_decoy_dd, "DecoyLib_DT": spectral_library_decoy_dt, "DecoyLib_TD": spectral_library_decoy_td, "FullLib": merged_spec_lib}

##### SCRIPT #####

if __name__ == "__main__":

    sl = main()
