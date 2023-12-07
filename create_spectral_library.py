#!/usr/bin/env python3

# MS ANNIKA SPECTRAL LIBRARY EXPORTER
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# version tracking
__version = "1.1.1"
__date = "2023-12-06"

# REQUIREMENTS
# pip install pandas
# pip install openpyxl
# pip install pyteomics

##### PARAMETERS #####

from config import *

######################

# import packages
import pandas as pd
from pyteomics import mgf, mass

from typing import Dict
from typing import List
from typing import Tuple
from typing import Set
from typing import Union
from typing import BinaryIO
import warnings

# reading spectra
# def read_spectra(filename: str | BinaryIO) -> Dict[int, Dict]:
# for backward compatibility >>
def read_spectra(filename: Union[str, BinaryIO]) -> Dict[int, Dict]:
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
    def parse_mod_str(mod_str, crosslinker):
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
        return string[:index] + character + string[:index]
    # end function

    mods_A = parse_mod_str(str(row["Modifications A"]), crosslinker)
    mods_B = parse_mod_str(str(row["Modifications B"]), crosslinker)

    # generate annotation for sequence A
    shift = 0
    mod_A_template_str = str(row["Sequence A"])
    for pos in mods_A.keys():
        current_mods = "(" + ", ".join(mods_A[pos]) + ")"
        mod_A_template_str = str_insert(mod_A_template_str, pos + shift, current_mods)
        shift += len(current_mods)

    # generate annotation for sequence B
    shift = 0
    mod_B_template_str = str(row["Sequence B"])
    for pos in mods_B.keys():
        current_mods = "(" + ", ".join(mods_B[pos]) + ")"
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
def get_IonMobility() -> float:
    """
    Dummy function.
    """
    return 0.0

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
         save_output: bool = True) -> pd.DataFrame:

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
    csms = pd.read_excel(csms_file)
    print("INFO: Done reading CSMs! Starting spectral library creation...")

    # columns
    linkId_s = list()
    ProteinID_s = list()
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

    # process CSMs
    for i, row in csms.iterrows():
        link_Id = get_linkId(row)
        ProteinID = get_ProteinID(row)
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
        IonMobility = get_IonMobility()
        fragments = get_fragment_values(row, spectra, crosslinker, modifications, ion_types, max_charge, match_tolerance)

        for k in fragments.keys():
            pep = fragments[k]
            for frag in pep:
                linkId_s.append(link_Id)
                ProteinID_s.append(ProteinID)
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

        if (i + 1) % 100 == 0:
            print("INFO: Processed " + str(i + 1) + " CSMs in total...")

    # generate dataframe
    df_dict = {"linkId": linkId_s,
               "ProteinID": ProteinID_s,
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
               "LossyFragment": LossyFragment_s}

    spectral_library = pd.DataFrame(df_dict)

    if save_output:
        # save spectral library
        spectral_library.to_csv(".".join(csms_file.split(".")[:-1]) + "_spectralLibrary.csv", index = True)

        print("SUCCESS: Spectral library created with filename:")
        print(".".join(csms_file.split(".")[:-1]) + "_spectralLibrary.csv")

    return spectral_library

##### SCRIPT #####

if __name__ == "__main__":

        sl = main()
