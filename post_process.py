#!/usr/bin/env python3

# SPECTRONAUT POST PROCESSING
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# version tracking
__version = "1.0.0"
__date = "2025-02-20"

# PARAMETERS

SPECTRONAUT_DELIM = "," #"\t" # delimiter in Spectronaut output file
SPECTRONAUT_MATCH_TOLERANCE = 0.05 # match tolerance in Da
SPECTRONAUT_FRAGMENT_MZ_COLUMN_NAME = "F.CalibratedMz" # which F Mz to use for matching
SPECTRONAUT_CSCORE_COLUMN_NAME = "PG.Cscore (Run-Wise)" # which Cscore to use for re-soring

# REQUIREMENTS
# pip install tqdm
# pip install pandas

# import packages
from tqdm import tqdm
import pandas as pd

from typing import Any
from typing import Dict
from typing import List

##################### POST PROCESS #####################

def get_mz_key(mz: float) -> float:
    #return str(int(mz * 1000))
    return mz
    
def get_fragment_key(mz: float) -> str:
    return f"{round(mz, 4):.4f}"

def get_key_spec_lib(row: pd.Series) -> str:
    # ModifiedPeptide
    # DAKQRIVDK_NGVKM[Oxidation]C[Carbamidomethyl]PR
    # PrecursorCharge
    # 4
    # linkId
    # A0A4V0IIP8_Q21298-126_45
    return f"{str(row['ModifiedPeptide']).strip()}.{int(row['PrecursorCharge'])}:{str(row['linkId']).strip()}"

def get_key_spectronaut(row: pd.Series) -> str:
    # EG.PrecursorId
    # AAHHADGLAKGLHETC[Carbamidomethyl]K_M[Oxidation]FIPKSHTK.5
    # PG.ProteinNames
    # P10771_P10771-113_46
    # ---> if PG.ProteinNames
    # P10771_P10771
    # use FG.Comment
    # 113_46
    if "-" in str(row["PG.ProteinNames"]):
        return f"{str(row['EG.PrecursorId']).strip()}:{str(row['PG.ProteinNames']).strip()}"
    return f"{str(row['EG.PrecursorId']).strip()}:{str(row['PG.ProteinNames']).strip()}-{str(row['FG.Comment']).strip()}"

def read_spectral_library(filename: str) -> Dict[str, Dict[str, Any]]:

    spec_lib = pd.read_csv(filename, low_memory = False)

    index = dict()
    for i, row in tqdm(spec_lib.iterrows(), total = spec_lib.shape[0], desc = "Creating index..."):
        key = get_key_spec_lib(row)
        if key in index:
            index[key]["rows"].append(row)

            if int(row["FragmentPepId"]) == 0:
                index[key]["total_ions_a"] += 1
            else:
                index[key]["total_ions_b"] += 1

            ion_mz = get_mz_key(float(row["FragmentMz"]))
            if ion_mz in index[key]["ions"]:
                index[key]["ions"][ion_mz].append(row)
            else:
                index[key]["ions"][ion_mz] = [row]
        else:
            index[key] = {"rows": [row],
                          "total_ions_a": 1 if int(row["FragmentPepId"]) == 0 else 0,
                          "total_ions_b": 1 if int(row["FragmentPepId"]) == 1 else 0,
                          "ions": {get_mz_key(float(row["FragmentMz"])): [row]}}

    return index

def generate_fragment_index(spectronaut: pd.DataFrame, index: dict) -> Dict[str, Dict[str, int]]:

    fragment_annotation = dict()
    for i, row in tqdm(spectronaut.iterrows(), total = spectronaut.shape[0], desc = "Generating fragment ion index..."):
        key = get_key_spectronaut(row)
        ion = float(row[SPECTRONAUT_FRAGMENT_MZ_COLUMN_NAME])
        ions = index[key]["ions"]
        for current_ion_mz in ions.keys():
            fragment_key_part = get_fragment_key(current_ion_mz)
            if round(current_ion_mz, 4) > round(ion - SPECTRONAUT_MATCH_TOLERANCE, 4) and round(current_ion_mz, 4) < round(ion + SPECTRONAUT_MATCH_TOLERANCE, 4):
                if key not in fragment_annotation:
                    fragment_annotation[key] = {"matched_number_ions_a": 0,
                                                "matched_number_ions_b": 0,
                                                "fragments": []}
                for current_ion in ions[current_ion_mz]:
                    if current_ion["FragmentPepId"] == 0:
                        fragment_key = fragment_key_part + "_0"
                        if fragment_key not in fragment_annotation[key]["fragments"]:
                            fragment_annotation[key]["matched_number_ions_a"] += 1
                            fragment_annotation[key]["fragments"].append(fragment_key)
                    else:
                        fragment_key = fragment_key_part + "_1"
                        if fragment_key not in fragment_annotation[key]["fragments"]:
                            fragment_annotation[key]["matched_number_ions_b"] += 1
                            fragment_annotation[key]["fragments"].append(fragment_key)

    return fragment_annotation

def annotate_spectronaut_result(filename: str) -> pd.DataFrame:

    spectronaut = pd.read_csv(filename, sep = SPECTRONAUT_DELIM, low_memory = False)
    filename_spec_lib = str(spectronaut["EG.Library"].at[0])
    index = read_spectral_library(filename_spec_lib)

    ## available columns in spec lib
    # COLUMN                                            EXAMPLE
    # linkId                                            P18948_P18948-516_509
    # ProteinID                                         P18948_P18948
    # Organism                                          Caenorhabditis elegans
    # StrippedPeptide                                   DAKQRIVDKNGVKMCPR
    # FragmentGroupId                                   DAKQRIVDK_NGVKMCPR-3_4:5
    # PrecursorCharge                                   5
    # PrecursorMz                                       429.823708
    # ModifiedPeptide                                   DAKQRIVDK_NGVKM[Oxidation]C[Carbamidomethyl]PR
    # IsotopeLabel                                      0
    # File                                              20240131_E0_NEO6_Mueller_MS_TechHub_IMP_THIDVJSSG004_Celegans_Nulcei_S3_DSG_SEC9_allng_FAIMS_001.raw
    # scanID                                            15497
    # run                                               20240131_E0_NEO6_THIDVJSSG004_Celegans_Nulcei_S1toS3_DSG_SEC_entrapment
    # searchID                                          MS Annika
    # crosslinkedResidues                               516_509
    # LabeledSequence                                   DAKQRIVDK_NGVKM[Oxidation]C[Carbamidomethyl]PR
    # iRT                                               -17.5246511627907
    # RT                                                32.1844
    # CCS                                               0
    # IonMobility                                       -50
    # FragmentCharge                                    1
    # FragmentType                                      b
    # FragmentNumber                                    2
    # FragmentPepId                                     0
    # FragmentMz                                        187.0715942
    # RelativeIntensity                                 0.115658057537211
    # FragmentLossType                                  None
    # CLContainingFragment                              False
    # LossyFragment                                     False
    # isDecoy                                           False
    # DecoyType                                         TT

    # for rescoring
    # a: matched number ions a
    # b: total number ions a
    # c: matched number ions b
    # d: total number ions b
    # e: relative match score a (a / b)
    # f: relative match score b (c / d)
    # g: partial c score a (c score * (a / a + c))
    # h: partial c score b (c score * (c / a + c))
    # i: composite relative match score min(e, f)
    # j: composite partial c score min(g, h)
    fragment_annotation = generate_fragment_index(spectronaut, index)

    # a
    def annotate_MatchedIonsA(row: pd.Series, fragment_annotation: dict) -> int:
        key = get_key_spectronaut(row)
        return fragment_annotation[key]["matched_number_ions_a"]

    tqdm.pandas(desc = "Annotating number of matched ions A...")
    spectronaut["MatchedIonsA"] = spectronaut.progress_apply(lambda row: annotate_MatchedIonsA(row, fragment_annotation), axis = 1)

    # b
    def annotate_TotalIonsA(row: pd.Series, index: dict) -> int:
        key = get_key_spectronaut(row)
        return index[key]["total_ions_a"]

    tqdm.pandas(desc = "Annotating number of total ions A...")
    spectronaut["TotalIonsA"] = spectronaut.progress_apply(lambda row: annotate_TotalIonsA(row, index), axis = 1)

    # c
    def annotate_MatchedIonsB(row: pd.Series, fragment_annotation: dict) -> int:
        key = get_key_spectronaut(row)
        return fragment_annotation[key]["matched_number_ions_b"]

    tqdm.pandas(desc = "Annotating number of matched ions B...")
    spectronaut["MatchedIonsB"] = spectronaut.progress_apply(lambda row: annotate_MatchedIonsB(row, fragment_annotation), axis = 1)

    # d
    def annotate_TotalIonsB(row: pd.Series, index: dict) -> int:
        key = get_key_spectronaut(row)
        return index[key]["total_ions_b"]

    tqdm.pandas(desc = "Annotating number of total ions B...")
    spectronaut["TotalIonsB"] = spectronaut.progress_apply(lambda row: annotate_TotalIonsB(row, index), axis = 1)

    # e
    tqdm.pandas(desc = "Annotating relative match score A...")
    spectronaut["RelativeMatchScoreA"] = spectronaut.progress_apply(lambda row: row["MatchedIonsA"] / row["TotalIonsA"], axis = 1)

    # f
    tqdm.pandas(desc = "Annotating relative match score B...")
    spectronaut["RelativeMatchScoreB"] = spectronaut.progress_apply(lambda row: row["MatchedIonsB"] / row["TotalIonsB"], axis = 1)

    # g
    def annotate_PartialCscoreA(row: pd.Series) -> float:
        cscore = row[SPECTRONAUT_CSCORE_COLUMN_NAME]
        partial = row["MatchedIonsA"] / (row["MatchedIonsA"] + row["MatchedIonsB"])
        return cscore * partial

    tqdm.pandas(desc = "Annotating partial Cscore A...")
    spectronaut["PartialCscoreA"] = spectronaut.progress_apply(lambda row: annotate_PartialCscoreA(row), axis = 1)

    # h
    def annotate_PartialCscoreB(row: pd.Series) -> float:
        cscore = row[SPECTRONAUT_CSCORE_COLUMN_NAME]
        partial = row["MatchedIonsB"] / (row["MatchedIonsA"] + row["MatchedIonsB"])
        return cscore * partial

    tqdm.pandas(desc = "Annotating partial Cscore B...")
    spectronaut["PartialCscoreB"] = spectronaut.progress_apply(lambda row: annotate_PartialCscoreB(row), axis = 1)

    # i
    tqdm.pandas(desc = "Annotating composite relative match score...")
    spectronaut["CompositeRelativeMatchScore"] = spectronaut.progress_apply(lambda row: min(row["RelativeMatchScoreA"], row["RelativeMatchScoreB"]), axis = 1)

    # j
    tqdm.pandas(desc = "Annotating composite partial Cscore...")
    spectronaut["CompositePartialCscore"] = spectronaut.progress_apply(lambda row: min(row["PartialCscoreA"], row["PartialCscoreB"]), axis = 1)

    # annotation of other properties
    def annotate_DecoyType(row: pd.Series, index: dict) -> str:
        key = get_key_spectronaut(row)
        return str(index[key][0]["DecoyType"]).strip()

    tqdm.pandas(desc = "Annotating DecoyType...")
    spectronaut["DecoyType"] = spectronaut.progress_apply(lambda row: annotate_DecoyType(row, index), axis = 1)

    return spectronaut

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        raise RuntimeError("No Spectronaut file was given!")

    filename = sys.argv[1]
    filename_o = filename + "_annotated.csv"

    df = annotate_spectronaut_result(filename)
    print("Writing annotated Spectronaut result to file...")
    df.to_csv(filename_o, index = False)
    print("Finished annotation!")
