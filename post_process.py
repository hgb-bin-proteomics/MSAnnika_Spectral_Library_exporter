#!/usr/bin/env python3

# SPECTRONAUT POST PROCESSING
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# version tracking
__version = "1.0.0"
__date = "2025-02-20"

# PARAMETERS
SPECTRONAUT_MATCH_TOLERANCE = 0.02 # match tolerance in Da
SPECTRONAUT_DELIM = "," #"\t" # delimiter in Spectronaut output file
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

def read_spectral_library(filename: str) -> Dict[str, List[pd.Series]]:

    spec_lib = pd.read_csv(filename, low_memory = False)

    index = dict()
    for i, row in tqdm(spec_lib.iterrows(), total = spec_lib.shape[0], desc = "Creating index..."):
        key = get_key_spec_lib(row)
        if key in index:
            index[key].append(row)
        else:
            index[key] = [row]

    return index

def annotate_spectronaut_result(filename: str) -> pd.DataFrame:

    spectronaut = pd.read_csv(filename, sep = "\t", low_memory = False)
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

    def annotate_DecoyType(row: pd.Series, index: dict) -> str:
        key = get_key_spectronaut(row)
        return str(index[key][0]["DecoyType"]).strip()

    tqdm.pandas(desc = "Annotating DecoyType...")
    spectronaut["DecoyType"] = spectronaut.progress_apply(lambda row: annotate_DecoyType(row, index), axis = 1)

    return spectronaut
