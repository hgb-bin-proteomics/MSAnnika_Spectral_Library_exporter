#!/usr/bin/env python3

# MS ANNIKA SPECTRAL LIBRARY EXPORTER - TESTS XI
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# check converted csm file
def test1_converter():

    import pytest
    import pandas as pd
    from create_spectral_library import main
    from config import CSMS_FILE as CSMS_FILE

    sl = main()
    csms = pd.read_csv(CSMS_FILE + ".converted.csv")

    assert csms.shape[0] == 2
    assert csms.shape[1] == 18

    assert str(csms.loc[0, "Sequence A"]) == "KIECFDSVEISGVEDR"
    assert str(csms.loc[0, "Modifications A"]) == "C4(Carbamidomethyl);K1(BS3)"
    assert str(csms.loc[0, "Sequence B"]) == "KIECFDSVEISGVEDR"
    assert str(csms.loc[0, "Modifications B"]) == "C4(Carbamidomethyl);K1(BS3)"
    assert str(csms.loc[0, "First Scan"]) == "19140"
    assert str(csms.loc[0, "Spectrum File"]) == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert str(csms.loc[0, "A in protein"]) == "574"
    assert str(csms.loc[0, "B in protein"]) == "574"
    assert str(csms.loc[0, "Crosslinker Position A"]) == "1"
    assert str(csms.loc[0, "Crosslinker Position B"]) == "1"
    assert str(csms.loc[0, "Accession A"]) == "Cas9"
    assert str(csms.loc[0, "Accession B"]) == "Cas9"
    assert str(csms.loc[0, "Charge"]) == "4"
    assert str(int(csms.loc[0, "m/z [Da]"])) == "976"
    assert str(csms.loc[0, "Crosslink Strategy"]) == "xi"
    assert str(int(csms.loc[0, "RT [min]"])) == "85"
    assert str(int(csms.loc[0, "Compensation Voltage"])) == "0"
    assert float(csms.loc[0, "Combined Score"]) == pytest.approx(27.268)

    assert str(csms.loc[1, "Sequence A"]) == "MIAKSEQEIGK"
    assert str(csms.loc[1, "Modifications A"]) == "K4(BS3)"
    assert str(csms.loc[1, "Sequence B"]) == "DFQFYKVR"
    assert str(csms.loc[1, "Modifications B"]) == "K6(BS3)"
    assert str(csms.loc[1, "First Scan"]) == "15867"
    assert str(csms.loc[1, "Spectrum File"]) == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert str(csms.loc[1, "A in protein"]) == "1024"
    assert str(csms.loc[1, "B in protein"]) == "972"
    assert str(csms.loc[1, "Crosslinker Position A"]) == "4"
    assert str(csms.loc[1, "Crosslinker Position B"]) == "6"
    assert str(csms.loc[1, "Accession A"]) == "Cas9"
    assert str(csms.loc[1, "Accession B"]) == "Cas9"
    assert str(csms.loc[1, "Charge"]) == "3"
    assert str(int(csms.loc[1, "m/z [Da]"])) == "825"
    assert str(csms.loc[1, "Crosslink Strategy"]) == "xi"
    assert str(int(csms.loc[1, "RT [min]"])) == "72"
    assert str(int(csms.loc[1, "Compensation Voltage"])) == "0"
    assert float(csms.loc[1, "Combined Score"]) == pytest.approx(24.473)

# check output shape target
def test1_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["TargetLib"]

    assert sl.shape[0] == 55 and sl.shape[1] == 32

# check values target
def test2_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["TargetLib"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "KIEC[Carbamidomethyl]FDSVEISGVEDR_KIEC[Carbamidomethyl]FDSVEISGVEDR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "False"
    assert str(sl.loc[0, "DecoyType"]) == "TT"

    assert str(sl.loc[54, "ModifiedPeptide"]) == "MIAKSEQEIGK_DFQFYKVR"
    assert str(sl.loc[54, "FragmentCharge"]) == "1"
    assert str(sl.loc[54, "FragmentType"]) == "b"
    assert str(sl.loc[54, "FragmentNumber"]) == "5"
    assert str(sl.loc[54, "FragmentPepId"]) == "1"
    assert str(sl.loc[54, "CLContainingFragment"]) == "False"
    assert str(sl.loc[54, "IsDecoy"]) == "False"
    assert str(sl.loc[54, "DecoyType"]) == "TT"

# check output shape decoy dd
def test3_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib"]

    assert sl.shape[0] == 54 and sl.shape[1] == 32

# check values decoy dd
def test4_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "DEVGSIEVSDFC[Carbamidomethyl]EIKR_DEVGSIEVSDFC[Carbamidomethyl]EIKR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "True"
    assert str(sl.loc[0, "DecoyType"]) == "DD"

    assert str(sl.loc[53, "ModifiedPeptide"]) == "GIEQESKAIMK_VKYFQFDR"
    assert str(sl.loc[53, "FragmentCharge"]) == "1"
    assert str(sl.loc[53, "FragmentType"]) == "b"
    assert str(sl.loc[53, "FragmentNumber"]) == "5"
    assert str(sl.loc[53, "FragmentPepId"]) == "1"
    assert str(sl.loc[53, "CLContainingFragment"]) == "True"
    assert str(sl.loc[53, "IsDecoy"]) == "True"
    assert str(sl.loc[53, "DecoyType"]) == "DD"

# check output shape decoy dt
def test5_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib_DT"]

    assert sl.shape[0] == 54 and sl.shape[1] == 32

# check values decoy dt
def test6_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib_DT"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "DEVGSIEVSDFC[Carbamidomethyl]EIKR_KIEC[Carbamidomethyl]FDSVEISGVEDR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "True"
    assert str(sl.loc[0, "DecoyType"]) == "DT"

    assert str(sl.loc[53, "ModifiedPeptide"]) == "GIEQESKAIMK_DFQFYKVR"
    assert str(sl.loc[53, "FragmentCharge"]) == "1"
    assert str(sl.loc[53, "FragmentType"]) == "b"
    assert str(sl.loc[53, "FragmentNumber"]) == "5"
    assert str(sl.loc[53, "FragmentPepId"]) == "1"
    assert str(sl.loc[53, "CLContainingFragment"]) == "False"
    assert str(sl.loc[53, "IsDecoy"]) == "True"
    assert str(sl.loc[53, "DecoyType"]) == "DT"

# check output shape decoy td
def test7_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib_TD"]

    assert sl.shape[0] == 55 and sl.shape[1] == 32

# check values decoy td
def test8_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib_TD"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "KIEC[Carbamidomethyl]FDSVEISGVEDR_DEVGSIEVSDFC[Carbamidomethyl]EIKR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "True"
    assert str(sl.loc[0, "DecoyType"]) == "TD"

    assert str(sl.loc[54, "ModifiedPeptide"]) == "MIAKSEQEIGK_VKYFQFDR"
    assert str(sl.loc[54, "FragmentCharge"]) == "1"
    assert str(sl.loc[54, "FragmentType"]) == "b"
    assert str(sl.loc[54, "FragmentNumber"]) == "5"
    assert str(sl.loc[54, "FragmentPepId"]) == "1"
    assert str(sl.loc[54, "CLContainingFragment"]) == "True"
    assert str(sl.loc[54, "IsDecoy"]) == "True"
    assert str(sl.loc[54, "DecoyType"]) == "TD"

# check output shape full
def test9_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["FullLib"]

    assert sl.shape[0] == 218 and sl.shape[1] == 32
    assert sl["DecoyType"].value_counts()["TT"] == 55
    assert sl["DecoyType"].value_counts()["TD"] == 55
    assert sl["DecoyType"].value_counts()["DD"] == 54
    assert sl["DecoyType"].value_counts()["DT"] == 54

# check values full
def test10_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["FullLib"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "KIEC[Carbamidomethyl]FDSVEISGVEDR_KIEC[Carbamidomethyl]FDSVEISGVEDR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "False"
    assert str(sl.loc[0, "DecoyType"]) == "TT"

    assert str(sl.loc[217, "ModifiedPeptide"]) == "MIAKSEQEIGK_VKYFQFDR"
    assert str(sl.loc[217, "FragmentCharge"]) == "1"
    assert str(sl.loc[217, "FragmentType"]) == "b"
    assert str(sl.loc[217, "FragmentNumber"]) == "5"
    assert str(sl.loc[217, "FragmentPepId"]) == "1"
    assert str(sl.loc[217, "CLContainingFragment"]) == "True"
    assert str(sl.loc[217, "IsDecoy"]) == "True"
    assert str(sl.loc[217, "DecoyType"]) == "TD"
