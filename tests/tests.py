#!/usr/bin/env python3

# MS ANNIKA SPECTRAL LIBRARY EXPORTER - TESTS
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# check output shape target
def test1_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["TargetLib"]

    assert sl.shape[0] == 12 and sl.shape[1] == 28

# check values target
def test2_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["TargetLib"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "True"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "b"
    assert str(sl.loc[2, "FragmentNumber"]) == "2"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "False"

    assert str(sl.loc[11, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[11, "FragmentCharge"]) == "1"
    assert str(sl.loc[11, "FragmentType"]) == "y"
    assert str(sl.loc[11, "FragmentNumber"]) == "5"
    assert str(sl.loc[11, "FragmentPepId"]) == "1"
    assert str(sl.loc[11, "CLContainingFragment"]) == "False"

# check output shape decoy dd
def test3_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib"]

    assert sl.shape[0] == 22 and sl.shape[1] == 28

# check values decoy dd
def test4_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "False"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "y"
    assert str(sl.loc[2, "FragmentNumber"]) == "3"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "True"

    assert str(sl.loc[21, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[21, "FragmentCharge"]) == "1"
    assert str(sl.loc[21, "FragmentType"]) == "y"
    assert str(sl.loc[21, "FragmentNumber"]) == "5"
    assert str(sl.loc[21, "FragmentPepId"]) == "1"
    assert str(sl.loc[21, "CLContainingFragment"]) == "True"

# check output shape decoy dt
def test5_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib_DT"]

    assert sl.shape[0] == 17 and sl.shape[1] == 28
    
# check values decoy dt
def test6_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib_DT"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "False"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "y"
    assert str(sl.loc[2, "FragmentNumber"]) == "3"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "True"

    assert str(sl.loc[16, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[16, "FragmentCharge"]) == "1"
    assert str(sl.loc[16, "FragmentType"]) == "y"
    assert str(sl.loc[16, "FragmentNumber"]) == "5"
    assert str(sl.loc[16, "FragmentPepId"]) == "1"
    assert str(sl.loc[16, "CLContainingFragment"]) == "False"
