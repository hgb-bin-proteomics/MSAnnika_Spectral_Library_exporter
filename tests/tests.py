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

    assert sl.shape[0] == 12 and sl.shape[1] == 32

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
    assert str(sl.loc[0, "IsDecoy"]) == "False"
    assert str(sl.loc[0, "DecoyType"]) == "TT"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "True"
    assert str(sl.loc[1, "IsDecoy"]) == "False"
    assert str(sl.loc[1, "DecoyType"]) == "TT"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "b"
    assert str(sl.loc[2, "FragmentNumber"]) == "2"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"
    assert str(sl.loc[2, "IsDecoy"]) == "False"
    assert str(sl.loc[2, "DecoyType"]) == "TT"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "False"
    assert str(sl.loc[3, "IsDecoy"]) == "False"
    assert str(sl.loc[3, "DecoyType"]) == "TT"

    assert str(sl.loc[11, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[11, "FragmentCharge"]) == "1"
    assert str(sl.loc[11, "FragmentType"]) == "y"
    assert str(sl.loc[11, "FragmentNumber"]) == "5"
    assert str(sl.loc[11, "FragmentPepId"]) == "1"
    assert str(sl.loc[11, "CLContainingFragment"]) == "False"
    assert str(sl.loc[11, "IsDecoy"]) == "False"
    assert str(sl.loc[11, "DecoyType"]) == "TT"

# check output shape decoy dd
def test3_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib"]

    assert sl.shape[0] == 22 and sl.shape[1] == 32

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
    assert str(sl.loc[0, "IsDecoy"]) == "True"
    assert str(sl.loc[0, "DecoyType"]) == "DD"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "False"
    assert str(sl.loc[1, "IsDecoy"]) == "True"
    assert str(sl.loc[1, "DecoyType"]) == "DD"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "y"
    assert str(sl.loc[2, "FragmentNumber"]) == "3"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"
    assert str(sl.loc[2, "IsDecoy"]) == "True"
    assert str(sl.loc[2, "DecoyType"]) == "DD"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "True"
    assert str(sl.loc[3, "IsDecoy"]) == "True"
    assert str(sl.loc[3, "DecoyType"]) == "DD"

    assert str(sl.loc[21, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[21, "FragmentCharge"]) == "1"
    assert str(sl.loc[21, "FragmentType"]) == "y"
    assert str(sl.loc[21, "FragmentNumber"]) == "5"
    assert str(sl.loc[21, "FragmentPepId"]) == "1"
    assert str(sl.loc[21, "CLContainingFragment"]) == "True"
    assert str(sl.loc[21, "IsDecoy"]) == "True"
    assert str(sl.loc[21, "DecoyType"]) == "DD"

# check output shape decoy dt
def test5_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib_DT"]

    assert sl.shape[0] == 17 and sl.shape[1] == 32

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
    assert str(sl.loc[0, "IsDecoy"]) == "True"
    assert str(sl.loc[0, "DecoyType"]) == "DT"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "False"
    assert str(sl.loc[1, "IsDecoy"]) == "True"
    assert str(sl.loc[1, "DecoyType"]) == "DT"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "y"
    assert str(sl.loc[2, "FragmentNumber"]) == "3"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"
    assert str(sl.loc[2, "IsDecoy"]) == "True"
    assert str(sl.loc[2, "DecoyType"]) == "DT"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "True"
    assert str(sl.loc[3, "IsDecoy"]) == "True"
    assert str(sl.loc[3, "DecoyType"]) == "DT"

    assert str(sl.loc[16, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[16, "FragmentCharge"]) == "1"
    assert str(sl.loc[16, "FragmentType"]) == "y"
    assert str(sl.loc[16, "FragmentNumber"]) == "5"
    assert str(sl.loc[16, "FragmentPepId"]) == "1"
    assert str(sl.loc[16, "CLContainingFragment"]) == "False"
    assert str(sl.loc[16, "IsDecoy"]) == "True"
    assert str(sl.loc[16, "DecoyType"]) == "DT"

# check output shape decoy td
def test7_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib_TD"]

    assert sl.shape[0] == 17 and sl.shape[1] == 32

# check values decoy td
def test8_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["DecoyLib_TD"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "True"
    assert str(sl.loc[0, "DecoyType"]) == "TD"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "True"
    assert str(sl.loc[1, "IsDecoy"]) == "True"
    assert str(sl.loc[1, "DecoyType"]) == "TD"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "b"
    assert str(sl.loc[2, "FragmentNumber"]) == "2"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"
    assert str(sl.loc[2, "IsDecoy"]) == "True"
    assert str(sl.loc[2, "DecoyType"]) == "TD"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "False"
    assert str(sl.loc[3, "IsDecoy"]) == "True"
    assert str(sl.loc[3, "DecoyType"]) == "TD"

    assert str(sl.loc[16, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[16, "FragmentCharge"]) == "1"
    assert str(sl.loc[16, "FragmentType"]) == "y"
    assert str(sl.loc[16, "FragmentNumber"]) == "5"
    assert str(sl.loc[16, "FragmentPepId"]) == "1"
    assert str(sl.loc[16, "CLContainingFragment"]) == "True"
    assert str(sl.loc[16, "IsDecoy"]) == "True"
    assert str(sl.loc[16, "DecoyType"]) == "TD"

# check output shape full
def test9_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["FullLib"]

    assert sl.shape[0] == 68 and sl.shape[1] == 32

# check values full
def test10_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()
    sl = sl["FullLib"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "False"
    assert str(sl.loc[0, "DecoyType"]) == "TT"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "True"
    assert str(sl.loc[1, "IsDecoy"]) == "False"
    assert str(sl.loc[1, "DecoyType"]) == "TT"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "b"
    assert str(sl.loc[2, "FragmentNumber"]) == "2"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"
    assert str(sl.loc[2, "IsDecoy"]) == "False"
    assert str(sl.loc[2, "DecoyType"]) == "TT"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "False"
    assert str(sl.loc[3, "IsDecoy"]) == "False"
    assert str(sl.loc[3, "DecoyType"]) == "TT"

    assert str(sl.loc[67, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[67, "FragmentCharge"]) == "1"
    assert str(sl.loc[67, "FragmentType"]) == "y"
    assert str(sl.loc[67, "FragmentNumber"]) == "5"
    assert str(sl.loc[67, "FragmentPepId"]) == "1"
    assert str(sl.loc[67, "CLContainingFragment"]) == "True"
    assert str(sl.loc[67, "IsDecoy"]) == "True"
    assert str(sl.loc[67, "DecoyType"]) == "TD"

# check reverse modification annotator
def test11_spectral_library_exporter():

    import pandas as pd
    from create_spectral_library import generate_decoy_csm_dd
    from create_spectral_library import get_ModifiedPeptide

    data = pd.read_excel("test_reverse_mods.xlsx")
    decoy = generate_decoy_csm_dd(data.iloc[0])

    assert str(get_ModifiedPeptide(decoy)) == "EPSYEKWKQEK_QSAIQSMMLEKATQVRVGGGTVDRLNM[Oxidation]M[Oxidation]SR"

# check filtering
def test12_spectral_library_exporter():

    import pytest
    import pandas as pd
    from create_spectral_library import filter_df_for_unique_residue_pairs

    data = pd.read_excel("test_filter.xlsx")
    filtered_data = filter_df_for_unique_residue_pairs(data)

    assert filtered_data.shape[0] == 2
    checked = 0
    for i, row in filtered_data.iterrows():
        if str(row["Modifications A"]).strip() == "K1(DSSO)":
            assert float(row["Combined Score"]) == pytest.approx(198.71)
            checked += 1
        if str(row["Modifications A"]).strip() == "K1(DSSO);M6(Oxidation)":
            assert float(row["Combined Score"]) == pytest.approx(8.71)
            checked += 1
    assert checked == 2

# check kmers calculation
def test13_test_kmers():

    from post_process import get_kmers

    unique_seq_positions = {1,2,3,7,8,11,10,15,16,17,18}
    assert get_kmers(unique_seq_positions) == [3,2,2,4]
    unique_seq_positions = {1,3,5}
    assert get_kmers(unique_seq_positions) == []
    unique_seq_positions = {0,1}
    assert get_kmers(unique_seq_positions) == [2]
    unique_seq_positions = {0,1,3,7,9}
    assert get_kmers(unique_seq_positions) == [2]
    unique_seq_positions = {0,1,3,7,8,9,15}
    assert get_kmers(unique_seq_positions) == [2,3]

## MZML

# check output shape target
def test14_spectral_library_exporter():

    from create_spectral_library import main

    sl = main(spectra_file=["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mzML"])
    sl = sl["TargetLib"]

    assert sl.shape[0] == 12 and sl.shape[1] == 32

# check values target
def test15_spectral_library_exporter():

    from create_spectral_library import main

    sl = main(spectra_file=["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mzML"])
    sl = sl["TargetLib"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "False"
    assert str(sl.loc[0, "DecoyType"]) == "TT"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "True"
    assert str(sl.loc[1, "IsDecoy"]) == "False"
    assert str(sl.loc[1, "DecoyType"]) == "TT"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "b"
    assert str(sl.loc[2, "FragmentNumber"]) == "2"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"
    assert str(sl.loc[2, "IsDecoy"]) == "False"
    assert str(sl.loc[2, "DecoyType"]) == "TT"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "False"
    assert str(sl.loc[3, "IsDecoy"]) == "False"
    assert str(sl.loc[3, "DecoyType"]) == "TT"

    assert str(sl.loc[11, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[11, "FragmentCharge"]) == "1"
    assert str(sl.loc[11, "FragmentType"]) == "y"
    assert str(sl.loc[11, "FragmentNumber"]) == "5"
    assert str(sl.loc[11, "FragmentPepId"]) == "1"
    assert str(sl.loc[11, "CLContainingFragment"]) == "False"
    assert str(sl.loc[11, "IsDecoy"]) == "False"
    assert str(sl.loc[11, "DecoyType"]) == "TT"

# check output shape decoy dd
def test16_spectral_library_exporter():

    from create_spectral_library import main

    sl = main(spectra_file=["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mzML"])
    sl = sl["DecoyLib"]

    assert sl.shape[0] == 22 and sl.shape[1] == 32

# check values decoy dd
def test17_spectral_library_exporter():

    from create_spectral_library import main

    sl = main(spectra_file=["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mzML"])
    sl = sl["DecoyLib"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "True"
    assert str(sl.loc[0, "DecoyType"]) == "DD"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "False"
    assert str(sl.loc[1, "IsDecoy"]) == "True"
    assert str(sl.loc[1, "DecoyType"]) == "DD"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "y"
    assert str(sl.loc[2, "FragmentNumber"]) == "3"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"
    assert str(sl.loc[2, "IsDecoy"]) == "True"
    assert str(sl.loc[2, "DecoyType"]) == "DD"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "True"
    assert str(sl.loc[3, "IsDecoy"]) == "True"
    assert str(sl.loc[3, "DecoyType"]) == "DD"

    assert str(sl.loc[21, "ModifiedPeptide"]) == "HGQQKR_HGQQKR"
    assert str(sl.loc[21, "FragmentCharge"]) == "1"
    assert str(sl.loc[21, "FragmentType"]) == "y"
    assert str(sl.loc[21, "FragmentNumber"]) == "5"
    assert str(sl.loc[21, "FragmentPepId"]) == "1"
    assert str(sl.loc[21, "CLContainingFragment"]) == "True"
    assert str(sl.loc[21, "IsDecoy"]) == "True"
    assert str(sl.loc[21, "DecoyType"]) == "DD"

# check output shape decoy dt
def test18_spectral_library_exporter():

    from create_spectral_library import main

    sl = main(spectra_file=["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mzML"])
    sl = sl["DecoyLib_DT"]

    assert sl.shape[0] == 17 and sl.shape[1] == 32

# check values decoy dt
def test19_spectral_library_exporter():

    from create_spectral_library import main

    sl = main(spectra_file=["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mzML"])
    sl = sl["DecoyLib_DT"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "True"
    assert str(sl.loc[0, "DecoyType"]) == "DT"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "False"
    assert str(sl.loc[1, "IsDecoy"]) == "True"
    assert str(sl.loc[1, "DecoyType"]) == "DT"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "y"
    assert str(sl.loc[2, "FragmentNumber"]) == "3"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"
    assert str(sl.loc[2, "IsDecoy"]) == "True"
    assert str(sl.loc[2, "DecoyType"]) == "DT"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "True"
    assert str(sl.loc[3, "IsDecoy"]) == "True"
    assert str(sl.loc[3, "DecoyType"]) == "DT"

    assert str(sl.loc[16, "ModifiedPeptide"]) == "HGQQKR_KQQGHR"
    assert str(sl.loc[16, "FragmentCharge"]) == "1"
    assert str(sl.loc[16, "FragmentType"]) == "y"
    assert str(sl.loc[16, "FragmentNumber"]) == "5"
    assert str(sl.loc[16, "FragmentPepId"]) == "1"
    assert str(sl.loc[16, "CLContainingFragment"]) == "False"
    assert str(sl.loc[16, "IsDecoy"]) == "True"
    assert str(sl.loc[16, "DecoyType"]) == "DT"

# check output shape decoy td
def test20_spectral_library_exporter():

    from create_spectral_library import main

    sl = main(spectra_file=["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mzML"])
    sl = sl["DecoyLib_TD"]

    assert sl.shape[0] == 17 and sl.shape[1] == 32

# check values decoy td
def test21_spectral_library_exporter():

    from create_spectral_library import main

    sl = main(spectra_file=["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mzML"])
    sl = sl["DecoyLib_TD"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "True"
    assert str(sl.loc[0, "DecoyType"]) == "TD"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "True"
    assert str(sl.loc[1, "IsDecoy"]) == "True"
    assert str(sl.loc[1, "DecoyType"]) == "TD"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "b"
    assert str(sl.loc[2, "FragmentNumber"]) == "2"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"
    assert str(sl.loc[2, "IsDecoy"]) == "True"
    assert str(sl.loc[2, "DecoyType"]) == "TD"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "False"
    assert str(sl.loc[3, "IsDecoy"]) == "True"
    assert str(sl.loc[3, "DecoyType"]) == "TD"

    assert str(sl.loc[16, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[16, "FragmentCharge"]) == "1"
    assert str(sl.loc[16, "FragmentType"]) == "y"
    assert str(sl.loc[16, "FragmentNumber"]) == "5"
    assert str(sl.loc[16, "FragmentPepId"]) == "1"
    assert str(sl.loc[16, "CLContainingFragment"]) == "True"
    assert str(sl.loc[16, "IsDecoy"]) == "True"
    assert str(sl.loc[16, "DecoyType"]) == "TD"

# check output shape full
def test22_spectral_library_exporter():

    from create_spectral_library import main

    sl = main(spectra_file=["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mzML"])
    sl = sl["FullLib"]

    assert sl.shape[0] == 68 and sl.shape[1] == 32

# check values full
def test23_spectral_library_exporter():

    from create_spectral_library import main

    sl = main(spectra_file=["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mzML"])
    sl = sl["FullLib"]

    assert str(sl.loc[0, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[0, "FragmentCharge"]) == "1"
    assert str(sl.loc[0, "FragmentType"]) == "y"
    assert str(sl.loc[0, "FragmentNumber"]) == "1"
    assert str(sl.loc[0, "FragmentPepId"]) == "0"
    assert str(sl.loc[0, "CLContainingFragment"]) == "False"
    assert str(sl.loc[0, "IsDecoy"]) == "False"
    assert str(sl.loc[0, "DecoyType"]) == "TT"

    assert str(sl.loc[1, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[1, "FragmentCharge"]) == "1"
    assert str(sl.loc[1, "FragmentType"]) == "b"
    assert str(sl.loc[1, "FragmentNumber"]) == "2"
    assert str(sl.loc[1, "FragmentPepId"]) == "0"
    assert str(sl.loc[1, "CLContainingFragment"]) == "True"
    assert str(sl.loc[1, "IsDecoy"]) == "False"
    assert str(sl.loc[1, "DecoyType"]) == "TT"

    assert str(sl.loc[2, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[2, "FragmentCharge"]) == "1"
    assert str(sl.loc[2, "FragmentType"]) == "b"
    assert str(sl.loc[2, "FragmentNumber"]) == "2"
    assert str(sl.loc[2, "FragmentPepId"]) == "0"
    assert str(sl.loc[2, "CLContainingFragment"]) == "True"
    assert str(sl.loc[2, "IsDecoy"]) == "False"
    assert str(sl.loc[2, "DecoyType"]) == "TT"

    assert str(sl.loc[3, "ModifiedPeptide"]) == "KQQGHR_KQQGHR"
    assert str(sl.loc[3, "FragmentCharge"]) == "1"
    assert str(sl.loc[3, "FragmentType"]) == "y"
    assert str(sl.loc[3, "FragmentNumber"]) == "3"
    assert str(sl.loc[3, "FragmentPepId"]) == "0"
    assert str(sl.loc[3, "CLContainingFragment"]) == "False"
    assert str(sl.loc[3, "IsDecoy"]) == "False"
    assert str(sl.loc[3, "DecoyType"]) == "TT"

    assert str(sl.loc[67, "ModifiedPeptide"]) == "KQQGHR_HGQQKR"
    assert str(sl.loc[67, "FragmentCharge"]) == "1"
    assert str(sl.loc[67, "FragmentType"]) == "y"
    assert str(sl.loc[67, "FragmentNumber"]) == "5"
    assert str(sl.loc[67, "FragmentPepId"]) == "1"
    assert str(sl.loc[67, "CLContainingFragment"]) == "True"
    assert str(sl.loc[67, "IsDecoy"]) == "True"
    assert str(sl.loc[67, "DecoyType"]) == "TD"

def test24_spectral_library_exporter():

    from create_spectral_library import main

    with pytest.raises(RuntimeError):
        sl = main(spectra_file=["20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_002.mzML"])
