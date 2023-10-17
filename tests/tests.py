#!/usr/bin/env python3

# MS ANNIKA SPECTRAL LIBRARY EXPORTER - TESTS
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# check output shape
def test1_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()

    assert sl.shape[0] == 12 and sl.shape[1] == 26

# check values
def test2_spectral_library_exporter():

    from create_spectral_library import main

    sl = main()

    value_tests = True

    if str(sl.loc[0, "FragmentCharge"]) != "1":
        value_tests = False
    if str(sl.loc[0, "FragmentType"]) != "y":
        value_tests = False
    if str(sl.loc[0, "FragmentNumber"]) != "1":
        value_tests = False
    if str(sl.loc[0, "FragmentPepId"]) != "0":
        value_tests = False
    if str(sl.loc[0, "CLContainingFragment"]) != "False":
        value_tests = False

    if str(sl.loc[1, "FragmentCharge"]) != "1":
        value_tests = False
    if str(sl.loc[1, "FragmentType"]) != "b":
        value_tests = False
    if str(sl.loc[1, "FragmentNumber"]) != "2":
        value_tests = False
    if str(sl.loc[1, "FragmentPepId"]) != "0":
        value_tests = False
    if str(sl.loc[1, "CLContainingFragment"]) != "True":
        value_tests = False

    if str(sl.loc[2, "FragmentCharge"]) != "1":
        value_tests = False
    if str(sl.loc[2, "FragmentType"]) != "b":
        value_tests = False
    if str(sl.loc[2, "FragmentNumber"]) != "2":
        value_tests = False
    if str(sl.loc[2, "FragmentPepId"]) != "0":
        value_tests = False
    if str(sl.loc[2, "CLContainingFragment"]) != "True":
        value_tests = False

    if str(sl.loc[3, "FragmentCharge"]) != "1":
        value_tests = False
    if str(sl.loc[3, "FragmentType"]) != "y":
        value_tests = False
    if str(sl.loc[3, "FragmentNumber"]) != "3":
        value_tests = False
    if str(sl.loc[3, "FragmentPepId"]) != "0":
        value_tests = False
    if str(sl.loc[3, "CLContainingFragment"]) != "False":
        value_tests = False

    if str(sl.loc[11, "FragmentCharge"]) != "1":
        value_tests = False
    if str(sl.loc[11, "FragmentType"]) != "y":
        value_tests = False
    if str(sl.loc[11, "FragmentNumber"]) != "5":
        value_tests = False
    if str(sl.loc[11, "FragmentPepId"]) != "1":
        value_tests = False
    if str(sl.loc[11, "CLContainingFragment"]) != "False":
        value_tests = False

    assert value_tests
