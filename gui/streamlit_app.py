#!/usr/bin/env python3

# MS ANNIKA SPECTRAL LIBRARY EXPORTER GUI
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

"""
#####################################################
##                                                 ##
##            -- STREAMLIT MAIN APP --             ##
##                                                 ##
#####################################################
"""

import json
import pandas as pd
import streamlit as st
from streamlit_util import *
from create_spectral_library import main as create_spectral_library

@st.cache_data
def dataframe_to_csv_stream(dataframe: pd.DataFrame) -> str:
    text = dataframe.to_csv(index = True).encode("utf-8")
    return text

# main page content
def main_page():

    title = st.title("MS Annika Spectral Library exporter")

    general_description = \
    """
    The MS Annika Spectral Library exporter allows you to generate a spectral library for
    [Spectronaut](https://biognosys.com/software/spectronaut/) from [MS Annika](https://github.com/hgb-bin-proteomics/MSAnnika)
    results.

    You'll need at least two input files:
    - One or more `.mgf` files containing the spectra associated with the CSMs reported by MS Annika
    - The CSMs reported by MS Annika exported in `.xlsx` format from Proteome Discoverer

    More information about the MS Annika Spectral Library exporter can be found here:

    [https://hgb-bin-proteomics.github.io/MSAnnika_Spectral_Library_exporter/](https://hgb-bin-proteomics.github.io/MSAnnika_Spectral_Library_exporter/)
    """
    description = st.markdown(general_description)

    header_1 = st.subheader("Data Import", divider = "rainbow")

    spectra_file = st.file_uploader("Upload one or more spectrum files:",
                                    type = ["mgf"],
                                    accept_multiple_files = True,
                                    help = "Upload one or more spectrum files to be analyzed in .mgf format.")

    csms_file = st.file_uploader("Upload the exported MS Annika CSMs from Proteome Discoverer:",
                                 type = ["xlsx"],
                                 help = "Upload an identification file that contains CSMs of the spectrum file(s) in .xlsx format.")

    header_2 = st.subheader("Parameters", divider = "rainbow")

    run_name = st.text_input("Name of the run [any descriptive text is allowed]:",
                             value = "Run 1",
                             help = "Identifier for the experiment/analysis/run. Any descriptive text is allowed.")

    crosslinker = st.text_input("Name of the crosslinker [as found in the modifications column]:",
                                value = "DSSO",
                                help = "The name of the used cross-linking reagent as it is denoted in the CSMs file. For example 'DSSO'.")

    modifications = st.text_input("Modifications:",
                                  value = \
                                  """{"Oxidation": [15.994915], "Carbamidomethyl": [57.021464], "DSSO": [54.01056, 85.98264, 103.99320]}""",
                                  help = "The post-translational modifications considered for identification. Given in valid .json format " +
                                         "(an example is given). The format is key-value pairs of 'modification name' and 'list of possible " +
                                         "(monoisotopic) modification masses'.")

    ion_types = st.text_input("Expected ion types:",
                              value = "b, y",
                              help = "The fragment ion types that are expected to occur in the mass spectra, delimited by commas.")

    max_charge = st.number_input("Maximum considered precursor ion charge:",
                                 value = 4,
                                 help = "The maximum considered precursor ion charge.")

    match_tolerance = st.number_input("Fragment ion match tolerance in Da:",
                                      value = 0.02,
                                      help = "The tolerance in Dalton considered for matching fragment ions.")

    iRT_m = st.number_input("The slope (m) of the iRT equation:",
                            value = 1.3066,
                            help = "The slope (m) of the iRT equation 'y = m * x + t'.")

    iRT_t = st.number_input("The y-intercept (t) of the iRT equation:",
                            value = 29.502,
                            help = "The y-intercept (t) of the iRT equation 'y = m * x + t'.")

    l1, l2, center_button, r1, r2 = st.columns(5)

    with center_button:
        run_analysis = st.button("Create Spectral Library!")

    if run_analysis:
        if spectra_file is not None and csms_file is not None:
            with st.spinner("MS Annika Spectral Library exporter is running..."):
                with st.expander("Show logging info:"):
                    with st_stdout("info"):
                        try:
                            parsed_modifications = json.loads(modifications)
                            parsed_ion_types = tuple([ion_type.strip() for ion_type in ion_types.split(",")])
                            result = create_spectral_library(spectra_file,
                                                             csms_file,
                                                             run_name,
                                                             crosslinker,
                                                             parsed_modifications,
                                                             parsed_ion_types,
                                                             int(max_charge),
                                                             float(match_tolerance),
                                                             float(iRT_m),
                                                             float(iRT_t),
                                                             True,
                                                             False)
                            st.session_state["result"] = result
                            status_1 = 0
                        except Exception as e:
                            this_e = st.exception(e)
                            status_1 = 1
            if status_1 == 0:
                res_status_1 = st.success("MS Annika Spectral Library exporter finished successfully!")
            else:
                res_status_1 = st.error("MS Annika Spectral Library exporter stopped prematurely! See log for more information!")
        else:
            res_status_1 = st.error("You need to specify a spectrum AND identifications file!")

    if "result" in st.session_state:
        results_header_1 = st.subheader("Generated Spectral Library:", divider = "rainbow")
        spectral_library = st.dataframe(st.session_state["result"])

    if "result" in st.session_state:
        results_header_2 = st.subheader("Download Results", divider = "rainbow")
        csv_1 = st.download_button(label = "Download spectral library!",
                                   data = dataframe_to_csv_stream(st.session_state["result"]),
                                   file_name = ".".join(csms_file.name.split(".")[:-1]) + "_spectralLibrary.csv",
                                   mime = "text/csv",
                                   help = "Download the generated spectral library in .csv format.")

# side bar and main page loader
def main():

    about_str = \
    """
    The MS Annika Spectral Library exporter generates a spectral library for Spectronaut from MS Annika results.
    """

    st.set_page_config(page_title = "MS Annika Spectral Library exporter",
                       page_icon = ":test_tube:",
                       layout = "wide",
                       initial_sidebar_state = "expanded",
                       menu_items = {"Get Help": "https://github.com/hgb-bin-proteomics/MSAnnika_Spectral_Library_exporter/discussions",
                                     "Report a bug": "https://github.com/hgb-bin-proteomics/MSAnnika_Spectral_Library_exporter/issues",
                                     "About": about_str}
                       )

    title = st.sidebar.title("MS Annika Spectral Library exporter")

    logo = st.sidebar.image("https://raw.githubusercontent.com/hgb-bin-proteomics/MSAnnika/master/logo/annika_logo.png",
                            caption = "The MS Annika logo: Two cross-linked peptides forming the letter 'A'.")

    doc = st.sidebar.markdown(about_str)

    citing_str = "**Citing:** If you are using the MS Annika Spectral Library exporter or MS Annika, please cite as described " + \
                 "[here](https://github.com/hgb-bin-proteomics/MSAnnika_Spectral_Library_exporter/tree/master#citing)."
    citing = st.sidebar.markdown(citing_str)

    contact_str = "**Contact:** [Micha Birklbauer](mailto:micha.birklbauer@fh-hagenberg.at), [Viktoria Dorfer](mailto:viktoria.dorfer@fh-hagenberg.at)"
    contact = st.sidebar.markdown(contact_str)

    license_str = "**License:** [MIT License](https://github.com/hgb-bin-proteomics/MSAnnika_Spectral_Library_exporter/blob/master/LICENSE.md)"
    license = st.sidebar.markdown(license_str)

    project_str = "**Project Page:** [GitHub](https://github.com/hgb-bin-proteomics/MSAnnika_Spectral_Library_exporter/)"
    project = st.sidebar.markdown(project_str)

    main_page()

if __name__ == "__main__":
    main()
