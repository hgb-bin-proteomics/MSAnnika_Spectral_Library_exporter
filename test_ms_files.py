#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.7"
# dependencies = [
#   "pyteomics[XML]",
# ]
# ///

# MS FILE TESTER
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# version tracking
__version = "1.0.0"
__date = "2025-09-02"

# REQUIREMENTS
# pip install pyteomics

import sys
from pyteomics import mgf, mzml
from typing import Dict, List, Tuple, Union, BinaryIO

# parse scan number from pyteomics mgf params
def parse_scannr(params: Dict, pattern: str = "\\.\\d+\\.", i: int = 0) -> Tuple[int, int]:
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

    if isinstance(filename, str) and filename.split(".")[-1].lower() == "mzml":
        print("Reading mzML file...")
        with mzml.read(filename) as reader:
            for spectrum in reader:
                scan_nr = int(str(spectrum["id"]).split("scan=")[1].split()[0])
                ms_level = int(spectrum["ms level"])
                if ms_level != 2:
                    continue
                if (
                    "scanList" not in spectrum
                    or "scan" not in spectrum["scanList"]
                    or len(spectrum["scanList"]["scan"]) != 1
                ):
                    raise RuntimeError(f"Can't get retention time for spectrum: {spectrum}")
                rt_in_min = 0.0
                try:
                    rt_in_min = float(spectrum["scanList"]["scan"][0]["scan start time"])
                except Exception as e:
                    pass
                rt_in_sec = rt_in_min * 60.0
                if "precursorList" not in spectrum:
                    raise RuntimeError(
                        f"[precursorList] No precursor for MS2 spectrum found: {spectrum}"
                    )
                if (
                    "precursor" not in spectrum["precursorList"]
                    or len(spectrum["precursorList"]["precursor"]) != 1
                ):
                    raise RuntimeError(
                        f"[precursor] No precursor for MS2 spectrum found: {spectrum}"
                    )
                for precursor in spectrum["precursorList"]["precursor"]:
                    if "selectedIonList" not in precursor:
                        raise RuntimeError(
                            f"[selectedIonList] No precursor for MS2 spectrum found: {spectrum}"
                        )
                    if (
                        "selectedIon" not in precursor["selectedIonList"]
                        or len(precursor["selectedIonList"]["selectedIon"]) != 1
                    ):
                        raise RuntimeError(
                            f"[selectedIon] No precursor for MS2 spectrum found: {spectrum}"
                        )
                    for ion in precursor["selectedIonList"]["selectedIon"]:
                        spectrum_dict = dict()
                        spectrum_dict["precursor"] = float(ion["selected ion m/z"])
                        spectrum_dict["charge"] = int(ion["charge state"])
                        spectrum_dict["rt"] = rt_in_sec
                        spectrum_dict["max_intensity"] = float(max(spectrum["intensity array"])) if len(spectrum["intensity array"]) > 0 else 0.0
                        peaks = dict()
                        for i, mz in enumerate(spectrum["m/z array"]):
                            peaks[float(mz)] = float(spectrum["intensity array"][i])
                        spectrum_dict["peaks"] = peaks
                        if scan_nr in result_dict:
                            raise RuntimeError(f"Spectrum with scan number {scan_nr} already exists!")
                        result_dict[scan_nr] = spectrum_dict
            reader.close()
    else:
        print("Reading MGF file...")
        with mgf.read(filename) as reader:
            for spectrum in reader:
                parser_result = parse_scannr(spectrum["params"])
                if parser_result[0] != 0:
                    raise RuntimeError(f"Could not parse scan number for spectrum {spectrum}. Please adjust PARSER_PATTERN in the config file!")
                scan_nr = parser_result[1]
                spectrum_dict = dict()
                spectrum_dict["precursor"] = float(spectrum["params"]["pepmass"][0])
                spectrum_dict["charge"] = int(spectrum["params"]["charge"][0])
                spectrum_dict["rt"] = float(spectrum["params"]["rtinseconds"]) if "rtinseconds" in spectrum["params"] else 0.0
                spectrum_dict["max_intensity"] = float(max(spectrum["intensity array"])) if len(spectrum["intensity array"]) > 0 else 0.0
                peaks = dict()
                for i, mz in enumerate(spectrum["m/z array"]):
                    peaks[float(mz)] = float(spectrum["intensity array"][i])
                spectrum_dict["peaks"] = peaks
                result_dict[scan_nr] = spectrum_dict
            reader.close()

    return result_dict

# read multiple spectra files
def read_multiple_spectra(filenames: List[str]) -> None:
    """
    Returns a dictionary that maps filenames to scan numbers to spectra:
    Dict[str -> Dict[int -> Dict["precursor"        -> float
                                 "charge"           -> int
                                 "max_intensity"    -> float
                                 "peaks"            -> Dict[m/z -> intensity]]
    """

    result_dict = dict()
    errors = list()

    for i, filename in enumerate(filenames):
        current_spectra_file = ".".join(filename.split(".")[:-1]).strip()
        try:
            result_dict = read_spectra(filename)
            print(f"INFO: Read all spectra successfully from file {filename}.")
        except Exception as e:
            print(f"Error while reading file: {filename}")
            print("Error details:")
            print(e)
            errors.append(filename)
        print(f"INFO: Read {i + 1}/{len(filenames)} files...")

    if len(errors) > 0:
        print("Found errors in the following file(s):")
        for error in errors:
            print(error)
        print("Exiting spectral library creation...")
        raise RuntimeError("Errors while reading spectra files!")

    print("Read all spectra files successfully!")
    return

def main():
    if len(sys.argv) < 2:
        return
    return read_multiple_spectra(sys.argv[1:])

if __name__ == "__main__":
    main()
