# Broadband Albedo Calculator

[![Documentation](https://img.shields.io/badge/Documentation-blue)](https://astrovidee.github.io/Broadband_albedo_Calculator/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16813647.svg)](https://doi.org/10.5281/zenodo.16813647)

Computes the broadband albedo of a surface from a stellar spectrum and a
surface reflectance spectrum. The results can be used as input to energy
balance models and global climate models. This is a Python version of the IDL
routine `ice_gcm.pro`.

![](BB_Calc.png)

## Requirements

Python 3 with NumPy, SciPy and Matplotlib.

## Installation

```bash
git clone https://github.com/astrovidee/Broadband_albedo_Calculator.git
cd Broadband_albedo_Calculator

```

## Usage

```bash
python Broadband_albedo_calculator_2_band.py   # two bands, split at 0.7 microns
python broadband_albedo_calculator_6_band.py   # six bands
```

The two-band version suits 1-D energy balance models and ExoCAM. The six-band
version suits ROCKE-3D. Each script prints the albedo in every band.

## Inputs

Both scripts read two text files from the working directory.

| File | Columns | Example included |
|---|---|---|
| Stellar spectrum | wavelength in meters, flux in any units | `sun.txt` |
| Surface reflectance | wavelength in microns, reflectance from 0 to 1 | `croyoconite.txt` |

To use your own spectra, replace these files or edit the two `open(...)` lines
near the top of the script. The flux is normalized, so its units do not matter.

## How to cite

If you use this code, please cite the paper and the software.

Venkatesan, V., Shields, A. L., Deitrick, R., Wolf, E. T., & Rushby, A. (2025).
A One-Dimensional Energy Balance Model Parameterization for the Formation of
CO2 Ice on the Surfaces of Eccentric Extrasolar Planets. *Astrobiology*, 25(1),
42–59. https://doi.org/10.1089/ast.2023.0103 (arXiv:2501.11667)

Venkatesan, V. (2025). Broadband Albedo Calculator (v1.0.0). Zenodo.
https://doi.org/10.5281/zenodo.16813647

## License and contact

MIT License. Questions and bug reports are welcome through GitHub issues.
