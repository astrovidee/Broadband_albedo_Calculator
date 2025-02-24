"""
Created on Sat June 5 21:25:32 2021
@author: Vidya Venkatesan (vidyav1@uci.edu)

Project description: This code is a python version of ice_gcm.pro which calculates the albedo from stellar 
spectrum and spectrum of reflectance. 
Note: any spectrum file can be used along with any
reflectance spectrum provided the minGrid and maxGrid parameters are changed to reflect the appropriate units. 
Units are not important (as long as the wavelength units are the same for each file) for the 
final albedo calculations since it is normalized by both the wavelength and the stellar spectrum.
This 6-band calculator can be used for models like ROCKE-3D.
Required Inputs:
    stellarfile = file containing two columns:
                         1 - wavelength, microns
                         2 - flux, in any units (this is normalized
                             for the albedo calculation)
    albedofile  = file containing three columns (as generate by the USGS):
                         1 - wavelength, in microns
                         2 - reflectance at that wavelength (0 - 1)
                         3 - error, if any
"""
# Importing all the libraries
import numpy as np
import scipy.interpolate as interpol
import matplotlib.pyplot as plt
stel_convert_wavl = 1.0e6
# Reading in the stellar file
with open('sun.txt', 'r') as f:
    lines = f.readlines()
    lamda = [float(line.split()[0]) * 10**6 for line in lines ] 
    flux = [float(line.split()[1])  for line in lines]
print(lamda)
# Reading in the Albedo file for surface
with open('croyoconite.txt', 'r') as f:
    lines = f.readlines()
    wave = [float(line.split()[0]) for line in lines]
    #wave = [x * 10**-3 for x in wave]  # Convert from nanometers to microns, if needed
    albedo = [float(line.split()[1]) for line in lines]

# Find the minimum and maximum in both stellar and albedo files and match them for interpolation
start_g = max(np.min(wave), np.min(lamda))
end_g = min(np.max(wave), np.max(lamda))

# Set up the grid
ngrid = 10000
n = np.arange(ngrid)
k = end_g - start_g
wavelengthgrid = start_g + k * n / (ngrid - 1)
dlamda = (np.max(wavelengthgrid) - np.min(wavelengthgrid)) / (len(wavelengthgrid) - 1)

# Interpolation
stellarInterpolate4 = np.interp(wavelengthgrid, lamda, flux)
albedoInterpolate = interpol.CubicSpline(wave, albedo)(wavelengthgrid)

# Normalization
flux4 = stellarInterpolate4 / np.max(stellarInterpolate4)

# More calculations for scaling
Snorm = 1360
total_sed4 = sum(stellarInterpolate4 * dlamda)
ScaleFac = Snorm / total_sed4
new_flux4 = stellarInterpolate4 * ScaleFac

# Define wavelength cutoffs in microns
cutoffs = [0.3, 0.77, 0.86, 1.25, 1.5, 2.2, 4.0]

# Create a dictionary to store albedo results for each range
albedo_results = {}

# Loop through each range defined by the cutoffs
for i in range(len(cutoffs) - 1):
    # Define the start and end of the current range
    range_start = cutoffs[i]
    range_end = cutoffs[i + 1]

    # Get indices for the current range
    indices = np.where((wavelengthgrid >= range_start) & (wavelengthgrid < range_end))

    # Calculate albedo for the current range
    albedo_sed_range = (albedoInterpolate[indices]) * (new_flux4[indices]) * dlamda
    total_albedo_sed_range = sum(albedo_sed_range)
    total_sed_range = sum(new_flux4[indices] * dlamda)
    albedo_range = total_albedo_sed_range / total_sed_range if total_sed_range != 0 else 0

    # Store the result with a descriptive key
    albedo_results[f"Albedo {range_start:.2f} - {range_end:.2f} μm"] = albedo_range

# Calculate albedo for the entire range of wavelengths
indices_all = np.where((wavelengthgrid >= cutoffs[0]) & (wavelengthgrid < cutoffs[-1]))
albedo_sed_all = (albedoInterpolate[indices_all]) * (new_flux4[indices_all]) * dlamda
total_albedo_sed_all = sum(albedo_sed_all)
total_sed_all = sum(new_flux4[indices_all] * dlamda)
albedo_all = total_albedo_sed_all / total_sed_all if total_sed_all != 0 else 0

# Store the total albedo result
albedo_results["Albedo Total (0.3 - 4.0 μm)"] = albedo_all

# Print the results for all ranges
for key, value in albedo_results.items():
    print(f"{key}: {value:.6f}")
