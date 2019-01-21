import numpy as np
from os import system
from tqdm import tqdm


# If FSPS has the wrong settings, switch it back and recompile

with open('/home/adam/Research/fsps/src/sps_vars.f90', 'r') as readfile:
    lines = readfile.readlines()
    if ('MILES 0' in lines[8]) or ('MIST 0' in lines[14]):
        system('fsps_mm')

import fsps



def output_ssp_spectra(taxis = np.arange(0.001,0.011,0.001), logz = 0.0, emissionlines = True):

    starpop = fsps.StellarPopulation(zcontinuous=1, add_neb_emission = emissionlines, logzsol = logz, imf_type = 1, sfh = 0, dust_type = 2)

    spectra = []
    wavelength = []

    for t in tqdm(taxis):

        wavelength, spectrum = starpop.get_spectrum(tage = t)

        spectra.append(spectrum)

    table = np.vstack((wavelength, spectra)).T

    np.savetxt('./synth_spectra/spectra_%.1f_Z_SSPs.txt' % np.exp(logz), table, header = '# lambda, followed by SSP spectra going from 1 to 10 Myr')


def output_csf_spectra(taxis = np.arange(0.001,0.011,0.001), logz = 0.0, sfr = 1., emissionlines = True):

    starpop = fsps.StellarPopulation(zcontinuous=1, add_neb_emission = emissionlines, logzsol = logz, imf_type = 1, sfh = 3, dust_type = 2)
    starpop.set_tabular_sfh(np.array([0,10]), np.array([sfr,sfr]))

    spectra = []
    wavelength = []

    for t in tqdm(taxis):

        wavelength, spectrum = starpop.get_spectrum(tage = t)

        spectra.append(spectrum)

    table = np.vstack((wavelength, spectra)).T

    np.savetxt('./synth_spectra/spectra_%.1f_Z_CSF.txt' % np.exp(logz), table, header = '# lambda, followed by CSF spectra going from 1 to 10 Myr')
