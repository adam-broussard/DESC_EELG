import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from photosim import photosim
from tqdm import tqdm

phot = photosim()

font = {'family':'Roboto', 'weight':'light'}


def sigmadev(redshift = 0.3, sfhtype = 3, metallicity = 1.0, imf_type = 1, timearray = np.array([-0.998, -0.995, -0.990, -0.980, -0.970, -0.960, 
    -0.940, -0.920, -0.900, -0.800, -0.600, -0.200, 0.200, 0.600, 1.000, 1.400, 1.800, 2.200, 
    2.600, 3.000, 3.400, 3.800, 4.200, 4.600, 5.000, 5.400, 5.800, 6.200, 6.600, 7.000, 7.400, 
    7.800, 8.200, 8.600, 9.000, 9.400, 9.800, 10.200, 10.600, 11.000])):

    emline_fnu = []
    noemline_fnu = []
    emline_wave = []
    noemline_wave = []

    # Generate the spectra with and without emission lines

    for (index, time) in enumerate(tqdm(timearray)):

        wave, l_nu = phot.find_spectrum(time, metallicity, imf_type, sfhtype, emline = True)
        new_wave, f_nu = phot.redshift(redshift, wave, l_nu)
        emline_fnu.append(f_nu)
        emline_wave.append(new_wave)

        wave, l_nu = phot.find_spectrum(time, metallicity, imf_type, sfhtype, emline = False)
        new_wave, f_nu = phot.redshift(redshift, wave, l_nu)
        noemline_fnu.append(f_nu)
        noemline_wave.append(new_wave)

    # Run each spectrum through the filters

    phot_diff = []

    for (index, (em_fnu, em_wave, noem_fnu, noem_wave)) in enumerate(zip(emline_fnu, emline_wave, noemline_fnu, noemline_wave)):

        em_phot_wave, em_phot_fnu, em_phot_err = phot.filters.get_photometry(em_wave, em_fnu)
        noem_phot_wave, noem_phot_fnu, noem_phot_err = phot.filters.get_photometry(noem_wave, noem_fnu)

        # em_phot_err and noem_phot_err should be the same, as should the wavelength arrays
        # Calculate the sigma difference between the photometry with and without emission lines

        phot_diff.append((em_phot_fnu - noem_phot_fnu) / em_phot_err)

    # Make phot_diff such that 1st dimension specifies filter, 2nd dimension specifies time

    phot_diff = np.vstack(phot_diff).T

    fig = plt.figure(figsize = (8,8))
    sp = fig.add_subplot(111)

    for (index, band_phot_diff) in enumerate(phot_diff):

        sp.plot(timearray+1., band_phot_diff, linewidth = 2, label = ['u', 'g', 'r', 'i', 'z', 'y'][index])

    sp.legend(loc = 'upper right', fontsize = 18)

    sp.set_ylabel('Sigma Difference', fontsize = 24, fontdict = font)
    sp.set_xlabel('Time (Gyr)', fontsize = 24, fontdict = font)



def sigmadev_diagnostic(redshift = 0.3, sfhtype = 3, metallicity = 1.0, imf_type = 1, timearray = np.array([-0.998, -0.995, -0.990, -0.980, -0.970, -0.960, 
    -0.940, -0.920, -0.900, -0.800, -0.600, -0.200, 0.200, 0.600, 1.000, 1.400, 1.800, 2.200, 
    2.600, 3.000, 3.400, 3.800, 4.200, 4.600, 5.000, 5.400, 5.800, 6.200, 6.600, 7.000, 7.400, 
    7.800, 8.200, 8.600, 9.000, 9.400, 9.800, 10.200, 10.600, 11.000])):

    emline_fnu = []
    noemline_fnu = []
    emline_wave = []
    noemline_wave = []

    # Generate the spectra with and without emission lines

    for (index, time) in enumerate(tqdm(timearray)):

        wave, l_nu = phot.find_spectrum(time, metallicity, imf_type, sfhtype, emline = True)
        new_wave, f_nu = phot.redshift(redshift, wave, l_nu)
        emline_fnu.append(f_nu)
        emline_wave.append(new_wave)

        wave, l_nu = phot.find_spectrum(time, metallicity, imf_type, sfhtype, emline = False)
        new_wave, f_nu = phot.redshift(redshift, wave, l_nu)
        noemline_fnu.append(f_nu)
        noemline_wave.append(new_wave)

    # Run each spectrum through the filters

    emline_phot_wave = []
    emline_phot_fnu = []
    emline_phot_err = []
    noemline_phot_wave = []
    noemline_phot_fnu = []
    noemline_phot_err = []

    phot_diff = []

    for (index, (em_fnu, em_wave, noem_fnu, noem_wave)) in enumerate(zip(emline_fnu, emline_wave, noemline_fnu, noemline_wave)):

        em_phot_wave, em_phot_fnu, em_phot_err = phot.filters.get_photometry(em_wave, em_fnu)
        noem_phot_wave, noem_phot_fnu, noem_phot_err = phot.filters.get_photometry(noem_wave, noem_fnu)

        emline_phot_wave.append(em_phot_wave)
        emline_phot_fnu.append(em_phot_fnu)
        emline_phot_err.append(em_phot_err)
        noemline_phot_wave.append(noem_phot_wave)
        noemline_phot_fnu.append(noem_phot_fnu)
        noemline_phot_err.append(noem_phot_err)

        # em_phot_err and noem_phot_err should be the same, as should the wavelength arrays
        # Calculate the sigma difference between the photometry with and without emission lines

        phot_diff.append((em_phot_fnu - noem_phot_fnu) / em_phot_err)

    # Make phot_diff such that 1st dimension specifies filter, 2nd dimension specifies time

    phot_diff = np.vstack(phot_diff).T
    emline_phot_wave = np.vstack(emline_phot_wave).T
    emline_phot_fnu = np.vstack(emline_phot_fnu).T
    emline_phot_err = np.vstack(emline_phot_err).T
    noemline_phot_wave = np.vstack(noemline_phot_wave).T
    noemline_phot_fnu = np.vstack(noemline_phot_fnu).T
    noemline_phot_err = np.vstack(noemline_phot_err).T

    phot_diff = (emline_phot_fnu - noemline_phot_fnu)/emline_phot_err

    fig = plt.figure(figsize = (8,8))
    sp = fig.add_subplot(111)

    # for (index, (em_wave, em_phot, noem_wave, noem_phot)) in enumerate(zip(emline_phot_wave, emline_phot_fnu, noemline_phot_wave, noemline_phot_fnu)):

        # sp.scatter(em_wave, em_phot, s = 15, edgecolor = 'None', c = 'b')
        # sp.scatter(noem_wave, noem_phot, s = 15, edgecolor = 'None', c = 'r')


    for (index, band_phot_diff) in enumerate(phot_diff):

        sp.plot(timearray+1., band_phot_diff, linewidth = 2, label = ['u', 'g', 'r', 'i', 'z', 'y'][index])

    sp.legend(loc = 'upper right', fontsize = 18)

    # sp.set_xscale('log')
    # sp.set_ylim(0,0.5e-27)

    sp.set_ylabel('Flux Density', fontsize = 24, fontdict = font)
    sp.set_xlabel('Wavelength (A)', fontsize = 24, fontdict = font)