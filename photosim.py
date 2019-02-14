import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from glob import glob
from filtersim import filtersim
from astropy.cosmology import Planck15 as cosmo
from os import system
from os.path import isfile
from tqdm import tqdm
from pandas import read_csv


with open('/home/adam/Research/fsps/src/sps_vars.f90', 'r') as readfile:
    lines = readfile.readlines()
    if ('MILES 0' in lines[8]) or ('MIST 0' in lines[14]):
        system('fsps_mm')

import fsps


font = {'family':'Roboto', 'weight':'light'}

Lsun = 3.848e33 # erg/s


class photosim:

    def __init__(self, inputfolder = './adam_synth_spectra/'):

        self.filters = filtersim()


    def redshift(self, redshift, wavelength, l_nu):

        if redshift == 0:
            return wavelength, l_nu
        else:

            new_wavelength = wavelength * (1.+redshift)
            flux = l_nu * (1+redshift) / (4*np.pi*np.power(cosmo.luminosity_distance(redshift).to('cm').value,2))

            return new_wavelength, flux


    def find_spectrum(self, tage, metallicity, imf_type, sfh_type, dust_type = 2, emline = True, increase_ssp = True, delay_csf = True):

        # params = starpop.params.all_params
        # metallicity = 10.**params['logzsol']
        # imf_type = params['imf_type']
        # sfh_type = params['sfh']
        # dust_type = params['dust_type']

        # Cache and retreive spectra

        fname_params = (tage, metallicity, imf_type, sfh_type, dust_type)

        fname = '%.5f_%.2f_%i_%i_%i.spec' % fname_params

        fstub = './cache/'

        if emline:
            fstub = fstub + 'emline/'
        else:
            fstub = fstub + 'noemline/'

        if isfile(fstub + fname):

            waves, l_nu = np.loadtxt(fstub + fname, unpack = True)

        else:

            starpop = fsps.StellarPopulation(zcontinuous=1, add_neb_emission = True, nebemlineinspec = emline, imf_type = imf_type, dust_type = dust_type, 
                sfh = sfh_type, logzsol = np.log10(metallicity), tage = tage) # Calzetti dust
            if sfh_type == 3:
                # Form stars at 1Msun/yr for 1Gyr, then spike to 200Msun/yr
                # starpop.set_tabular_sfh(np.array([0,0.999,1,1.1]), np.array([1,1,10,10]))
                starpop.set_tabular_sfh(np.array([0,13]), np.array([1,1]))
                if delay_csf:
                    tage = tage + 1.
           
            waves, l_nu = starpop.get_spectrum(tage = tage)

            l_nu = l_nu * Lsun

            np.savetxt(fstub + fname, np.vstack((waves, l_nu)).T, fmt = '%i   %.6e')

        if increase_ssp and sfh_type == 0:
            l_nu = l_nu * 10.**7

        return waves, l_nu

            





    def gencat(self, cat_input = './EAZY_runs/cat_input.param', cat_output = './EAZY_runs/cat.dat'):


        gal_id, redshift, age, sfh, metal, imf = read_csv(cat_input,
            header = None, comment = '#', delimiter = '\s+').values.T

        sfh_type = np.zeros(len(sfh))
        imf_type = np.zeros(len(imf))

        imf_type[imf == 'CHABRIER'] = 1
        imf_type[imf == 'KROUPA'] = 2

        sfh_type[sfh == 'CSF'] = 3

        with open(cat_output, 'w') as writefile:

            writefile.write('# ')
            writefile.write('id'.rjust(4))
            writefile.write('  ')

            for x in range(len(self.filters.keys)):

                writefile.write(('f_' + self.filters.keys[x]).ljust(15))
                writefile.write(('e_' + self.filters.keys[x]).ljust(15))

            # writefile.write('z_spec')

            writefile.write('\n')

            for x in tqdm(xrange(len(gal_id))):

                wavelengths, spec_l_nu = self.find_spectrum(age[x], metal[x], imf_type[x], sfh_type[x])
                spec_l_lambda = spec_l_nu * (wavelengths**2.) / 3.e10 

                shifted_wavelengths, spec_flux = self.redshift(redshift[x], wavelengths, spec_l_nu)

                phot_wave, phot_flux, phot_err = self.filters.get_photometry(shifted_wavelengths, spec_flux)

                writefile.write('   %03i' % gal_id[x])
                writefile.write('  ')

                for y in range(len(phot_wave)):

                    writefile.write(('%.6e' % phot_flux[y]).ljust(15))
                    writefile.write(('%.6e' % phot_err[y]).ljust(15))

                # writefile.write('-1.000')
                writefile.write('\n')





    def save_temp(self, sfhtype = 'CSF', metallicity = 1.0, time = 0.001, renorm = True, savefp = '/home/adam/Research/eazy-photoz/templates/AdamTemps/'):

        fname = savefp + sfhtype +'_%iMyr_Z_%.1f.dat' % (time*1000, metallicity)

        if sfhtype == 'SSP':
            starpop = fsps.StellarPopulation(zcontinuous=1, add_neb_emission = True, imf_type = 1, dust_type = 2, 
                sfh = 0, logzsol = np.log10(metallicity)) # chabrier IMF and Calzetti dust
        elif sfhtype == 'CSF':
            starpop = fsps.StellarPopulation(zcontinuous=1, add_neb_emission = True, imf_type = 1, dust_type = 2, 
                sfh = 3, logzsol = np.log10(metallicity))
            starpop.set_tabular_sfh(np.array([0,1]), np.array([1,1]))

        wavelength, l_nu = self.find_spectrum(starpop, time)

        if renorm:
            l_nu = 2 * l_nu / max(l_nu) # Templates in EAZY seem to be normalized so they peak around 2, so maybe this will help?

        with open(fname, 'w') as writefile:

            for x in range(len(wavelength)):
                writefile.write(('%.5e' % wavelength[x]).ljust(20))
                writefile.write('%.5e' % l_nu[x] + '\n')



























    # def spec_test(self):

    #     fig = plt.figure(figsize = (8,8))
    #     sp = fig.add_subplot(111)

    #     sp.plot(self.wavelength, np.log10(self.l_nu[0][0]/Lsun), linewidth = 1.5, color = 'steelblue', label = 'Z=%.1f, t=1' % self.metallicity[0])
    #     sp.plot(self.wavelength, np.log10(self.l_nu[0][1]/Lsun), linewidth = 1.5, color = 'navy', label = 'Z=%.1f, t=2' % self.metallicity[0])
    #     sp.plot(self.wavelength, np.log10(self.l_nu[1][0]/Lsun), linewidth = 1.5, color = 'limegreen', label = 'Z=%.1f, t=1' % self.metallicity[1])
    #     sp.plot(self.wavelength, np.log10(self.l_nu[1][1]/Lsun), linewidth = 1.5, color = 'forestgreen', label = 'Z=%.1f, t=2' % self.metallicity[1])


    #     sp.set_xscale('log')
    #     # sp.set_yscale('log')

    #     sp.set_ylim(-7, -4)
    #     sp.set_xlim(10**3, 10**5)

    #     sp.set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontdict = font, fontsize = 24)
    #     sp.set_ylabel(r'$\log_{10} L_\nu$ $(\mathrm{L_\odot/cm^2/Hz})$', fontdict = font, fontsize = 24)

    #     sp.legend(loc = 'upper left')


    # def spec_test2(self):

    #     # Plots a few spectra from FSPS, but in F_lambda units for comparison to Byler+2017

    #     fig = plt.figure(figsize = (8,8))
    #     sp = fig.add_subplot(111)

    #     sp.plot(self.wavelength, np.log10(self.l_lambda[0][0]/Lsun), linewidth = 1.5, color = 'steelblue', label = 'Z=%.1f, t=1' % self.metallicity[0])
    #     sp.plot(self.wavelength, np.log10(self.l_lambda[0][1]/Lsun), linewidth = 1.5, color = 'navy', label = 'Z=%.1f, t=2' % self.metallicity[0])
    #     sp.plot(self.wavelength, np.log10(self.l_lambda[1][0]/Lsun), linewidth = 1.5, color = 'limegreen', label = 'Z=%.1f, t=1' % self.metallicity[1])
    #     sp.plot(self.wavelength, np.log10(self.l_lambda[1][1]/Lsun), linewidth = 1.5, color = 'forestgreen', label = 'Z=%.1f, t=2' % self.metallicity[1])


    #     sp.set_xscale('log')
    #     # sp.set_yscale('log')

    #     sp.set_ylim(-6, 0)
    #     sp.set_xlim(10**3, 10**5)

    #     sp.set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontdict = font, fontsize = 24)
    #     sp.set_ylabel(r'$\log_{10} L_\lambda$ $(\mathrm{L_\odot/cm^2/\AA})$', fontdict = font, fontsize = 24)

    #     sp.legend(loc = 'upper right')


    # def spec_test_all(self):

    #     fig = plt.figure(figsize = (8,8))
    #     sp = fig.add_subplot(111)

    #     x = 1

    #     colors = [mpl.cm.jet(int(place)) for place in np.linspace(0, mpl.cm.jet.N, len(self.l_lambda[x]))]

    #     # for x in range(len(self.l_lambda)):
    #     for y in range(0,len(self.l_lambda[x]),2):

    #         sp.plot(self.wavelength, np.log10(self.l_lambda[x][y]/Lsun), linewidth = 1.5, color = colors[y], label = '%i Myr' % (y+1))


    #     sp.set_xscale('log')
    #     # sp.set_yscale('log')

    #     sp.set_ylim(-5.5, 0)
    #     sp.set_xlim(10**3, 10**5)

    #     sp.set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontdict = font, fontsize = 24)
    #     sp.set_ylabel(r'$\log_{10} L_\nu$ $(\mathrm{L_\odot/cm^2/Hz})$', fontdict = font, fontsize = 24)

    #     sp.legend(loc = 'upper right')




    # def phot_test(self, redshift = 0., oldphot = False):

    #     fig = plt.figure(figsize = (8,8))
    #     sp = fig.add_subplot(111)

    #     wavelength, spec_flux = self.redshift(redshift, self.wavelength, self.l_nu[0][0])

    #     if oldphot:
    #         phot_wave, phot_flux, phot_err = self.filters.get_photometry_old(wavelength, spec_flux)
            
    #     else:
    #         phot_wave, phot_flux, phot_err = self.filters.get_photometry(wavelength, spec_flux)

    #     colors = [mpl.cm.inferno(int(x)) for x in np.linspace(0, mpl.cm.viridis.N, len(self.filters.keys))]

    #     sp.plot(wavelength, np.log10(spec_flux), linewidth = 2, color = 'k')

    #     for x in range(len(phot_wave)):
    #         sp.scatter(phot_wave[x], np.log10(phot_flux[x]), s = 40, color = colors[x], edgecolors = 'k', label = self.filters.keys[x])

    #     # sp.set_xscale('log')
    #     # sp.set_yscale('log')

    #     # sp.set_ylim(19, 23.5)
    #     sp.set_xlim(3000, 12000)
    #     sp.set_ylim(25,35)

    #     sp.set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontdict = font, fontsize = 24)
    #     sp.set_ylabel(r'$\log_{10}L_\nu$ $(\mathrm{erg/s/Hz})$', fontdict = font, fontsize = 24)

    #     sp.legend(loc = 'upper left')



    # def plot_metallicity(self, metallicity = 0):


    #     fig = plt.figure(figsize = (8,8))
    #     sp = fig.add_subplot(111)

    #     colors = [mpl.cm.viridis(int(x)) for x in np.linspace(0, mpl.cm.viridis.N, len(self.l_nu[metallicity]))]

    #     for x in range(len(self.l_nu[metallicity])):

    #         sp.plot(self.wavelength, self.l_nu[metallicity][x], color = colors[x], label = 'Z=%.1f, t=%i' % (self.metallicity[metallicity], x+1))


    #     sp.set_xscale('log')
    #     sp.set_yscale('log')

    #     sp.set_ylim(10**-19, 10**-9)
    #     sp.set_xlim(10**2, 10**6)

    #     sp.legend(loc = 'upper left')





