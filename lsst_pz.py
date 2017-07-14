import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from glob import glob

# import FSPS - this lets us generate galaxy spectral energy distributions 
# for a wide range of parameter values

import fsps

# import speclite - since FSPS doesn't have the LSST filter curves yet, 
# use this package to multiply spectra by the filter curves to get SEDs

import speclite

import speclite.filters as filts

centers_lsst = np.array([3664.37, 4807.02, 6209.82, 7542.84, 8700.52, 9633.00])

noNewLine = '\x1b[1A\x1b[1M'

class lsst_filts:

    def __init__(self, z = 1.0, age = 0.0001, filter_fp = '../Tools/speclite/speclite/data/filters/lsst2016-*'):

        filterfiles = glob(filter_fp)

        self.filter_wave = []
        self.filter_resp = []

        for x in filterfiles:

            tempfile = np.transpose(np.loadtxt(x, skiprows = 18))

            self.filter_wave.append(tempfile[0]*10.)

            self.filter_resp.append(tempfile[1])

        self.filters = filts.load_filters('lsst2016-*')

        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

        self.gen_sed(z = z, age = age)



    def gen_sed(self, z = 1.0, age = 0.0001):

        print 'Generating stellar population...'

        self.starpop = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,sfh=0, imf_type=1, logzsol=0.0, dust_type=2, dust1=0.0, dust2=0.0, add_neb_emission=True,add_neb_continuum=True)

        print noNewLine + 'Generating spectrum...'

        lam, spec = self.starpop.get_spectrum(tage = age)

        muJ = self.convert_to_microjansky(spec, z)

        return lam, muJ



    def convert_to_microjansky(self, spec, z):
        return spec *1e6 * 1e23*3.48e33/(4*np.pi*3.086e+24*3.086e+24*self.cosmo.luminosity_distance(z).value*self.cosmo.luminosity_distance(z).value)

    def calc_sed_lsst(self, mags):
        sed = np.zeros((6,))   
        for i in range(6):
            # sed[i] = np.power(10,(mags[0][i] + 48.600)*(-2/5))/(3.34e4*centers_lsst[i]*centers_lsst[i]/1e23)
            sed[i] = 10**(-.4*mags[0][i]-8.9)*10**6
        return sed


    def filtersample(self, inputwave, inputspec, filter_wave, filter_resp):

        norm_filter = filter_resp / np.trapz(filter_resp, x = filter_wave)

        filter_convolved = norm_filter * np.interp(filter_wave, inputwave, inputspec)

        return np.trapz(filter_convolved, x = filter_wave)


    def lsst_sample(self, inputwave, inputspec):

        measurement = []

        for x in xrange(len(self.filter_resp)):

            measurement.append(self.filtersample(inputwave, inputspec, self.filter_wave[x], self.filter_resp[x]))

        return np.array(measurement)
        

    def plot_spec(self, subplot, z):

        pass



    def plot_all_spec(self, zmin, zmax, zstep, saveprefix = './Figures/z_'):

        fig = plt.figure(figsize = (8,8))

        subplot = fig.add_subplot(111)

        for z in arange(zmin, zmax, zstep):

            self.plot_spec(subplot, z)

            save(saveprefix + str(round(z, 2)) + '.png')

            subplot.cla()

        

