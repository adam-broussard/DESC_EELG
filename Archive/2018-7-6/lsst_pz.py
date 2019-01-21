import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from glob import glob

# import FSPS - this lets us generate galaxy spectral energy distributions 
# for a wide range of parameter values

# import fsps

# import speclite - since FSPS doesn't have the LSST filter curves yet, 
# use this package to multiply spectra by the filter curves to get SEDs

# import speclite

# import speclite.filters as filts

# ADD FILTER WIDTHS AS BARS IN PLOTS
# FIX FILTER CONVOLUTIONS

# filter_centers = np.array([3664.37, 4807.02, 6209.82, 7542.84, 8700.52, 9633.00])
filter_centers = np.array([4807.02, 8700.52, 6209.82, 9633.00, 3664.37, 7542.84])
filter_edges = [[4094.0,5520.0], [8190.0,9211.0], [5531.0,6889.0], [9206.0,10060.0], [3353.0,3975.0], [6920.0,8166.0]]
filter_half_width = np.array([ 712.98,  510.48,  679.18,  427.  ,  310.63,  623.16])

noNewLine = '\x1b[1A\x1b[1M'
font = {'family':'Roboto', 'weight':'light'}

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

        self.lam, self.spec = self.starpop.get_spectrum(tage = age)

        self.muJ = self.convert_to_microjansky(self.spec, z)



    def convert_to_microjansky(self, spec, z):
        return spec *1e6 * 1e23*3.48e33/(4*np.pi*3.086e+24*3.086e+24*self.cosmo.luminosity_distance(z).value*self.cosmo.luminosity_distance(z).value)

    def calc_sed_lsst(self, mags):
        sed = np.zeros((6,))   
        for i in range(6):
            # sed[i] = np.power(10,(mags[0][i] + 48.600)*(-2/5))/(3.34e4*filter_centers[i]*filter_centers[i]/1e23)
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

        ydata = self.lsst_sample(self.lam * (1+z), self.muJ)

        subplot.plot(self.lam * (1+z), self.muJ, label = 'Spectrum', zorder = 1)

        for y in xrange(len(filter_half_width)):
            subplot.errorbar(filter_centers[y], ydata[y], xerr = filter_half_width[y], fmt = 'o')
        # subplot.scatter(filter_centers, self.lsst_sample(self.lam * (1+z), self.muJ), c = 'r', s = 50, zorder = 3)

        subplot.set_xlabel('Wavelength [$\AA$]', fontsize = 24, fontdict = font)
        subplot.set_ylabel('Flux [$\mu$J]', fontsize = 24, fontdict = font)
        subplot.set_title('LSST Filters', fontsize = 28, fontdict = font)

        subplot.text(0.02, 0.98, r'$z = ' + str(round(z,2)) + '$', fontsize = 24, family = 'Roboto', weight = 'light', transform = subplot.transAxes, bbox = dict(facecolor = 'w', edgecolor = 'k'))

        subplot.set_xlim(0, 12000)
        subplot.set_ylim(0, 4.5*10**-8)




    def plot_all_spec(self, zmin, zmax, zstep, saveprefix = './Figures/z_'):

        fig = plt.figure(figsize = (8,8))

        subplot = fig.add_subplot(111)

        for z in np.arange(zmin, zmax, zstep):

            self.plot_spec(subplot, z)

            plt.savefig(saveprefix + str(round(z, 2)).ljust(4, '0') + '.png')

            subplot.cla()

        plt.close()







