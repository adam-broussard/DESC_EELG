import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from glob import glob

noNewLine = '\x1b[1A\x1b[1M'
font = {'family':'Roboto', 'weight':'light'}

class filtersim:

    def __init__(self, filter_fp = './lsst_filts/*.dat', manualkeys = ['u', 'g', 'r', 'i', 'z', 'y']):

        filterfiles = glob(filter_fp)

        self.wavelength = {}
        self.response = {}
        self.norm_response = {}

        for fname in filterfiles:

            tempwave, tempresp = np.loadtxt(fname, unpack = True)
            self.wavelength[fname.split('/')[-1][-5]] = tempwave*10.
            self.response[fname.split('/')[-1][-5]] = tempresp

            self.norm_response[fname.split('/')[-1][-5]] = tempresp/np.trapz(tempresp, x = tempwave*10.)

        if manualkeys:
            self.keys = manualkeys
        else:
            self.keys = self.wavelength.keys()

        self.wavelength_centers = {}

        for index, key in enumerate(self.keys):
            
            self.wavelength_centers[key] = np.trapz(self.norm_response[key] * self.wavelength[key], x = self.wavelength[key])

        # From the LSST Sience book pg 21 https://arxiv.org/pdf/0912.0201.pdf
        # Units are erg/s/cm^2/Hz

        self.error = {'u':np.power(10., (26.3 + 48.6)/(-2.5)),
                        'g':np.power(10., (27.5 + 48.6)/(-2.5)),
                        'r':np.power(10., (27.7 + 48.6)/(-2.5)),
                        'i':np.power(10., (27.0 + 48.6)/(-2.5)),
                        'z':np.power(10., (26.2 + 48.6)/(-2.5)),
                        'y':np.power(10., (24.9 + 48.6)/(-2.5))}




    def plot_filters(self):

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        colors = [mpl.cm.jet(int(place)) for place in np.linspace(0, mpl.cm.jet.N, len(self.keys))]

        for index, key in enumerate(self.keys):

            sp.plot(self.wavelength[key], self.response[key], color = colors[index], linewidth = 2, label = key)

        # sp.set_xscale('log')
        sp.set_ylim(-0.05, 1.0)
        sp.legend(loc = 'upper right', fontsize = 16)

        sp.set_xlabel(r'$\lambda\,\,(\mathrm{\AA})$', fontdict = font, fontsize = 24)
        sp.set_ylabel('Throughput', fontdict = font, fontsize = 24)
        sp.set_title('LSST Throughput Curves', fontdict = font, fontsize = 28)




    def get_photometry_old(self, wavelength, flux_density, wavestep = 0.1):

        wavelengths = [self.wavelength_centers[key] for key in self.keys]
        phot = []

        interp_x = np.arange(3000, 12000, wavestep)

        for index, key in enumerate(self.keys):

            interp_y = np.interp(interp_x, wavelength, flux_density)

            filter_interp_y = np.interp(interp_x, self.wavelength[key], self.norm_response[key])

            phot.append(np.trapz(filter_interp_y * interp_y, x = interp_x))

            phot_err = [self.error[key] for key in self.keys]

        return wavelengths, phot, phot_err



    def get_photometry(self, wavelength, flux_density, wavestep = 0.1):

        # This is slightly different from the above because LSST filters are photon count filters

        c = 3.e10 # speed of light in cm/s

        f_lambda = flux_density * c/(wavelength**2.)

        wavelengths = [self.wavelength_centers[key] for key in self.keys]
        phot = []

        interp_x = np.arange(3000, 12000, wavestep)

        for index, key in enumerate(self.keys):

            interp_y = np.interp(interp_x, wavelength, f_lambda)

            filter_interp_y = np.interp(interp_x, self.wavelength[key], self.norm_response[key])

            phot.append(np.trapz(filter_interp_y * interp_y * interp_x, x = interp_x)/np.trapz(filter_interp_y * c / interp_x, x = interp_x))

            phot_err = [self.error[key] for key in self.keys]

        return wavelengths, phot, phot_err



    def output_eazy_filters(self, fname = './EAZY_runs/Inputs/lsst.filters.res'):

        with open(fname, 'w') as writefile:

            for index, key in enumerate(self.keys):

                writefile.write(str(len(self.wavelength[key])).ljust(10) + 'LSST ' + key 
                    + '-band Filter; lam_ctr = %.2f\n' % self.wavelength_centers[key])

                for x in range(len(self.wavelength[key])):

                    writefile.write('%6i' % x)
                    writefile.write('  %.2f'.ljust(12) % self.wavelength[key][x])
                    writefile.write('%.4f\n' % self.response[key][x])


    def output_filter_translate(self, fname = './EAZY_runs/Inputs/zphot.translate'):

        with open(fname, 'w') as writefile:

            for index, key in enumerate(self.keys):

                writefile.write('f_' + key + '  F%i\n' % (index+1))
                writefile.write('e_' + key + '  E%i\n' % (index+1))


















