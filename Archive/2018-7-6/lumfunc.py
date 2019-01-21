from matplotlib import pyplot as plt
import numpy as np
import matplotlib as mpl
from pandas import read_csv
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
import deep2reader

font = {'family':'Roboto', 'weight':'light'}


c = 3e8 # in m/s

class lumfunc:

    def __init__(self, filename = './protoDC2_sed.dat'):

        (self.z,
            self.sed_1000_246, self.sed_11467_1710, self.sed_1246_306, self.sed_13177_1966, self.sed_15143_2259, 
            self.sed_1552_381, self.sed_17402_2596, self.sed_1933_474, self.sed_2407_591, self.sed_2998_186, 
            self.sed_3184_197, self.sed_3381_209, self.sed_3590_222, self.sed_3812_236, self.sed_4048_251, 
            self.sed_4299_266, self.sed_4565_283, self.sed_4848_300, self.sed_5148_319, self.sed_5467_339, 
            self.sed_5806_360, self.sed_6166_382, self.sed_6548_406, self.sed_6954_431, self.sed_7385_458, 
            self.sed_7843_486, self.sed_8329_517, self.sed_8846_549, self.sed_9395_583, self.sed_9978_1489,
            self.Ha_rest, self.Hb_rest, self.OII3726_rest, self.OII3729_rest, self.OIII4959_rest, self.OIII5007_rest) = np.transpose(read_csv(filename, header = None, comment = '#', delimiter = '\s+').values)

        self.Ha_sch_L_star = 10**42.56
        self.Ha_sau_L_star = None

        self.lam_Ha = 6563e-10     # All in meters
        self.lam_Hb = 4861e-10
        self.lam_OII = 3728.5e-10
        self.lam_OIII = 4975e-10

        self.nu_Ha = c/self.lam_Ha
        self.nu_Hb = c/self.lam_Hb
        self.nu_OII = c/self.lam_OII
        self.nu_OIII = c/self.lam_OIII

        # self.Ha_rest = self.Ha_rest*self.nu_Ha

        # self.OII3726_rest = self.OII3726_rest * self.nu_OII
        # self.OII3729_rest = self.OII3729_rest * self.nu_OII

        # self.OIII4959_rest = self.OIII4959_rest * self.nu_OIII
        # self.OIII5007_rest = self.OIII5007_rest * self.nu_OIII

        self.medz = round(np.median(self.z),1)

        l_lambda_Ha = self.sed_6548_406 * 3.0e10 * ((6548. + 203.) * 1.0e-8 )**-2.
        self.EW_Ha = 1.e8 * self.Ha_rest / l_lambda_Ha

        l_lambda_Hb = self.sed_4848_300 * 3.0e10 * ((4848. + 150.) * 1.0e-8 )**-2.
        self.EW_Hb = 1.e8 * self.Hb_rest / l_lambda_Hb

        l_lambda_OII = self.sed_3590_222 * 3.0e10 * ((3590. + 111.) * 1.0e-8 )**-2.
        self.EW_OII = 1.e8 * (self.OII3729_rest + self.OII3726_rest) / l_lambda_OII

        l_lambda_OIII = self.sed_4848_300 * 3.0e10 * ((4848. + 150.) * 1.0e-8)**-2.
        self.EW_OIII = 1.e8 * (self.OIII5007_rest) / l_lambda_OIII

        self.deep2 = deep2reader.d2r()



    def schechter(self, L, phi_star, L_star, alpha):

        norm_L = L/L_star

        return phi_star * norm_L**alpha * np.exp(-norm_L)


    def saunders(self, L, phi_star, L_star, alpha, sigma = 0.54):

        # Sigma from Comparat et al. (2016)

        norm_L = L/L_star

        exponent = -1 * (np.log10(1+norm_L) / (np.sqrt(2)*sigma))**2

        return phi_star * norm_L**alpha * np.exp(exponent)
 


    def Ha_lumfunc(self, L, z = 1.0, norm = True, model = 'schechter'):

        # z = 0.84, Sobral et al. (2012)

        if z != 1.0:

            print 'H-alpha luminosity function not configured for variable redshift.'

            return 0

        if model == 'schechter':

            phi_star = 10**-2.61 # Mpc^-3
            L_star = 10**42.56 # erg/s
            alpha = -1.62


            lumdist = self.schechter(L, phi_star, L_star, alpha)

        else:

            print 'Invalid model.'


        if norm:

            lumspace = 10**np.arange(np.log10(10**-2.*L_star), np.log10(10**2*L_star), 0.001)

            norm_const = np.trapz(self.Ha_lumfunc(lumspace, z = z, norm = False, model = model),lumspace)

        else:

            norm_const = 1.


        return lumdist/norm_const



    def OII_lumfunc(self, L, z = 1.0, norm = True, model = 'schechter'):

        # From Comparat et al. (2016)

        if model == 'schechter':

            phi_star = 10**-2.4 * (1.+z)**-.73
            L_star = 10**41.1 * (1.+z)**2.33 # erg/s
            alpha = -1.46

            lumdist = self.schechter(L, phi_star, L_star, alpha)

        elif model == 'saunders':

            phi_star = 10**-1.95 * (1.+z)**-.07
            L_star = 10**40.1 * (1.+z)**1.92 # erg/s
            alpha = -1.12

            lumdist = self.saunders(L, phi_star, L_star, alpha)

        else:

            print 'Invalid model.'


        if norm:

            lumspace = 10**np.arange(np.log10(10**-2.*L_star), np.log10(10**2*L_star), 0.001)

            norm_const = np.trapz(self.OII_lumfunc(lumspace, z = z, norm = False, model = model),lumspace)

        else:

            norm_const = 1.


        return lumdist/norm_const





    def OIII_lumfunc(self, L, z = 1.0, norm = True, model = 'schechter'):

        # From Comparat et al. (2016)

        if model == 'schechter':

            phi_star = 10**-3.41 * (1.+z)**-.76
            L_star = 10**41.42 * (1.+z)**3.91 # erg/s
            alpha = -1.83

            lumdist = self.schechter(L, phi_star, L_star, alpha)

        elif model == 'saunders':

            phi_star = 10**-2.91 * (1.+z)**-.22
            L_star = 10**40.81 * (1.+z)**3.31 # erg/s
            alpha = -1.81

            lumdist = self.saunders(L, phi_star, L_star, alpha)

        else:

            print 'Invalid model'


        if norm:

            lumspace = 10**np.arange(np.log10(10**-2.*L_star), np.log10(10**2*L_star), 0.001)

            norm_const = np.trapz(self.OIII_lumfunc(lumspace, z = z, norm = False, model = model),lumspace)

        else:

            norm_const = 1.


        return lumdist/norm_const




    def plot_lumfunc(self, model = 'schechter', lims = [33,44]):

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        lumspace = np.logspace(lims[0],lims[1],1000)

        # Plot Ha
        sp.plot(lumspace, self.Ha_lumfunc(lumspace, model = model, norm = True), color = 'r', linewidth = 2, label = r'H$\alpha$')

        # Plot OII
        sp.plot(lumspace, self.OII_lumfunc(lumspace, model = model, norm = True), color = 'g', linewidth = 2, label = '[OII]')

        # Plot OIII
        sp.plot(lumspace, self.OIII_lumfunc(lumspace, model = model, norm = True), color = 'b', linewidth = 2, label = '[OIII]')


        # Ha Histogram
        sp.hist(self.Ha_rest, bins = np.logspace(lims[0],lims[1],25), histtype = 'step', color = 'r', normed = True)

        # OII Histogram
        sp.hist(self.OII3729_rest, bins = np.logspace(lims[0],lims[1],25), histtype = 'step', color = 'g', normed = True)

        # OIII Histogram
        sp.hist(self.OIII5007_rest, bins = np.logspace(lims[0],lims[1],25), histtype = 'step', color = 'b', normed = True)

        sp.legend(loc = 'upper right', fontsize = 20)

        sp.set_xscale('log')
        sp.set_yscale('log')

        sp.set_ylim(10**-55, 10**-33)

        sp.set_xlabel(r'Luminosity ($L$)', family = 'Roboto', weight = 'light', fontsize = 24)
        sp.set_ylabel('Normalized Frequency', family = 'Roboto', weight = 'light', fontsize = 24)

        sp.set_title('Luminosity Functions for LSST Emission Lines', family = 'Roboto', weight = 'light', fontsize = 28)



    def plot_ew_dist(self):

        l_lambda_Ha = self.sed_6548_406 * 3.0e10 * ((6548. + 203.) * 1.0e-8 )**-2.
        EW_Ha = 1.e8 * self.Ha_rest / l_lambda_Ha

        l_lambda_OII = self.sed_3590_222 * 3.0e10 * ((3590. + 111.) * 1.0e-8 )**-2.
        EW_OII = 1.e8 * (self.OII3729_rest + self.OII3726_rest) / l_lambda_OII

        l_lambda_OIII = self.sed_4848_300 * 3.0e10 * ((4848. + 150.) * 1.0e-8)**-2.
        EW_OIII = 1.e8 * (self.OIII5007_rest) / l_lambda_OIII

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        sp.hist(EW_Ha, bins = np.logspace(-4, 4, 40), histtype = 'step', color = 'r', range = [0,10000], label = r'$H\alpha$')
        sp.hist(EW_OII, bins = np.logspace(-4, 4, 40), histtype = 'step', color = 'g', range = [0,10000], label = r'[OII]$_{3726 + 3729}$')
        sp.hist(EW_OIII, bins = np.logspace(-4, 4, 40), histtype = 'step', color = 'b', range = [0,10000], label = r'[OIII]$_{5007}$')

        sp.set_xscale('log')
        sp.set_xlim(0,10000)

        sp.legend(loc = 'upper left', fontsize = 16)

        sp.set_xlabel(r'Equivalent Width ($\mathrm{\AA}$)', family = 'Roboto', weight = 'light', fontsize = 24)
        sp.set_ylabel('Frequency', family = 'Roboto', weight = 'light', fontsize = 24)

        sp.set_title(r'Equivalent Width Distribution for $z \approx ' + '%.1f$' % np.median(self.z), family = 'Roboto', weight = 'light', fontsize = 28)




    def plot_ha_ew(self):

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        bins = np.logspace(-4, 4, 40)

        sp.hist(self.EW_Ha, bins = bins, histtype = 'step', color = 'r', label = r'$H\alpha$')

        sp.set_xscale('log')
        sp.set_xlim(min(bins),max(bins))

        sp.set_xlabel(r'Equivalent Width ($\mathrm{\AA}$)', family = 'Roboto', weight = 'light', fontsize = 24)
        sp.set_ylabel('Frequency', family = 'Roboto', weight = 'light', fontsize = 24)

        sp.set_title(r'H$\alpha$ Equivalent Width Distribution for $z \approx ' + '%.1f$' % np.median(self.z), family = 'Roboto', weight = 'light', fontsize = 28)




    def plot_ew_frac(self):

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        # bins = np.logspace(-2,2, 40)
        bins = np.linspace(-2, 2, 40)

        sp.hist(np.log10(self.EW_OIII/self.EW_Hb), bins = bins, histtype = 'step', linewidth = 2, color = 'b', normed = True, label = r'[OIII]/H$\beta$')
        sp.hist(np.log10(self.EW_OIII/self.EW_OII), bins = bins, histtype = 'step', linewidth = 2, color = 'g', normed = True, label = r'[OIII]/[OII]')
        sp.hist(np.log10(self.EW_OII/self.EW_Ha), bins = bins, histtype = 'step', linewidth = 2, color = 'r', normed = True, label = r'[OII]/H$\alpha$')

        sp.hist(np.log10(self.deep2.OIIIHb[self.deep2.OIIIHb > 0]), bins = bins, histtype = 'stepfilled', color = 'b', alpha = 0.4, normed = True, label = r'DEEP2 [OIII]/H$\beta$')
        sp.hist(np.log10(self.deep2.OIIIOII[self.deep2.OIIIOII > 0]), bins = bins, histtype = 'stepfilled', color = 'g', alpha = 0.4, normed = True, label = r'DEEP2 [OIII]/[OII]')
        
        # sp.set_xscale('log')
        sp.set_xlim(min(bins),max(bins))

        sp.legend(loc = 'upper left', fontsize = 16)

        sp.set_xlabel(r'Equivalent Width Ratio', family = 'Roboto', weight = 'light', fontsize = 24)
        sp.set_ylabel('Frequency', family = 'Roboto', weight = 'light', fontsize = 24)

        sp.set_title(r'Equivalent Width Ratio Distribution for $z \approx ' + '%.1f$' % np.median(self.z), family = 'Roboto', weight = 'light', fontsize = 28)





    def plot_lum_frac(self):

        l_lambda_Ha = self.sed_6548_406 * 3.0e10 * ((6548. + 203.) * 1.0e-8 )**-2.
        EW_Ha = 1.e8 * self.Ha_rest / l_lambda_Ha

        l_lambda_OII = self.sed_3590_222 * 3.0e10 * ((3590. + 111.) * 1.0e-8 )**-2.
        EW_OII = 1.e8 * (self.OII3729_rest + self.OII3726_rest) / l_lambda_OII

        l_lambda_OIII = self.sed_4848_300 * 3.0e10 * ((4848. + 150.) * 1.0e-8)**-2.
        EW_OIII = 1.e8 * (self.OIII5007_rest + self.OIII4959_rest) / l_lambda_OIII

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        sp.hist(EW_OIII/EW_OII, bins = np.logspace(-1, 2, 40), histtype = 'step', color = 'teal', range = [0,10000], label = r'[OIII]$_{4959 + 5007}$/[OII]$_{3726 + 3729}$')
        sp.hist(EW_Ha/EW_OII, bins = np.logspace(-1, 2, 40), histtype = 'step', color = 'brown', range = [0,10000], label = r'$H\alpha$/[OII]$_{3726 + 3729}$')
        sp.hist(EW_Ha/EW_OIII, bins = np.logspace(-1, 2, 40), histtype = 'step', color = 'purple', range = [0,10000], label = r'$H\alpha$/[OIII]$_{4959 + 5007}$')

        sp.set_xscale('log')
        sp.set_xlim(10**-1,10**2)
        sp.set_ylim(0,3250)

        sp.legend(loc = 'upper left', fontsize = 14)

        sp.set_xlabel(r'Equivalent Width Ratio', family = 'Roboto', weight = 'light', fontsize = 24)
        sp.set_ylabel('Frequency', family = 'Roboto', weight = 'light', fontsize = 24)

        sp.set_title(r'Equivalent Width Ratios for $z \approx' + ' %.1f$' % np.median(self.z), family = 'Roboto', weight = 'light', fontsize = 28)
