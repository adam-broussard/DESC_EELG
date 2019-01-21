import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from glob import glob
import os
import filtersim
import photosim



font = {'family':'Roboto', 'weight':'light'}
phot = photosim.photosim()


def diagnostic_plots(root_fp = './EAZY_runs/Inputs/', input_fp = 'Output/', output_fp = 'diagplots/', default_cat = True, logrange = [-30,-25.5]):



    if default_cat:
        metallicities = [0, ]*10 + [1,]*10
        time = list(range(10)) + list(range(10))


    filters = filtersim.filtersim()

    if not os.path.isdir(root_fp + output_fp):
        os.makedirs(root_fp + output_fp)

    temp_sed_list = sorted(glob(root_fp + input_fp + '*.temp_sed'))
    obs_sed_list = sorted(glob(root_fp + input_fp + '*.obs_sed'))

    ID, z = np.loadtxt(root_fp + input_fp + 'photz.zout', unpack = True, usecols = [0,2])

    temp_waves = []
    temp_fluxes = []
    obs_waves = []
    obs_fluxes = []
    cat_fluxes = []

    for x in range(len(temp_sed_list)):

        thiswave, thistempflux = np.loadtxt(temp_sed_list[x], unpack = True, usecols = [0,1])

        temp_waves.append(thiswave)
        temp_fluxes.append(thistempflux)

        thiswave, thiscatflux, fit_err, thisobsflux = np.loadtxt(obs_sed_list[x], unpack = True, usecols = [0,1,3,4])

        obs_waves.append(thiswave)
        obs_fluxes.append(thisobsflux)
        cat_fluxes.append(thiscatflux)

    for x in range(len(obs_waves)):

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        sp.plot(temp_waves[x], temp_fluxes[x], linewidth = 2, color = 'k', label = 'Fit')
        sp.scatter(obs_waves[x], cat_fluxes[x], s = 100, facecolor = 'darkred', edgecolors = 'None')
        sp.scatter(obs_waves[x], obs_fluxes[x], s = 100, facecolor = 'None', edgecolors = 'k', linewidth = 2)

        if default_cat:

            orig_wave, orig_flux = phot.redshift(0.3, phot.wavelength, phot.l_nu[metallicities[x]][time[x]])
            sp.plot(orig_wave, orig_flux, linewidth = 2, color = 'r', label = 'Simulated')
            sp.errorbar(obs_waves[x], cat_fluxes[x], yerr = fit_err, color = 'darkred', fmt = ' ', zorder = 1)


        index = np.where(ID == int(temp_sed_list[x].split('_')[-2].split('.')[0]))[0]

        sp.set_xlim(3000, 12000)
        sp.set_ylim(10**logrange[0], 10**logrange[1])

        sp.set_yscale('log')

        sp.text(0.98, 0.08, r'z$_\mathrm{orig}$ = ' + '%.2f' % 0.3, ha = 'right', va = 'bottom', fontdict = font, fontsize = 36, transform = sp.transAxes)
        sp.text(0.98, 0.02, r'z$_\mathrm{fit}$ = ' + '%.2f' % z[index], ha = 'right', va = 'bottom', fontdict = font, fontsize = 36, transform = sp.transAxes)
        # sp.text(0.02, 0.02, '%03i' % ID[index], ha = 'left', va = 'bottom', fontdict = font, fontsize = 36, transform = sp.transAxes)
        
        sp.set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontdict = font, fontsize = 24)
        sp.set_ylabel(r'Flux Density ($\mathrm{erg\,s^{-1}\,cm^{-2}\,Hz^{-1}}$)', fontdict = font, fontsize = 24)

        sp.legend(loc = 'upper right', fontsize = 18)

        plt.savefig(root_fp + output_fp + '%i.png' % ID[index], bbox_inches = 'tight')
        plt.close()



def plot_templates(specf = '../eazy-photoz/templates/eazy_v1.2_dusty_modified.spectra.param'):

    fnames = list(np.loadtxt(specf, delimiter = '   ', usecols = 2, dtype = str))

    for (x, fname) in enumerate(fnames):

        pass
