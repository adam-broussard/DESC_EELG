import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from glob import glob
import os
import filtersim
import photosim
from pandas import read_csv
from tqdm import tqdm



font = {'family':'Roboto', 'weight':'light'}
phot = photosim.photosim()


def eazy_store(input_fp = './EAZY_runs/', output_fp = './Outputs/'):

    files = sorted(glob('./EAZY_runs/*'))

    # x = 1

    # while os.path.isdir(output_fp + 'run%03i/' % x):
    #     x += 1

    # outputdir = output_fp + 'run%03i/' % x

    if not os.path.isdir(output_fp + '_EAZY/'):
        os.makedirs(output_fp + '_EAZY/')
    else:
        os.system('rm ' + output_fp + '_EAZY/*')

    os.system('cp -r ' + input_fp + 'Output/* ' + output_fp + '_EAZY/')
    os.system('cp ' + input_fp + 'cat.dat ' + output_fp + '_EAZY/')
    os.system('cp ' + input_fp + 'cat_input.param ' + output_fp + '_EAZY/')

    return output_fp



def diagnostic_plots(outputdir = './Outputs/', diagnostic_text = True):
    
    if outputdir[-1] != '/':
        outputdir = outputdir + '/'

    eazy_store(output_fp = outputdir)

    if not os.path.isdir(outputdir + 'Figures/'):
        os.mkdir(outputdir + 'Figures/')
    else:
        os.system('rm ' + outputdir + 'Figures/*')

    filters = filtersim.filtersim()

    temp_sed_list = sorted(glob(outputdir + '_EAZY/*.temp_sed'))
    obs_sed_list = sorted(glob(outputdir + '_EAZY/*.obs_sed'))

    z_fit = np.loadtxt(outputdir + '_EAZY/photz.zout', unpack = True, usecols = [2])

    ID, z_orig, age, sfh, metallicity, imf = read_csv(outputdir + '_EAZY/cat_input.param',
            header = None, comment = '#', delimiter = '\s+').values.T 
    sfh_type = np.zeros(len(sfh)) # SSP
    sfh_type[sfh == 'CSF'] = 3

    imf_type = np.zeros(len(imf)) # SALPETER
    imf_type[imf == 'CHABRIER'] = 1
    imf_type[imf == 'KROUPA'] = 2

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

    for x in tqdm(xrange(len(obs_waves))):

        index = np.where(ID == int(temp_sed_list[x].split('_')[-2].split('.')[0]))[0]

        if not len(index) == 0:

            index = index[0]

            fig = plt.figure(figsize = (8,8))
            sp = fig.add_subplot(111)

            sp.plot(temp_waves[index], temp_fluxes[index], linewidth = 2, color = 'k', label = 'Fit')
            sp.scatter(obs_waves[index], cat_fluxes[index], s = 100, facecolor = 'darkred', edgecolors = 'None')
            sp.scatter(obs_waves[index], obs_fluxes[index], s = 100, facecolor = 'None', edgecolors = 'k', linewidth = 2)

            # NEEDS UPDATING FROM HERE ON

            sp.errorbar(obs_waves[index], cat_fluxes[index], yerr = fit_err, color = 'darkred', fmt = ' ', zorder = 1)
            orig_wave, orig_flux = phot.redshift(z_orig[index], *phot.find_spectrum(age[index], metallicity[index], imf_type[index], sfh_type[index]))
            sp.plot(orig_wave, orig_flux, linewidth = 2, color = 'r', label = 'Simulated')


            sp.set_xlim(3000, 12000)
            # sp.set_ylim(10**-31, 10**-24)
            sp.set_ylim(-0.25*10**-27, 0.5*10**-27)
            # sp.set_yscale('log')


            # sp.autoscale(enable = True, axis = 'y', tight = None)

            sp.text(0.98, 0.08, r'z$_\mathrm{orig}$ = ' + '%.2f' % z_orig[index], ha = 'right', va = 'bottom', fontdict = font, fontsize = 36, transform = sp.transAxes)
            sp.text(0.98, 0.02, r'z$_\mathrm{fit}$ = ' + '%.2f' % z_fit[index], ha = 'right', va = 'bottom', fontdict = font, fontsize = 36, transform = sp.transAxes)
            # sp.text(0.02, 0.02, '%03i' % ID[index], ha = 'left', va = 'bottom', fontdict = font, fontsize = 36, transform = sp.transAxes)
            
            sp.set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontdict = font, fontsize = 24)
            sp.set_ylabel(r'Flux Density ($\mathrm{erg\,s^{-1}\,cm^{-2}\,Hz^{-1}}$)', fontdict = font, fontsize = 24)

            sp.legend(loc = 'upper right', fontsize = 18)

            if diagnostic_text:

                sp.text(0.02, 0.02, 'ID: %03i\nAge: %.3fGyr\nZ/Z$_\odot$: %.3f\nSFH: %i\nIMF: %i' % (ID[index], age[index], metallicity[index], sfh_type[index], imf_type[index]),
                    fontdict = font, fontsize = 14, ha = 'left', va = 'bottom', transform = sp.transAxes)

            plt.savefig(outputdir + 'Figures/%i.png' % ID[index], bbox_inches = 'tight')
            plt.close()



def plot_templates(specf = None, outputdir = './Outputs/'):

    if specf == None:
        keys, vals = np.loadtxt('./EAZY_runs/zphot.param', dtype = str, unpack = True)
        specf = vals[keys == 'TEMPLATES_FILE'][0]
        specf = './EAZY_runs/' + specf

    fnames = list(np.loadtxt(specf, delimiter = '   ', usecols = [1], dtype = str))

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir + 'Figures/')
    else:
        os.system('rm ' + outputdir + 'Figures/templates.png')

    fig = plt.figure(figsize = (8,8))
    sp = fig.add_subplot(111)

    for (x, fname) in enumerate(fnames):

        lam, f_lam = np.loadtxt('./EAZY_runs/' + fname, unpack = True)

        # Find the closest wavelength to 6500 AA

        closest = np.argmin(abs(lam - 10000.))

        # sp.plot(lam, (f_lam / f_lam[closest]) + x, linewidth = 2, label = fname.split('/')[-1][:-4])
        sp.plot(lam, (f_lam / max(f_lam)) + (len(fnames)-x-1), linewidth = 2, label = fname.split('/')[-1][:-4])

    sp.set_xlabel(r'$\lambda$ ($\AA$)', fontdict = font, fontsize = 24)
    sp.set_ylabel(r'F$_\lambda$', fontdict = font, fontsize = 24)
    sp.legend(loc = 'upper right', fontsize = 8)
    sp.set_ylim(-0.5, 10.5)

    sp.set_xscale('log')

    plt.savefig(outputdir + 'Figures/templates.png', bbox_inches = 'tight')
    plt.close()



def runeazy(dir = './EAZY_runs/'):

    current_directory = os.getcwd()
    os.chdir(dir)
    os.system('rm ./Output/*')
    os.system('/home/adam/Research/eazy-photoz/src/eazy > ./Output/log.txt')
    os.chdir(current_directory)


def plot_sigma_deviation(dir = './Outputs/_EAZY/'):

    obs_sed_list = sorted(glob(dir + '*.obs_sed'))

    ID, z_orig, age, sfh, metallicity, imf = read_csv(dir + 'cat_input.param',
        header = None, comment = '#', delimiter = '\s+').values.T 

    agelist = []
    sigmalist = []

    for (x, obsfile) in enumerate(obs_sed_list):

        fileid = int(obsfile[-11:-8])
        agelist.append(age[ID == fileid])

        obs_flux, fit_err, temp_flux = np.loadtxt(obsfile, usecols = [1,3,4], unpack = True)

        sigmalist.append(abs((obs_flux - temp_flux) / fit_err)) #-2 is for z band
        # Try taking out the -2 here to see if it works with all bands!

    agelist = np.array(agelist)
    sigmalist = np.array(sigmalist)

    agelist = agelist + 1.
    # print agelist
    # print sigmalist

    # agelist, sigmalist = (list(t) for t in zip(*sorted(zip(agelist, sigmalist))))


    print 'Warning: Modifying CSF starting age...'

    fig = plt.figure(figsize = (8,8))
    sp = fig.add_subplot(111)

    for (x,sigma) in enumerate(sigmalist.T):
        sp.plot(agelist, sigma, linewidth = 2, label = ['u', 'g', 'r', 'i', 'z', 'y'][x])
    
    sp.set_xlabel('Time (Gyr)', fontsize = 24, fontdict = font)
    sp.set_ylabel(r'$\sigma$ Deviations', fontsize = 24, fontdict = font)
    sp.legend(loc = 'upper right', fontsize = 18)




