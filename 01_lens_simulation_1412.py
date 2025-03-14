
#parent: lsst_lens_281124

# -*- coding: utf-8 -*-
"""
 This script creates mock lenses with parameters in their given ranges.

Original file is located at
    https://colab.research.google.com/github/lenstronomy/lenstronomy-tutorials/blob/main/Notebooks/LensModeling/modeling_a_simple_Einstein_ring.ipynb

# Simple ring parameter constraints

This notebook provides test cases for the precision of a simple lens model (with simplified assumptions). This is for show-casing and to assess the uncertainty limit in how well the parameters of this model can be constrained.
"""

import numpy as np
import os
import time
import corner
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import datetime
import lenstronomy
from astropy.io import fits
from streaks import create_streak

c = datetime.now()
# Displays Time
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
print("current time :", current_time)

if not os.path.exists("/data/a.saricaoglu/lenstronomy/lenstronomy-data/mock_lenses/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M')): 
    os.makedirs("/data/a.saricaoglu/lenstronomy/lenstronomy-data/mock_lenses/"+  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M')) 

#In this part, we create and write on the .fits files.

#PrimaryHDU will have the image data, with image related headers. Data is empty for now.
hdu_pr = fits.PrimaryHDU()
#ImageHDU is where we keep the mock lens data, with lens related headers.
# Write to a FITS file
hdu_pr.writeto("/data/a.saricaoglu/lenstronomy/lenstronomy-data/mock_lenses/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M') +"_simulation.fits", overwrite=True)
# Verify the file by reading it back
hdulist = fits.open("/data/a.saricaoglu/lenstronomy/lenstronomy-data/mock_lenses/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M') +"_simulation.fits", mode='update')
print(hdulist[0].header)

noisy_min = []
noisy_max = []
image_min = []
image_max = []

lens_number = 100
"""## simulation choices"""

for i in range(0,lens_number):
    # data specifics
    background_rms = .005  #  background noise per pixel
    exp_time = 30  #  exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
    numPix = 100  #  cutout pixel size per axis
    pixel_scale = 0.2  #  pixel size in arcsec (area per pixel = pixel_scale**2)
    fwhm = 0.7  # full width at half maximum of PSF
    psf_type = 'GAUSSIAN'  # 'GAUSSIAN', 'PIXEL', 'NONE'

    # lensing quantities
    lens_model_list = ['SIE', 'SHEAR']
    kwargs_spemd = {'theta_E': np.random.uniform(4, 12.0), 'center_x': 0.05, 'center_y': 0, 'e1': np.random.uniform(-0.2, 0.2), 'e2': np.random.uniform(-0.2, 0.2)}  # parameters of the deflector lens model
    kwargs_shear = {'gamma1': 0.0, 'gamma2': -0.05}  # shear values to the source plane

    kwargs_lens = [kwargs_spemd, kwargs_shear]
    from lenstronomy.LensModel.lens_model import LensModel
    lens_model_class = LensModel(lens_model_list)


    # Sersic parameters in the initial simulation for the source
    kwargs_sersic = {'amp': 16, 'R_sersic': 0.1, 'n_sersic': 1, 'e1': np.random.uniform(-0.2, 0.2), 'e2': np.random.uniform(-0.2, 0.2), 'center_x': 0.1, 'center_y': 0}
    source_model_list = ['SERSIC_ELLIPSE']
    kwargs_source = [kwargs_sersic]


    from lenstronomy.LightModel.light_model import LightModel
    source_model_class = LightModel(source_model_list)


    kwargs_sersic_lens = {'amp': 16, 'R_sersic': 0.6, 'n_sersic': 2, 'e1': -0.1, 'e2': 0.1, 'center_x': 0.05, 'center_y': 0}

    lens_light_model_list = ['SERSIC_ELLIPSE']
    kwargs_lens_light = [kwargs_sersic_lens]
    lens_light_model_class = LightModel(lens_light_model_list)

    # import main simulation class of lenstronomy
    from lenstronomy.Util import util
    from lenstronomy.Data.imaging_data import ImageData
    from lenstronomy.Data.psf import PSF
    import lenstronomy.Util.image_util as image_util
    from lenstronomy.ImSim.image_model import ImageModel

    # generate the coordinate grid and image properties (we only read out the relevant lines we need)
    _, _, ra_at_xy_0, dec_at_xy_0, _, _, Mpix2coord, _ = util.make_grid_with_coordtransform(numPix=numPix, deltapix=pixel_scale, center_ra=0, center_dec=0, subgrid_res=1, inverse=False)


    kwargs_data = {'background_rms': background_rms,  # rms of background noise
                'exposure_time': exp_time,  # exposure time (or a map per pixel)
                'ra_at_xy_0': ra_at_xy_0,  # RA at (0,0) pixel
                'dec_at_xy_0': dec_at_xy_0,  # DEC at (0,0) pixel
                'transform_pix2angle': Mpix2coord,  # matrix to translate shift in pixel in shift in relative RA/DEC (2x2 matrix). Make sure it's units are arcseconds or the angular units you want to model.
                'image_data': np.zeros((numPix, numPix))  # 2d data vector, here initialized with zeros as place holders that get's overwritten once a simulated image with noise is created.
                }

    data_class = ImageData(**kwargs_data)
    # generate the psf variables
    kwargs_psf = {'psf_type': 'GAUSSIAN', 'fwhm': fwhm, 'pixel_size': pixel_scale, 'truncation': 3}

    # if you are using a PSF estimate from e.g. a star in the FoV of your exposure, you can set
    #kwargs_psf = {'psf_type': 'PIXEL', 'pixel_size': deltaPix, 'kernel_point_source': 'odd numbered 2d grid with centered star/PSF model'}


    psf_class = PSF(**kwargs_psf)
    kwargs_numerics = {'supersampling_factor': 1, 'supersampling_convolution': False}

    imageModel = ImageModel(data_class, psf_class, lens_model_class=lens_model_class,
                            source_model_class=source_model_class, lens_light_model_class=lens_light_model_class,
                            kwargs_numerics=kwargs_numerics)

    # generate image
    image_model = imageModel.image(kwargs_lens, kwargs_source, kwargs_lens_light=kwargs_lens_light, kwargs_ps=None)

    poisson = image_util.add_poisson(image_model, exp_time=exp_time)
    bkg = image_util.add_background(image_model, sigma_bkd=background_rms)
    #streak = image_util.add_layer2image(image_model, 50,50, create_streak(100,100,1000))
    image_noisy = image_model + bkg + poisson
    #image_streaky = image_noisy + streak


    # display the initial simulated image
    import matplotlib as mpl
    cmap = mpl.cm.get_cmap("viridis").copy()
    cmap.set_bad(color='k', alpha=1.)
    cmap.set_under('k')

    v_min = -4
    v_max = 1

    f, axex = plt.subplots(1, 2,  figsize=(10, 6), sharex=False, sharey=False, gridspec_kw={'width_ratios': [1, 1.24]})


    noisy_min.append(np.min(image_noisy))
    noisy_max.append(np.max(image_noisy))
    image_min.append(np.min(image_model))
    image_max.append(np.max(image_model))
    # f, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=False, sharey=False)
    axex[0].matshow((image_model), origin='lower', cmap=cmap, vmin=np.percentile(image_model, 5), vmax=np.percentile(image_model, 95))
    #f.colorbar(axex[0].images[0],  ax=axex[0], label='Intensity', shrink=0.9)
    axex[1].matshow((image_noisy), origin='lower', cmap=cmap,vmin=np.percentile(image_noisy, 5), vmax=np.percentile(image_noisy, 95))
    #f.colorbar(axex[1].images[0], ax=axex[1], label='Intensity', shrink=0.9)
    #axex[2].matshow((image_streaky), origin='lower', cmap=cmap,vmin=np.percentile(image_noisy, 5), vmax=np.percentile(image_noisy, 95))
    f.colorbar(axex[1].images[0],  ax=axex[1], shrink=0.7)
    axex[0].set_xlabel("x")
    axex[0].set_ylabel("y")
    axex[1].set_xlabel("x")
    #axex[2].set_xlabel("x")

    # title1 = {k: round(v, 2) for k, v in kwargs_light_source[0].items()}
    # title2 = {k: round(v, 2) for k, v in kwargs_light_lens[0].items()}
    # axex[0].set_title("z-src :" + "%.2f" % z_source + "   Src :"+"\n z-model :" + "%.2f" % redshift_list[0] + "   Lens :"  )
    # axex[1].set_title( str(title1) + "\n" + str(title2))
    f.savefig("/data/a.saricaoglu/lenstronomy/lenstronomy-data/mock_lenses/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M') +"/lens_models_images_" + str(i) + "_" + current_time + ".png",  bbox_inches='tight')
    plt.close()

    image_hdu1 = fits.ImageHDU(image_model)
    image_hdu1.header['LENS_ID'] = str(i)
    image_hdu1.header['TYPE'] = 'Image'

    image_hdu2 = fits.ImageHDU(image_noisy)   
    image_hdu2.header['LENS_ID'] = str(i)
    image_hdu2.header['TYPE'] = 'Image with background and poisson noise'

    # image_hdu3 = fits.ImageHDU(image_streaky)   
    # image_hdu3.header['LENS_ID'] = str(i)
    # image_hdu3.header['TYPE'] = 'Image with background and poisson noise and streak'
    # Write the table to a new FITS file
    hdulist.append(image_hdu2)
    #hdulist.append(image_hdu3)

    # # """## Model fitting
    # # in the blocks above we simulated a mock lens with noise properties. From these cells, we only require the kwargs_data and kwargs_psf arguments to perform the modeling. If you have real data, you can leave out the image simulation and directly read in the data, PSF and noise properties into the keyword arguement list. Make sure the units are correct. Further information on the settings are available in the ImageData() and PSF() classes in the lenstronomy.Data module.
    # # """

    # # # lens models
    # for im in ['noisy_', 'streaky_']:
    #     if im == 'noisy_':
    #         # update the image in data class with the without streak added image
    #         data_class.update_data(image_noisy)
    #         kwargs_data['image_data'] = image_noisy
    #     else:
    #         # update the image in data class with the streak added image
    #         data_class.update_data(image_streaky)
    #         kwargs_data['image_data'] = image_streaky

    #     fixed_lens = []
    #     kwargs_lens_init = []
    #     kwargs_lens_sigma = []
    #     kwargs_lower_lens = []
    #     kwargs_upper_lens = []

    #     fixed_lens.append({})  # for this example, we fix the power-law index of the lens model to be isothermal
    #     kwargs_lens_init.append({'theta_E': 0.7, 'e1': 0., 'e2': 0.,
    #                             'center_x': 0., 'center_y': 0.})
    #     kwargs_lens_sigma.append({'theta_E': .2, 'e1': 0.05, 'e2': 0.05,
    #                             'center_x': 0.05, 'center_y': 0.05})
    #     kwargs_lower_lens.append({'theta_E': 3, 'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10}) #original params: 'theta_E': 0.01, 'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10
    #     kwargs_upper_lens.append({'theta_E': 13., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10}) #original params: theta_E': 10., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10

    #     fixed_lens.append({'ra_0': 0, 'dec_0': 0})
    #     kwargs_lens_init.append({'gamma1': 0., 'gamma2': 0.0})
    #     kwargs_lens_sigma.append({'gamma1': 0.1, 'gamma2': 0.1})
    #     kwargs_lower_lens.append({'gamma1': -0.2, 'gamma2': -0.2})
    #     kwargs_upper_lens.append({'gamma1': 0.2, 'gamma2': 0.2})

    #     lens_params = [kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens]


    #     fixed_source = []
    #     kwargs_source_init = []
    #     kwargs_source_sigma = []
    #     kwargs_lower_source = []
    #     kwargs_upper_source = []


    #     fixed_source.append({})
    #     kwargs_source_init.append({'R_sersic': 0.2, 'n_sersic': 1, 'e1': 0, 'e2': 0, 'center_x': 0., 'center_y': 0, 'amp': 16})
    #     kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.1, 'e1': 0.05, 'e2': 0.05, 'center_x': 0.2, 'center_y': 0.2, 'amp': 10})
    #     kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.001, 'n_sersic': .5, 'center_x': -10, 'center_y': -10, 'amp': 0})
    #     kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10, 'n_sersic': 5., 'center_x': 10, 'center_y': 10, 'amp': 100})

    #     source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]


    #     fixed_lens_light = []
    #     kwargs_lens_light_init = []
    #     kwargs_lens_light_sigma = []Automatically generated by Colab.

    #     fixed_lens_light.append({})
    #     kwargs_lens_light_init.append({'R_sersic': 0.5, 'n_sersic': 2, 'e1': 0, 'e2': 0, 'center_x': 0., 'center_y': 0, 'amp': 16})
    #     kwargs_lens_light_sigma.append({'n_sersic': 1, 'R_sersic': 0.3, 'e1': 0.05, 'e2': 0.05, 'center_x': 0.1, 'center_y': 0.1, 'amp': 10})
    #     kwargs_lower_lens_light.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.001, 'n_sersic': .5, 'center_x': -10, 'center_y': -10, 'amp': 0})
    #     kwargs_upper_lens_light.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10, 'n_sersic': 5., 'center_x': 10, 'center_y': 10, 'amp': 100})

    #     lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma, fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]

    #     kwargs_params = {'lens_model': lens_params,
    #                     'source_model': source_params,
    #                     'lens_light_model': lens_light_params}

    #     kwargs_likelihood = {'source_marg': False}
    #     kwargs_model = {'lens_model_list': lens_model_list, 'source_light_model_list': source_model_list, 'lens_light_model_list': lens_light_model_list}

    #     multi_band_list = [[kwargs_data, kwargs_psf, kwargs_numerics]]
    #     # if you have multiple  bands to be modeled simultaneously, you can append them to the mutli_band_list
    #     kwargs_data_joint = {'multi_band_list': multi_band_list,
    #                         'multi_band_type': 'single-band'  # 'multi-linear': every imaging band has independent solutions of the surface brightness, 'joint-linear': there is one joint solution of the linear coefficients demanded across the bands.
    #                         }
    #     kwargs_constraints = {'linear_solver': True}  # optional, if 'linear_solver': False, lenstronomy does not apply a linear inversion of the 'amp' parameters during fitting but instead samples them.

    #     from lenstronomy.Workflow.fitting_sequence import FittingSequence
    #     fitting_seq = FittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)

    #     fitting_kwargs_list = [['PSO', {'sigma_scale': 1., 'n_particles': 200, 'n_iterations': 200}],
    #                         ['MCMC', {'n_burn': 200, 'n_run': 600, 'n_walkers': 200, 'sigma_scale': .1}]
    #             ]

    #     chain_list = fitting_seq.fit_sequence(fitting_kwargs_list)
    #     kwargs_result = fitting_seq.best_fit()

    #     """## analyse MCMC chain"""

    #     from lenstronomy.Plots import chain_plot
    #     from lenstronomy.Plots.model_plot import ModelPlot

    #     modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat",
    #                         linear_solver=kwargs_constraints.get('linear_solver', True))

    #     f2, axes2 = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    #     modelPlot.data_plot(ax=axes2[0,0])
    #     modelPlot.model_plot(ax=axes2[0,1])
    #     modelPlot.normalized_residual_plot(ax=axes2[0,2], v_min=-6, v_max=6)
    #     modelPlot.source_plot(ax=axes2[1, 0], deltaPix_source=0.01, numPix=100)
    #     modelPlot.convergence_plot(ax=axes2[1, 1], v_max=1)
    #     modelPlot.magnification_plot(ax=axes2[1, 2])
    #     f2.tight_layout()
    #     f2.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    #     f2.savefig("/data/a.saricaoglu/lenstronomy/lenstronomy-data/mock_lenses/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M') +"/"+ im + "lens_models_pstprcsI_" + str(i) + "_" + current_time + ".png",  bbox_inches='tight')
    #     plt.close()
    #     #plt.show()

    #     f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    #     modelPlot.decomposition_plot(ax=axes[0,0], text='Lens light', lens_light_add=True, unconvolved=True)
    #     modelPlot.decomposition_plot(ax=axes[1,0], text='Lens light convolved', lens_light_add=True)
    #     modelPlot.decomposition_plot(ax=axes[0,1], text='Source light', source_add=True, unconvolved=True)
    #     modelPlot.decomposition_plot(ax=axes[1,1], text='Source light convolved', source_add=True)
    #     modelPlot.decomposition_plot(ax=axes[0,2], text='All components', source_add=True, lens_light_add=True, unconvolved=True)
    #     modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
    #     f.tight_layout()
    #     f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    #     f.savefig("/data/a.saricaoglu/lenstronomy/lenstronomy-data/mock_lenses/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M') +"/"+ im + "lens_models_pstprcsII_" + str(i) + "_" + current_time + ".png",  bbox_inches='tight')
    #     plt.close()
    #     #plt.show()
    #     print(kwargs_result)

    #     # the results of the MCMC chain

    #     sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[1]

    #     param_class = fitting_seq.param_class
    #     param_truths = param_class.kwargs2args(kwargs_lens=kwargs_lens, kwargs_source=kwargs_source, kwargs_lens_light=kwargs_lens_light)

    #     for j in range(len(chain_list)):
    #         ff, axx = chain_plot.plot_chain_list(chain_list, j)
    #         ff.savefig("/data/a.saricaoglu/lenstronomy/lenstronomy-data/mock_lenses/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M') +"/"+ im + "lens_models_chainplot_" + str(j) + "_" + str(i) +  "_" + current_time + ".png",  bbox_inches='tight')
    #         plt.close()
    #     print("number of non-linear parameters in the MCMC process: ", len(param_mcmc))
    #     print("parameters in order: ", param_mcmc)
    #     print("number of evaluations in the MCMC process: ", np.shape(samples_mcmc)[0])
    #     n_sample = len(samples_mcmc)
    #     print(n_sample)
    #     samples_mcmc_cut = samples_mcmc[int(n_sample*1/2.):]
    #     # if not samples_mcmc == []:
    #     n, num_param = np.shape(samples_mcmc_cut)

    #     plot = corner.corner(samples_mcmc_cut[:,:], labels=param_mcmc[:], show_titles=True, truths=param_truths, title_kwargs={"fontsize": 18}, figsize=(30,30))
    #     plot.savefig("/data/a.saricaoglu/lenstronomy/lenstronomy-data/mock_lenses/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M') +"/"+ im + "lens_models_corner_" + str(i) + "_" + current_time + ".png",  bbox_inches='tight')
    #     plt.close()
    #     print(f'Exposure time: {str(exp_time)} image min ave {np.mean(image_min)} image max ave {np.mean(image_max)} noisy image min ave  {np.mean(noisy_min)} noisy image max ave {np.mean(noisy_max)}')
hdulist.close()
