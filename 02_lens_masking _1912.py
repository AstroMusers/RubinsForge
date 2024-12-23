# This version chooses 20 cutouts.
from lsst.daf.butler import Butler
import os
from lsst.meas.algorithms.maskStreaks import MaskStreaksTask
from lsst.meas.algorithms.maskStreaks import MaskStreaksTask
import lsst.afw.geom as afwGeom
from lsst.pipe.tasks.calexpCutout import CalexpCutoutTask
from lsst.obs.base import FitsRawFormatterBase
from lsst.afw.fits import Fits
from lsst.geom import Box2I, Point2I, Extent2I
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import read_lens as gs
import datetime as dt
import streaks as ads
import cutouts as ct

matplotlib.rcParams['font.size'] = 19

c = dt.datetime.now()
if not os.path.exists("/data/a.saricaoglu/lsst_pipeline/Files/" +  str(c.strftime("%m.%d"))+ "/" + c.strftime('%H%M')): 
    directoryf = "/data/a.saricaoglu/lsst_pipeline/Files/"  +  str(c.strftime("%m.%d"))+ "/" + c.strftime('%H%M')
    os.makedirs(directoryf) 
if not os.path.exists("/data/a.saricaoglu/lsst_pipeline/Plots/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M')): 
    directoryp = "/data/a.saricaoglu/lsst_pipeline/Plots/"  +  str(c.strftime("%m.%d"))+ "/" + c.strftime('%H%M')
    os.makedirs(directoryp) 
directoryf = "/data/a.saricaoglu/lsst_pipeline/Files/"  +  str(c.strftime("%m.%d"))
directoryp = "/data/a.saricaoglu/lsst_pipeline/Plots/"  +  str(c.strftime("%m.%d"))


# repo_path = os.path.join(os.environ['RC2_SUBSET_DIR'], 'SMALL_HSC')
# butler = Butler(repo_path)

# registry = butler.registry
# for col in registry.queryCollections():
#     print(col)
# for ref in registry.queryDatasets('raw', collections='HSC/raw/all', instrument='HSC'):
#     print(ref.dataId)

collection = f"u/{os.environ['USER']}/single_frame"
butler = Butler('SMALL_HSC', collections=collection, instrument='HSC')

# Initialize the MaskStreaks task
config = MaskStreaksTask.ConfigClass()
#config.clusterMinimumDeviation = 5              #def. = 2  -Allowed deviation (in pixels) from a straight line for a detected line",
#config.clusterMinimumSize = 500                #def. = 50 -Allowed deviation (in pixels) from a straight line for a detected line",
config.nSigma = 2                      #def. = 2  -Number of sigmas from center of kernel to include in voting procedure
config.invSigma = 0.1             #def. = 10.**-1 -"Inverse of the Moffat sigma parameter (in units of pixels) describing the profile of the streak",
#config.dChi2Tolerance = 0.5      
mask_streaks_task = MaskStreaksTask(config=config)

registry = butler.registry
# for col in registry.queryCollections():
#     print(col)

detected = 0
detected_original = 0
n = 1
for ref in butler.registry.queryDatasets('calexp', physical_filter='HSC-G'):

    print(ref.dataId)
    # data_id = {"instrument": "HSC", "visit": 11694, "detector": 50}  # Example data ID
    # try:
    calexp_image = butler.get("calexp",  dataId=ref.dataId)

    mask_results_or =  mask_streaks_task.find(calexp_image.maskedImage)
    mask_or = mask_results_or.mask
    if np.sum(mask_or) != 0 :
        detected_original = detected_original + 1

    image_array = calexp_image.maskedImage.image.array
    width = calexp_image.getWidth()
    height = calexp_image.getHeight()
    center_x = width // 2
    center_y = height // 2
    center = Point2I(center_x, center_y)

    cutouts_factory= []
    cutout_size = 100  # Each cutout is 100x100 pixels

    # Define the offsets for the top 6 boxes relative to the center
    top_offsets = [

        (-800, 1400),  # Top-left
        (0, 1600),     # Top-center
        (800, 1400),   # Top-right
        (-800, -1400),  # Top-left
        (0, -1600),     # Top-center
        (800, -1400),   # Top-right
        (-600, 1100),  # Top-left
        (0, 1200),     # Top-center
        (600, 1100),   # Top-right
        (-600, -1100),  # Top-left
        (0, -1200),     # Top-center
        (600, -1100),   # Top-right
        (-400, 700),  # Top-left
        (0, 800),     # Top-center
        (400, 700),   # Top-right
        (-400, -700),  # Top-left
        (0, -800),     # Top-center
        (400, -700),   # Top-right
        (-200, 200),  # Top-left
        (0, 300),     # Top-center
        (200, 200),   # Top-right
        (-200, -200),  # Top-left
        (0, -300),     # Top-center
        (200, -200),   # Top-right
    ]

    # Define the center of the full image (assumed to be at pixel (1500, 1500))

    # Create bounding boxes for the top 6 cutouts
    bboxes = []
    for dx, dy in top_offsets:
        # Calculate the lower-left corner of the cutout
        lower_left = Point2I(center.getX() + dx - cutout_size // 2,
                                center.getY() + dy - cutout_size // 2)
        # Create a bounding box of size 100x100 pixels
        bbox = Box2I(lower_left, Extent2I(cutout_size, cutout_size))
        bboxes.append(bbox)


    
    f, ax = plt.subplots(2, 3, figsize=(28, 19), sharex=False, sharey=False)
    im0 = ax[0][0].imshow(image_array, origin='lower', cmap='viridis', vmin=np.percentile(image_array, 5), vmax=np.percentile(image_array, 95))
    plt.colorbar(im0,label='Pixel Value')
    ax[0][0].set_title('Original Image Data \n (calexp.maskedImage.array)')
    ax[0][0].set_xlabel('X Pixel')
    ax[0][0].set_ylabel('Y Pixel')

    print('image min max', np.min(calexp_image.maskedImage.image.array), np.max(calexp_image.maskedImage.image.array))

    streak = ads.get_streak(4176,2048,10e6)
    streaky_image = calexp_image.maskedImage.image.array + streak
    calexp_image.maskedImage.image.array[:,:] = streaky_image
    streaky_mask = calexp_image.maskedImage.mask.array + streak
    calexp_image.maskedImage.mask.array[:,:] = streaky_mask    
    #image_array = cutout.image.array + streak 
    im1 = ax[0][1].imshow(calexp_image.maskedImage.image.array, origin='lower', cmap='viridis', vmin=np.percentile(image_array, 5), vmax=np.percentile(image_array, 95))
    plt.colorbar(im1, label='Pixel Value')
    ax[0][1].set_title('Streak Added \n (calexp.maskedImage.image.array)')
    ax[0][1].set_xlabel('X Pixel')
    ax[0][1].set_ylabel('Y Pixel')
    # plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/preProcessed_lens_sim2" + str(n)+ ".png")
    mask_results = mask_streaks_task.find(calexp_image.maskedImage)
    
    mask = mask_results.mask

    if np.sum(mask) != 0 :
        detected = detected + 1

    config.detectedMaskPlane = 'STREAK'
    
    masked = mask_streaks_task.run(calexp_image.maskedImage)
    mask_plane_dict = calexp_image.mask.getMaskPlaneDict()
    streak_bit = 1 << mask_plane_dict['STREAK']
    mask2 = (calexp_image.maskedImage.mask.array & streak_bit) != 0
    masked_image_array1 = calexp_image.maskedImage.image.array.copy()
    masked_image_array1[mask2] = np.nan 
    masked_image_array2 = np.ma.masked_array(calexp_image.maskedImage.image.array,mask)
    im2 = ax[0][2].imshow(mask, origin='lower', cmap='viridis')
    plt.colorbar(im2, label='Pixel Value')
    ax[0][2].set_title('Streak Mask \n (mask_results.mask, after "find")')
    ax[0][2].set_xlabel('X Pixel')
    ax[0][2].set_ylabel('Y Pixel')
    im3 = ax[1][0].imshow(mask2, origin='lower', cmap='viridis')
    plt.colorbar(im3, label='Pixel Value')
    ax[1][0].set_title('Streak Mask \n (calexp.maskedImage.mask, after "run")')
    ax[1][0].set_xlabel('X Pixel')
    ax[1][0].set_ylabel('Y Pixel')    
    im4 = ax[1][1].imshow(masked_image_array2, origin='lower', cmap='viridis', vmin=np.percentile(image_array, 5), vmax=np.percentile(image_array, 95))
    plt.colorbar(im4, label='Pixel Value')
    ax[1][1].set_title('Streak Masked \n (np.ma.(calexp.maskedImage.image.array,mask))')
    ax[1][1].set_xlabel('X Pixel')
    ax[1][1].set_ylabel('Y Pixel')
    im5 = ax[1][2].imshow(masked_image_array1, origin='lower', cmap='viridis',vmin=np.percentile(image_array, 5), vmax=np.percentile(image_array, 95))
    plt.colorbar(im5, label='Pixel Value ')
    ax[1][2].set_title('Streak Masked \n  (calexp.maskedImage.image.array)')
    ax[1][2].set_xlabel('X Pixel')
    ax[1][2].set_ylabel('Y Pixel')
    
    j = 1    # Overlay the bounding boxes

    for bbox in bboxes:
        ax[1][2].add_patch(plt.Rectangle((bbox.getMinX(), bbox.getMinY()),
                                        bbox.getWidth(), bbox.getHeight(),
                                        edgecolor='red', facecolor='none', lw=2))
        x_min = bbox.getMinX()
        y_min = bbox.getMinY()
        width = bbox.getWidth()
        height = bbox.getHeight()

        upper_right_x = x_min + width
        upper_right_y = y_min + height

        # Add the label in the upper-right corner
        ax[1][2].text(upper_right_x, upper_right_y, str(j), color='red', fontsize=15,
            ha='right', va='bottom', fontweight='bold') 
        j = j + 1
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/preProcessed_lens_sim3" + str(n)+ ".png",bbox_inches='tight')


    f2, ax2 = plt.subplots(4, 6, figsize=(20, 13), sharex=True, sharey=True)
    # Define the center (in pixel coordinates) and size of the cutout
    print(type(streak), np.shape(streak))
    print(type(image_array), np.shape(image_array))


    #cutouts_factory= ct.extract_cutouts(calexp_image, top_offsets,100)
    cutouts_method= ct.extract_cutouts_with_task(n, calexp_image, top_offsets, cutout_size).getMaskedImages()

    print('cutout len ',len(cutouts_method))
    cutouts_factory= []
    m = 0
    lens, noisy_lens = gs.get_simulated_lens(n)
    k = 0
    p = 1
    for bbox in bboxes:

        cutout = calexp_image.Factory(calexp_image, bbox)
        cutouts_factory.append(cutout)
        lensed_image = cutout.maskedImage.image.array + lens
        cutout.maskedImage.image.array[:,:] = lensed_image
        lensed_mask = cutout.maskedImage.mask.array + lens
        cutout.maskedImage.mask.array[:,:] = lensed_mask 
        im1 = ax2[k][m].imshow(cutout.maskedImage.image.array, cmap='viridis', origin='lower')
        ax2[k][m].set_title(f'Cutout {p}')
        #plt.colorbar(im1, label='Pixel Value ', ax=ax2[k][m], shrink=0.9)
        m = m + 1
        p = p + 1
        if p == 7:
            k = 1
            m = 0
        if p == 13:
            k = 2
            m = 0
        if p == 19:
            k = 3
            m = 0
    #ax2[0][3].axis('off')
    #f2.suptitle('Cutouts from the Full Image')
    f2.colorbar(im1, label='Pixel Value ', ax=ax2.ravel().tolist(),  shrink=1)
    f2.savefig(directoryp +  "/" + c.strftime('%H%M') + "/lens_cutouts" + str(n)+ ".png" , bbox_inches='tight')

    f4, ax4 = plt.subplots(4, 6, figsize=(20, 13), sharex=True, sharey=True)
    m = 0
    k = 0
    p = 1
    for cutout in cutouts_method:        
        lensed_image = cutout.image.array + lens
        cutout.image.array[:,:] = lensed_image
        lensed_mask = cutout.mask.array + lens
        cutout.mask.array[:,:] = lensed_mask 
        mask_plane_dict = cutout.mask.getMaskPlaneDict()
        streak_bit = 1 << mask_plane_dict['STREAK']
        mask2 = (cutout.mask.array & streak_bit) != 0
        cutout.image.array[mask2] = np.nan
        if np.sum(mask2) != 0:
            print('cutout streak detected in cutout no ', p)
        im1 = ax4[k][m].imshow(cutout.image.array, cmap='viridis', origin='lower')
        ax4[k][m].set_title(f'Cutout {p}')
        #plt.colorbar(im1, label='Pixel Value ', ax=ax2[k][m], shrink=0.9)
        cutout.writeFits(directoryf +  "/" + c.strftime('%H%M') + "/lens_cutouts" + str(n) + "_cut_" + str(k))
        m = m + 1
        p = p + 1
        if p == 7:
            k = 1
            m = 0
        if p == 13:
            k = 2
            m = 0
        if p == 19:
            k = 3
            m = 0

    f4.colorbar(im1, label='Pixel Value ', ax=ax4.ravel().tolist(),  shrink=1)
    f4.savefig(directoryp +  "/" + c.strftime('%H%M') + "/lens_cutouts_masked" + str(n)+ ".png" , bbox_inches='tight')      


    f3, ax3 = plt.subplots(1, 2, figsize=(11, 4), sharex=True, sharey=True,gridspec_kw={'width_ratios': [1, 1]})  
    im0 = ax3[0].imshow(lens, cmap='viridis', origin='lower')
    im1 = ax3[1].imshow(noisy_lens, cmap='viridis', origin='lower')
    f3.suptitle('Injected Lens')
    f3.colorbar(im1, label='Pixel Value ', ax=ax3.ravel().tolist(),  shrink=0.9)
    f3.savefig(directoryp +  "/" + c.strftime('%H%M') + "/lens_injected" + str(n)+ ".png" , bbox_inches='tight')


    # if np.sum(mask_results.mask) != 0:
    #     f, ax = plt.subplots(1, 2, figsize=(25, 8), sharex=False, sharey=False)
    #     im0 = ax[0][0]imshow(image_array, origin='lower', cmap='viridis', vmin=np.percentile(image_array, 5), vmax=np.percentile(image_array, 95))
    #     plt.colorbar(im0,label='Pixel Value')
    #     ax[0][0]set_title('LSST Image Data')
    #     ax[0][0]set_xlabel('X Pixel')
    #     ax[0][0]set_ylabel('Y Pixel')
    #     #mask = mask_results.mask
    #     m = mask_streaks_task.run(calexp_image)
    #     masked = calexp_image.image.array
    #     #masked_image_array = cutout.mask.getArray()
    #     # Inspect the results
    #     #   print("Number of streaks masked:", mask_results.mask)
    #     # plt.figure(figsize=(10, 10))
    #     # plt.imshow(image_array, origin='lower', cmap='viridis', vmin=np.percentile(image_array, 5), vmax=np.percentile(image_array, 95))
    #     # plt.colorbar(label='Pixel Value')
    #     # plt.title('LSST Image Data')
    #     # plt.xlabel('X Pixel')
    #     # plt.ylabel('Y Pixel')
    #     # plt.show()
    #     im1 = ax[0][1].imshow(masked, origin='lower',cmap='viridis', vmin=np.percentile(image_array, 5), vmax=np.percentile(image_array, 95))
    #     plt.colorbar(im1, label='Pixel Value')
    #     ax[0][1].set_title('LSST Mask')
    #     ax[0][1].set_xlabel('X Pixel')
    #     ax[0][1].set_ylabel('Y Pixel')

    #     # im2 = ax[0][2].imshow(m.mask.getArray(), origin='lower', cmap='viridis', vmin=np.percentile(image_array, 5), vmax=np.percentile(image_array, 95))
    #     # plt.colorbar(im2, label='Pixel Value')
    #     # ax[0][2].set_title('LSST Masked Image')
    #     # ax[0][2].set_xlabel('X Pixel')
    #     # ax[0][2].set_ylabel('Y Pixel')
    #     plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Processed_lens_sim" + str(n)+ ".png")
        #plt.show()

    n = n + 1
    print('detected streaks ', detected, 'detected streaks in original image ', detected_original)
    if n == 5:
           exit()
    # # except:
    #    print('no image')
