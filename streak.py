import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import random
import os
from datetime import datetime
from astropy.io import fits

def min_max_normalize(arr):
    """Normalize the array using Min-Max normalization."""
    arr_min = np.min(arr)
    arr_max = np.max(arr)
    normalized_arr = (arr - arr_min) / (arr_max - arr_min)
    return normalized_arr



def create_streak(xdim, ydim, amplitude):
    while True:    
        # Step 1: Create the xdim x ydim grid
        maxval = max(xdim,ydim)
        grid_size_x = int(maxval) 
        grid_size_y = int(maxval)
        # Streak amplitude
        Z = np.zeros((grid_size_x, grid_size_y))
        x = np.linspace(0,xdim,grid_size_x)
        y = np.linspace(0,ydim,grid_size_y)
        X, Y = np.meshgrid(x,y)
        # print(np.shape(X))
        # print(np.shape(Y))
        # Step 2: Define random center point on the grid
        x_center = random.randint(grid_size_x//2-xdim//4, grid_size_x//2+xdim//4)
        y_center = random.randint(grid_size_y//2-ydim//4, grid_size_y//2+ydim//4)
        x_shift = random.randint(-xdim//4, xdim//4)
        y_shift = random.randint(-ydim//4, ydim//4)
        #print(x_center, y_center)
        # Step 3: Create the infinite line with a width of 3 pixels
        line_width = 1
        angle = np.random.randint(0,360)  # You can change the angle here

        Yp =  rotate(Y.transpose(), angle, reshape=False)
        Xp =  rotate(X, angle, reshape=False)
        # Xp = X.transpose()
        # Yp = Y
        xsrc = x_center
        ysrc = y_center
        Z = np.exp(-0.5*(((Yp - ysrc)**2 + (Xp - xsrc)**2)/(0.025*300)**2))   
        # # Create an infinite line by setting pixels along the x-axis and a width
        # for i in range(grid_size_x):
        #     for w in range(-grid_size_x, grid_size_x):
        #         if 0 <= grid_size_x + w < grid_size_x and 0 <= i < grid_size_x:
        #             Z[grid_size_x + w,i] = np.exp(-0.5*((Y[i, grid_size_x+w]-y[grid_size_x + w])**2)/(0.1/np.sqrt(8))**2)    # Draw the line horizontally at the center
        #             # previously np.exp(-0.5*((Y[i, grid_size_x+w]-y[grid_size_x + w])**2)/(0.025*3/np.sqrt(8))**2) 
        #             #print(Z[i, x_center + w])
        # Step 4: Rotate the line by a certain angle around the center point
        # x_shift = random.randint(-xdim//4, xdim//4)
        # y_shift = random.randint(-ydim//4, ydim//4)
        Z = Z[grid_size_x//2-xdim//2:grid_size_x//2+xdim//2 + 1,grid_size_y//2-ydim//2: grid_size_y//2+ydim//2]
        Z = min_max_normalize(Z)

        if np.max(Z) == 1:
            break
    Z = Z * amplitude

    # f, ax = plt.subplots(1, 3, figsize=(25, 8), sharex=False, sharey=False)    


    # im1 = ax[1].imshow(mask, origin='lower', cmap='gray')
    # plt.colorbar(im1, label='Pixel Value')
    # ax[1].set_title('LSST Mask')
    # ax[1].set_xlabel('X Pixel')
    # ax[1].set_ylabel('Y Pixel')
    # plt.scatter(x_center, y_center)
    # plt.imshow(Z[grid_size_x//2-xdim//2:grid_size_x//2+xdim//2 - 1,grid_size_y//2-ydim//2: grid_size_y//2+ydim//2 - 1], cmap='gray', origin='lower')
    # plt.show()
    # Move the grid so that the line center is at the origin, rotate, then move it back
    #Z_rotated = rotate(Z, angle, reshape=False)
    #Z_shifted = np.roll(Z_rotated, shift=(x_shift, y_shift), axis=(1, 0))
    #print('streak min max,', np.min(Z), np.max(Z))

    # # # Step 5: Plot the result using plt.imshow()
    # plt.figure(figsize=(6, 6))
    # # # plt.xlim((50,150))
    # # # plt.ylim((50,150))

    # plt.scatter(x_center, y_center)
    # plt.imshow(Z, cmap='gray', origin='lower')
    # # # plt.imshow(Z_rotated, cmap='gray', origin='lower')
    # plt.title(f'Streak Rotated by {angle}')

    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.colorbar(label='Streak Intensity')
    # # # plt.xticks(np.linspace(50, 150, 6), np.linspace(0, 100, 6).astype(int))  # 6 ticks from 0 to 100
    # # # plt.yticks(np.linspace(50, 150, 6), np.linspace(0, 100, 6).astype(int))  # 6 ticks
    # # # Show the plot
    # plt.show()
    #print(np.shape(Z))

    return Z
# s = create_streak(4176,2048,1)
# print(np.max(s), np.min(s))
def save_streak(Z):
    # Displays Time
    #In this part, we create and write on the .fits files.
    file = "/data/a.saricaoglu/lsst_pipeline/streaks.fits"
    #PrimaryHDU will have the image data, with image related headers. Data is empty for now.
    hdu_pr = fits.PrimaryHDU()
    #ImageHDU is where we keep the mock lens data, with lens related headers.
    # Write to a FITS file
    hdu_pr.writeto(file, overwrite=True)
    # Verify the file by reading it back
    hdulist = fits.open(file, mode='append')
    image_hdu = fits.ImageHDU(Z)
    hdulist.append(image_hdu)
    hdulist.close()

def create_and_save(n, xdim, ydim, amplitude):
    for i in range(0,n):
        s = create_streak(xdim, ydim, amplitude)
        save_streak(s)

def get_streak(xdim, ydim, amplitude):
    file = "/data/a.saricaoglu/lsst_pipeline/streaks.fits"
    if os.path.exists(file) :
        hdu = fits.open(file)
        streak = hdu[random.randint(0, len(hdu))].data
        return streak
    else:
        print('No streak file exists, creating streaks...')
        print('Streak simulation start; ',datetime.now().strftime("%d.%m.%y, %H:%M"))
        create_and_save(1000, xdim, ydim, amplitude)
        print('Streak simulation end; ',datetime.now().strftime("%d.%m.%y, %H:%M"))
        hdu = fits.open(file)
        streak = hdu[random.randint(0, len(hdu))].data
        return streak


