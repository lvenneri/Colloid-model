from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import os
import matplotlib.patches as patches
import cv2
from os.path import join as pjoin  # tool to make paths
import imutils
import pandas as pd
import scipy.io


def bin_image(im, factor, interp_down='INTER_AREA', interp_up='INTER_LINEAR', maintain_dim=False):
    interp_dic = {"INTER_LINEAR": cv2.INTER_LINEAR,
                  "INTER_CUBIC": cv2.INTER_CUBIC,
                  "INTER_AREA": cv2.INTER_AREA,
                  "INTER_NEAREST": cv2.INTER_NEAREST,
                  "INTER_LANCZOS4": cv2.INTER_LANCZOS4,
                  "INTER_LINEAR_EXACT": cv2.INTER_LINEAR_EXACT}

    im = cv2.resize(im, (0, 0), fx=1.0 / float(factor[0]), fy=1.0 / float(factor[1]),
                    interpolation=interp_dic[interp_down])
    if maintain_dim:
        im = cv2.resize(im, (0, 0), fx=float(factor[0]), fy=float(factor[1]),
                        interpolation=interp_dic[interp_up])

    return im


def noise_scaling(ims, save_loc,startTime, sampling_IR=[20, 5], height_fraction=.5, pixel_reduction_factor=16,):
    """
    This script takes the fts file from the camera and processes by
    - median/average of each burst
    - selecting ROI w/ rotation and cropping (approximate)
    - saving the image 3d matrix, and saving the lineouts to a csv with the time index

    :param ims: array of images, [images, x dim,y dim]
    :param save_loc: where to save, probably same directory fine
    :param startTime: datetime from the labview csv file in the same directory as this fts file
    :param sampling_IR: [ # frames, time between spurts in minutes]
    :param height_fraction: fraction of horizontal pixels to use in the thin film array, about the center. good way
    to get rid of border effects
    :param pixel_reduction_factor: binning in the z direction. if you have 300 pixels, and set a factor of 10,
    will have 30 pixels via averaging of adjecent pixels

    :return:
    """

    # number of whole bursts in the data set
    whole_bursts = int((ims.shape[0] / sampling_IR[0]))
    whole_bursts_end_idx = int((ims.shape[0] / sampling_IR[0])) * sampling_IR[0]
    # capture time starting with 0 for each burst, in ms
    time = np.linspace(0, 1, whole_bursts) * whole_bursts * sampling_IR[1] * 60 * 1000

    time = pd.to_datetime(time, unit='ms',origin=startTime[0])  # DD:HH:MM:SS.fff
    print(time)
    # reshape array be [burst time, burst length, im x dimensions, im y dimensions]
    ims = ims[0:whole_bursts_end_idx].reshape(-1, sampling_IR[0], ims.shape[1], ims.shape[2])

    # take a median across the images (prefer median or sigma clipped mean to avoid problems from skipped frames (
    # ~1/10000)
    ims = np.median(ims, axis=(1))

    force_select = False
    if force_select or not os.path.isfile(save_loc + '.npz'):
        fig, ax = plt.subplots(figsize=(10, 5))
    global coords
    coords = []

    def onclick(event):
        global ix, iy
        ix, iy = event.xdata, event.ydata

        global coords
        coords.append((event.xdata, event.ydata))
        rect = patches.Rectangle((event.xdata, event.ydata), 1, 1, linewidth=.5,
                                 edgecolor='r',
                                 facecolor='none')
        # Add the patch to the Axes
        ax.add_patch(rect)
        plt.draw()

        if len(coords) == 4:
            import time
            fig.savefig(save_loc + '.png')
            np.savez(save_loc, coords=coords)

            time.sleep(1)

            fig.canvas.mpl_disconnect(cid)

            plt.close()

        return coords

    # s0_ratio, s0_rang, s0_o1, s0_o2

    if not force_select and os.path.isfile(save_loc + '.npz'):
        coords = np.load(save_loc + '.npz')['coords']
        grab_ROI = 1
    else:
        grab_ROI = 0

    if 1 == 1:  # SHITE IMAGE CONDITION
        if grab_ROI == 0:  # get coordinates
            dispim = ims[0, :, :]
            vmin, vmax = np.nanpercentile(dispim.ravel(), 5), np.nanpercentile(dispim.ravel(), 95)
            imaa = ax.imshow(dispim, vmin=vmin, vmax=vmax, origin='lower', interpolation='nearest')  # ,cmap='jet')
            ax.axis('off')
            plt.colorbar(imaa, orientation='horizontal', aspect=30, fraction=.06, shrink=1.0, pad=0.0,
                         label='Counts')

            plt.title('Select line path and height: Left - Right - Bottom - Top')
            cid = fig.canvas.mpl_connect('button_press_event', onclick)

            plt.show()

            # plt.title('Select 4 points (bottom side), Box Edge = ' + str(box * 2 + 1))

            # global coords  # contains the coordinates at which to eval
            # print coords
            grab_ROI = 1

        # extract analytics at those coordinates and add to array
        # cycle through ims and through rois
        im_lineouts = []
        im_matrix = []
        for i in range(ims.shape[0]):
            im = ims[i, :, :]

            # rotate the image
            angle = np.arcsin((coords[0][1] - coords[1][1]) / (coords[0][0] - coords[1][0])) / np.pi * 180
            # then crop the middle section to 50% of the height coords[4][1]:coords[3][1]
            scale_cm_per_pixel = .6521 / (coords[3][1] - coords[2][1])

            height = int((coords[3][1] - coords[2][1]) * height_fraction / 2)

            im_rot = imutils.rotate_bound(im, angle)[int(coords[2][1]) + height:int(coords[3][1]) - height,
                     int(coords[0][0]):int(coords[1][0])]

            scale_cm_per_pixel_binned = .6521*height_fraction/ im_rot.shape[0]*pixel_reduction_factor
            # print(im.shape)
            # print(im_rot.shape)
            # print(height)

            im_matrix.append(im_rot)
            im_rot = bin_image(im_rot, factor=[pixel_reduction_factor, im_rot.shape[0]])
            im_lineouts.append(im_rot)

        im_matrix = np.stack(im_matrix, axis=2)
        im_lineouts = np.vstack(im_lineouts)
        plt.imshow(im_lineouts, interpolation="nearest")

        # grab the z-direction distnace along the flow path
        # 0,1,2,3,4 / N * length
        # 0,1,2,3,4 * length/pixel; length/pixel =
        pixel_spatial = np.round(np.arange(im_lineouts.shape[1]) * scale_cm_per_pixel_binned, 3)
        print("_____ PIXEL SPATIAL _____")
        print(pixel_spatial)
        data_table = pd.DataFrame(im_lineouts, columns=pixel_spatial, index=time)

    # save data points including data_loc in csv
    data_table.to_csv(save_loc + '.csv', index=True)
    scipy.io.savemat(save_loc + '.mat', mdict={'source': im_matrix})

    return data_table

def grabfits(datapath,startTtime_ms=0):
    # startTime_ms, is the initial true time in ms, which you might grab from the LabView data.

    file_loc = [pjoin(datapath, f) for f in os.listdir(datapath) if f.endswith('.fts')][0]
    csv_loc = [f for f in os.listdir(datapath) if f.endswith('.csv') and not f.startswith('Rec-') and not f.startswith('coll')][0]

    save_loc = file_loc.split('.fts')[0]

    # fits_image_filename = fits.util.get_testdata_filepath('test0.fits')
    hdul = fits.open(file_loc)
    data = hdul[0].data
    print('Import size', data.shape)

    average_vals = np.mean(data, axis=(1, 2))

    # average_vals2 = np.mean(data_2, axis=(1, 2))


    startTime = pd.to_datetime(pd.DataFrame({'year': [csv_loc.split('_')[0]],
                       'month': [csv_loc.split('_')[1]],
                       'day': [csv_loc.split('_')[2]],
                              'hour': [csv_loc.split('_')[3]],
                              'minute': [csv_loc.split('_')[4]],
                              'second': [csv_loc.split('_')[5].split('.')[0]]
                       }))
    print("START TIME", startTime[0])

    data_table = noise_scaling(data, save_loc, startTime=startTime,sampling_IR=[20, 5], height_fraction=.5, pixel_reduction_factor=16)

    # if 1 == 0:
    #     plt.figure()
    #     plt.subplot(211)
    #     plt.plot(average_vals)
    #     plt.subplot(212)
    #     plt.plot(average_vals2)
    #
    #     plt.subplot(121)
    #     plt.imshow(data[0, :, :])
    #     plt.subplot(122)
    #     plt.imshow(data_1[0, 0, :, :] - data[0, :, :])
    #     plt.figure()
    #     plt.subplot(211)
    #     plt.plot(average_vals)
    #     plt.subplot(212)
    #     plt.plot(average_vals2)
    #     plt.show()


    hdul.close()


if __name__ == "__main__":
    datapath ="/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210801_isoparaffin_5W"
    grabfits(datapath)
    plt.show()

