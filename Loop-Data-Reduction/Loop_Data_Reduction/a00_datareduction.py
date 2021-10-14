from a01_IR_camera_data_import import grabfits
from a02_data_interlacer import interlace_main
from a03_data_process import data_proces

"""
Mother script to do sequential processing on an IR + csv from the Lubrizol Flow Loop for a whole set of specified 
data folders in data_paths or just a few.
"""
if __name__ == "__main__":

    datapath = "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210930_Water_5W"

    # run single dataset
    if 1 == 1:
        grabfits(datapath)
        interlace_main(datapath)
        data_proces(datapath)

    # or run multiple datasets
    elif 1 == 1:


        datapaths = ["/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210726_isoparaffin_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210727_isoparaffin_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210728_isoparaffin_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210729_isoparaffin_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210730_isoparaffin_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210801_isoparaffin_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210802_isoparaffin_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210804_isoparaffin_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210805_OS622339H_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210806_OS622339H_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210808_OS622339H_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210809_OS625475_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210810_OS625475_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210811_OS625475_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210812_OS625475_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210813_OS625475_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210814_OS625475_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210815_OS625475_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210817_isoparaffin_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210824_isoparaffin_5W",
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210826_isoparaffin_5W"
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210924_water_5W"
                     "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210925_water_5W"

                     ]

        for datapath in datapaths:
            print('____________________________________________' + datapath)
            # grabfits(datapath)
            interlace_main(datapath)
            data_proces(datapath)
