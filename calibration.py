##
# name: calibration.py
#
# purpose: tbd
#
# date last edited: 3/23/17
#
# author: Grayson Petter
#
##

from astropy.io import ascii
from astropy.table import Table
import os
import sewpy
import sys
from pyraf import iraf
import matplotlib.pyplot as plt
import math
from shutil import copyfile

if len(sys.argv) > 1:
    (oct_list_filename, sept_list_filename) = (sys.argv[1], sys.argv[2])
else:
    oct_list_filename = 'oct.txt'
    sept_list_filename = 'sept.txt'
oct_list_filename = 'oct.txt'
sept_list_filename = 'sept.txt'

cur_dir = os.getcwd()
sept = cur_dir + "/Sept/"
oct = cur_dir + "/Oct/"

# list of october fits filenames
oct_list = []

# contains october Sextractor catalogues
oct_extracted = []

# edited october catalogue
octclean = []
septlist = []
septextracted = []

# lists of xy coordinates used for xyxymatch
octxlist = []
octylist = []
septxlist = []
septylist = []
septmatchlinelist = []
octmatchlinelist = []
septcatlinelist = []
octcatlinelist = []
oct_total = []
octlineno = []
septlineno = []
linelist = []


RA = []
Dec = []

# tune these values as necessary, one for each chip
class_stars = [.96, .967, .96, .98]
seeing_fwhm = [.954, .556, .795,.2067]
aperture_list = [0.6, .8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5,
                5.2, 5.4, 5.6, 5.8, 6]

flux_aper = []
fwhm = []

# create work paths
if not os.path.exists(sept):
    os.makedirs(sept)
if not os.path.exists(oct):
    os.makedirs(oct)

fimg = open(oct_list_filename, "r")
lines = fimg.readlines()
for x in lines:
    x = x.replace('\n', '')
    if x.endswith('.fits'):
        oct_list.append(x)
fimg.close()

fimg = open(sept_list_filename, "r")
lines = fimg.readlines()
for x in lines:
    x = x.replace('\n', '')
    if x.endswith('.fits'):
        septlist.append(x)
fimg.close()


def extractSept():
    sew = sewpy.SEW(
        params=["NUMBER", "X_IMAGE", "Y_IMAGE", "FLAGS", "CLASS_STAR", "ALPHA_SKY", "DELTA_SKY", "FWHM_IMAGE",
                "FLUX_APER(28)"], config={"MAG_ZEROPOINT": 29.1653, "PIXEL_SCALE": 0.159, "SEEING_FWHM": 0.7795454545,
                                          "FILTER_NAME": "gauss_4.0_7x7.conv", "PHOT_APERTURES": (aperture_list)})
    for num in range(0, septlist.__len__()):
        out = sew(sept + septlist[num])
        fullname = sept + septlist[num].split(".fits")[0] + "_sexout.txt"
        septextracted.append(fullname)
        ascii.write(out["table"], fullname, overwrite=True)


def extractOct():

    for num in range(0, oct_list.__len__()):
        sew = sewpy.SEW(params=["NUMBER", "X_IMAGE", "Y_IMAGE", "FLAGS", "CLASS_STAR", "FWHM_IMAGE", "FLUX_APER(28)"],
                        config={"MAG_ZEROPOINT": 29.1653, "PIXEL_SCALE": 0.159, "SEEING_FWHM": seeing_fwhm[num],
                                "FILTER_NAME": "gauss_2.0_3x3.conv", "PHOT_APERTURES": (aperture_list)})
        out = sew(oct + oct_list[num])
        name = oct + oct_list[num].split(".fits")[0]
        cleanname = name + "_clean.txt"
        sexname = name + "_sexout.txt"
        octclean.append(cleanname)
        oct_extracted.append(sexname)
        ascii.write(out["table"], cleanname, overwrite=True)
        copyfile(cleanname, sexname)




def makeplot():
    for num in range(0, octclean.__len__()):
        f = open(octclean[num], "r+")
        lines = f.readlines()
        f.seek(0)
        for x in lines:
            if x.startswith('N') == False and float(x.split(' ')[4]) > class_stars[num] and float(
                    x.split(' ')[23]) > 0 and int(x.split(' ')[3]) == 0:
                flux_aper.append(x.split(' ')[23])
                fwhm.append(x.split(' ')[5])
                f.write(x)
        f.truncate()
        f.close()

        for i in range(0, flux_aper.__len__()):
            flux_aper[i] = -2.5 * (math.log10(float(flux_aper[i])))
        plt.figure(num)
        plt.plot(flux_aper, fwhm, 'ro')
        plt.title('c%s Flux vs. FWHM' % (num+1))
        plt.xlabel('-2.5*log_10(flux)')
        plt.ylabel('FWHM')
        plt.savefig(('%s_flux vs. fwhm.png' % ('c' + str(num + 1))), bbox_inches='tight')
        plt.show()
        del flux_aper[:]
        del fwhm[:]
        plt.close()



def removeNearObjs():
    for num in range(0, octclean.__len__()):
        # get xy for extracted sources
        thefile = octclean[num]
        f = open(thefile, "r")
        lines = f.readlines()
        for x in lines:
            octxlist.append(x.split(' ')[1])
            octylist.append(x.split(' ')[2])
        full = oct + 'c%s_xyoct_clean.txt' % (num + 1)
        ascii.write([octxlist, octylist], full, names=('#x', 'y'), overwrite=True)
        del octxlist[:]
        del octylist[:]
        str_out = oct + 'c%s_self_match.txt' % (num + 1)
        # match catalogue with itself, separation of 6''/(0.15923263854301059061''/pixel) = 37.6807170622 pixels
        iraf.xyxymatch(input=full, reference=full, output=str_out, tolerance=5, separation=37.6807170622, verbose="no")
        f.close()
        f = open(str_out)
        lines = f.readlines()
        for x in lines:
            x = x.strip()
            if not x or x.startswith('#'):
                continue
            linelist.append(int(x.split()[5]))
        f.close()
        f = open(thefile, "r+")
        lines = f.readlines()
        f.seek(0)
        i = 1
        for x in lines:
            if i in linelist:
                f.write(x)
            i += 1
        f.truncate()
        del linelist[:]
        f.close()

def findAperture():
    for num in range(0, octclean.__len__()):
        f = open(octclean[num], 'r')
        lines = f.readlines()
        for x in lines:
            x=x
def matchOctSept():
    for num in range(0, septextracted.__len__()):
        f = open(septextracted[num], "r")
        lines = f.readlines()
        for x in lines:
            if x.startswith('N') == False and int(x.split(' ')[3]) == 0 and float(x.split(' ')[4]) > .7:
                septxlist.append(x.split(' ')[1])
                septylist.append(x.split(' ')[2])
                septlineno.append('#%s' % x.split(' ')[0])
        full = sept + 'xysept_full.txt'
        ascii.write([septxlist, septylist, septlineno], full, names=('#x', 'y', 'line no.'), overwrite=True)
    for num in range(0, oct_extracted.__len__()):
        f = open(oct_extracted[num], "r")
        lines = f.readlines()
        for x in lines:
            if x.startswith('N') == False and int(x.split(' ')[3]) == 0 and float(x.split(' ')[4]) > .7:
                octxlist.append(x.split(' ')[1])
                octylist.append(x.split(' ')[2])
                octlineno.append('#%s' % x.split(' ')[0])
        octxy = oct + 'c%s_xyoct_full.txt' % (num + 1)
        ascii.write([octxlist, octylist, octlineno], octxy, names=('#x', 'y', 'line no.'), overwrite=True)
        del octxlist[:]
        del octylist[:]
        del octlineno[:]
        str_out = 'c%s_oct_sept_match.txt' % (num + 1)
        iraf.xyxymatch(input=octxy, reference=(sept + 'xysept_full.txt'), output=str_out, tolerance=30, verbose="no")
        f = open(str_out, "r")
        lines = f.readlines()
        for x in lines:
            x = x.strip()
            if not x or x.startswith('#'):
                continue
            septmatchlinelist.append(int(x.split()[4]))
            octmatchlinelist.append(int(x.split()[5]))
        f.close()
        septmatchlinelist.sort(), octmatchlinelist.sort()
        f = open((sept+'xysept_full.txt'), 'r')
        lines = f.readlines()
        i = 1
        for x in lines:
            if i in septmatchlinelist:
                septcatlinelist.append(int(x.split('#')[1]))
            i += 1
        f.close()
        f = open((octxy), 'r')
        lines = f.readlines()
        i = 1
        for x in lines:
            if i in octmatchlinelist:
                octcatlinelist.append(int(x.split('#')[1]))
            i += 1
        f.close()
        septcatlinelist.sort(), octcatlinelist.sort()
        f = open(septextracted[0], 'r')
        lines = f.readlines()
        i = 1
        for x in lines:
            if i in septcatlinelist:
                RA.append(float(x.split(' ')[5]))
                Dec.append(float(x.split(' ')[6]))
            i += 1
        f.close()
        #print linelist, linelist2, RA, Dec
        f = open(oct_extracted[num], 'r')
        lines = f.readlines()
        j, k = 1, 0
        for x in lines:
            if x.startswith("N"):
                continue
            elif j in octcatlinelist:
                oct_total.append('%s %s %s\n' % (x.rstrip('\n'), RA[k], Dec[k]))
                k += 1
            j += 1
        del septmatchlinelist[:], octmatchlinelist[:], septcatlinelist[:], octcatlinelist[:], RA[:], Dec[:]
        f.close()
    f = open('oct_total.txt', 'w')
    for item in oct_total:
        f.write('%s\n' % item)


def extractionComment():
    print "First, Sextractor is run on September data. The pixel scale is given as 0.159 asec/pixel, and the FWHM is " \
          "measured to be 4.9 pixels. Therefore the SEEING_FWHM is calculated to be 0.779545 asec. This value along " \
          "with the zeropoint estimate is used to configure Sextractor. The parameters retrieved are x and y" \
          "coordinates, flags, star class, RA & Dec, FWHM_Image, and fluxes from apertures of size 0.4'' to 6 '', " \
          "counting by 0.2'' \n \n The same is done on October data, though RAs and Decs are excluded. Catalogs " \
          "are written to files for both October and September data, though FWHM is measured to be 2.4 pixels in " \
          "October data. FWHM is adjusted accordingly. \n"


def plotComment():
    print "Now sources are selected from the October catalogues which throw no flags, and have a CLASS_STAR  " \
          "above a certain value specified. The flux of that source through the 4'' aperture and the FWHM are " \
          "stored. Then for each chip, -2.5log_10(flux) is plotted against FWHM. Looking at this plot, observe " \
          "whether or not your CLASS_STAR value selected objects on the stellar sequence. Adjust as needed."
extractSept()
extractOct()
makeplot()
removeNearObjs()
matchOctSept()