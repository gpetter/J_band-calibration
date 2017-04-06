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
import numpy
import fileinput
import time


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
large_aperture = []
flux_list = []
octplot = []

RA = []
Dec = []

# tune these values as necessary, one for each chip
class_star = float(.982)
seeing_fwhm = [0.4452, 0.4452, 0.477, .4452]
aperture_list = [0.6, .8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5,
                 5.2, 5.4, 5.6, 5.8, 6, 6.2, 6.4, 6.6, 6.8]

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
                "FLUX_APER(32)"], config={"MAG_ZEROPOINT": 29.1653, "PIXEL_SCALE": 0.159, "SEEING_FWHM": 0.4452,
                                          "FILTER_NAME": "gauss_4.0_7x7.conv", "CHECKIMAGE_TYPE": 'SEGMENTATION',
                                          "PHOT_APERTURES": (aperture_list)})
    for num in range(0, septlist.__len__()):
        out = sew(sept + septlist[num])
        fullname = sept + septlist[num].split(".fits")[0] + "_sexout.txt"
        septextracted.append(fullname)
        ascii.write(out["table"], fullname, overwrite=True)


def extractOct():
    for num in range(0, oct_list.__len__()):
        sew = sewpy.SEW(params=["NUMBER", "X_IMAGE", "Y_IMAGE", "FLAGS", "CLASS_STAR", "FWHM_IMAGE", "FLUX_APER(32)"],
                        config={"MAG_ZEROPOINT": 29.1653, "PIXEL_SCALE": 0.159, "SEEING_FWHM": seeing_fwhm[num],
                                "FILTER_NAME": "gauss_2.0_3x3.conv", "CHECKIMAGE_TYPE": 'SEGMENTATION', "CHECKIMAGE_NAME": ('c%s_check.fits' % (num + 1)), "PHOT_APERTURES": (aperture_list)})
        out = sew(oct + oct_list[num])
        name = oct + oct_list[num].split(".fits")[0]
        cleanname = name + "_clean.txt"
        sexname = name + "_sexout.txt"
        plotname = name + '_plot.txt'
        octclean.append(cleanname)
        oct_extracted.append(sexname)
        octplot.append(plotname)
        ascii.write(out["table"], cleanname, overwrite=True)
        copyfile(cleanname, sexname)

        for line in fileinput.input(cleanname, inplace=1):
            if (not line.startswith('N')) and (int(line.split(' ')[3]) == 0) and (
                float(line.split(' ')[1]) < 2250) and (
                        float(line.split(' ')[1]) > 100) and (float(line.split(' ')[2]) < 2250) and (
                        float(line.split(' ')[2]) > 100) and float(line.split(' ')[23] > 0):
                sys.stdout.write(line)
        fileinput.close()
        copyfile(cleanname, plotname)


def makeplot():
    for num in range(0, octplot.__len__()):
        thefile = octplot[num]
        f = open(thefile, 'r')
        lines = f.readlines()
        for x in lines:
            flux_aper.append((-2.5 * (math.log10(float(x.split(' ')[23])))))
            fwhm.append(x.split(' ')[5])
        # f.close()

        plt.figure(num)
        plt.plot(flux_aper, fwhm, 'ro')
        plt.title('c%s Flux vs. FWHM' % (num + 1))
        plt.xlabel('-2.5*log_10(flux)')
        plt.ylabel('FWHM')
        axes = plt.gca()
        axes.set_xlim([-14, -2])
        axes.set_ylim([0, 20])
        del flux_aper[:], fwhm[:]
        for line in fileinput.input(thefile, inplace=1):
            if line.startswith('N') == False and ((-2.5 * (math.log10(float(line.split(' ')[23])))) < -7) and (
                        (-2.5 * (math.log10(float(line.split(' ')[23])))) > -12):
                flux_aper.append((-2.5 * (math.log10(float(line.split(' ')[23])))))
                fwhm.append(float(line.split(' ')[5]))
                sys.stdout.write(line)
        fileinput.close()
        med = median(fwhm)
        sig = numpy.std(fwhm)
        print sig
        del flux_aper[:], fwhm[:]
        for line in fileinput.input(thefile, inplace=1):
            if line.startswith('N') == False and float(line.split(' ')[5]) < (med +  sig) and float(line.split(' ')[5]) > (
                    med - sig):
                flux_aper.append((-2.5 * (math.log10(float(line.split(' ')[23])))))
                fwhm.append(line.split(' ')[5])
                sys.stdout.write(line)
        fileinput.close()
        plt.plot(flux_aper, fwhm, 'bo')
        del flux_aper[:], fwhm[:]

        plt.savefig(('%s_flux vs. fwhm.png' % ('c' + str(num + 1))), bbox_inches='tight')
        # plt.show()
        plt.clf()
        plt.close()


def removeNearObjs():
    for num in range(0, oct_extracted.__len__()):
        # get xy for extracted sources
        thefile = oct_extracted[num]
        f = open(thefile, "r")
        lines = f.readlines()
        for x in lines:
            if not x.startswith("N"):
                octxlist.append(x.split(' ')[1])
                octylist.append(x.split(' ')[2])
        full = oct + 'c%s_xyoct_extracted.txt' % (num + 1)
        ascii.write([octxlist, octylist], full, names=('#x', 'y'), overwrite=True)
        del octxlist[:], octylist[:]
        str_out = oct + 'c%s_self_match.txt' % (num + 1)
        # match catalogue with itself, separation of 6''/(0.15923263854301059061''/pixel) = 37.6807170622 pixels
        iraf.xyxymatch(input=full, reference=full, output=str_out, tolerance=5, separation=39, verbose="no")
        f.close()
        f = open(str_out)
        lines = f.readlines()
        for x in lines:
            x = x.strip()
            if not x or x.startswith('#'):
                continue
            linelist.append(int(x.split()[5])-1)
        f.close()

        tempx = []
        tempy = []
        templineno = []
        f = open(thefile, "r")
        lines = f.readlines()
        i = 1
        for x in lines:
            if x.startswith("N"):
                continue
            elif i in linelist:
                templineno.append(x.split(' ')[0])
            i += 1
        del linelist[:]
        f.close()
        for line in fileinput.input(octplot[num], inplace=1):
            if line.split(' ')[0] in templineno:
                tempx.append(line.split(' ')[1])
                tempy.append(line.split(' ')[2])
                sys.stdout.write(line)
        fileinput.close()
        f = open(('c%s.xml' % (num + 1)), 'w')
        for y in range(0, tempx.__len__()):
            f.write('circle(%s,%s,%s)' % (tempx[y], tempy[y], 39) + '\n')
        f.close()


def findAperture():
    for num in range(0, octplot.__len__()):
        plt.figure(num)
        f = open(octplot[num], 'r')
        lines = f.readlines()
        for x in lines:
            for i in range(6, 38):
                flux_list.append((float(x.split(' ')[i])) / (float(x.split(' ')[37])))
            plt.plot(aperture_list, flux_list, 'ro')
            del flux_list[:]
        plt.title('F(<R)/F(6'') vs. R')
        axes = plt.gca()
        axes.set_ylim([0, 1.2])
        plt.xlabel('Aperture radius (asec)')
        plt.ylabel('Flux(<6 asec)/Flux(6 asec)')
        plotname = 'c%s_apertureplot.png' % (num + 1)
        plt.savefig(plotname, bbox_inches='tight')
        plt.clf()
        plt.close()


def matchOctSept():
    f = open(septextracted[0], "r")
    lines = f.readlines()
    tempx1, tempx2, tempx3, tempy1, tempy2, tempy3, templine1, templine2, templine3 = [], [], [], [], [], [], [], [], []
    for x in lines:
        if x.startswith('N') == False and int(x.split(' ')[3]) == 0 and float(x.split(' ')[4]) > .8 and \
                        int(x.split(' ')[1]) > 2200 and int(x.split(' ')[2]) > 2200:
            septxlist.append(x.split(' ')[1])
            septylist.append(x.split(' ')[2])
            septlineno.append('#%s' % x.split(' ')[0])
        elif x.startswith('N') == False and int(x.split(' ')[3]) == 0 and float(x.split(' ')[4]) > .8 and \
                        int(x.split(' ')[1]) > 2200 and int(x.split(' ')[2]) < 2200:
            tempx1.append(x.split(' ')[1])
            tempy1.append(x.split(' ')[2])
            templine1.append('#%s' % x.split(' ')[0])
        elif x.startswith('N') == False and int(x.split(' ')[3]) == 0 and float(x.split(' ')[4]) > .8 and \
                        int(x.split(' ')[1]) < 2200 and int(x.split(' ')[2]) < 2200:
            tempx2.append(x.split(' ')[1])
            tempy2.append(x.split(' ')[2])
            templine2.append('#%s' % x.split(' ')[0])
    for num in range(0, 4):
        full = sept + '%s_xysept_full.txt' % (num + 1)
        ascii.write([septxlist, septylist, septlineno], full, names=('#x', 'y', 'line no.'), overwrite=True)
    for num in range(0, octplot.__len__()):
        octfile = octplot[num]
        f = open(octfile, "r")
        lines = f.readlines()
        for x in lines:
            if x.startswith('N') == False and int(x.split(' ')[3]) == 0 and float(x.split(' ')[4]) > .8:
                octxlist.append(x.split(' ')[1])
                octylist.append(x.split(' ')[2])
                octlineno.append('#%s' % x.split(' ')[0])
        octxy = oct + 'c%s_xyoct_full.txt' % (num + 1)
        ascii.write([octxlist, octylist, octlineno], octxy, names=('#x', 'y', 'line no.'), overwrite=True)
        del octxlist[:], octylist[:], octlineno[:]
        str_out = 'c%s_oct_sept_match.txt' % (num + 1)
        iraf.xyxymatch(input=octxy, reference=(sept + 'xysept_full.txt'), output=str_out, tolerance=500, verbose="no")
        f = open(str_out, "r")
        lines = f.readlines()
        for x in lines:
            x = x.strip()
            if not x or x.startswith('#'):
                continue
            septmatchlinelist.append(int(x.split()[4]) - 1)
            octmatchlinelist.append(int(x.split()[5]) - 1)
        f.close()
        # septmatchlinelist.sort(), octmatchlinelist.sort()
        f = open((sept + 'xysept_full.txt'), 'r')
        lines = f.readlines()
        i = 1
        for x in lines:
            if x.startswith("#"):
                continue
            elif i in septmatchlinelist:
                septcatlinelist.append(int(x.split('#')[1]))
            i += 1
        f.close()
        f = open((octxy), 'r')
        lines = f.readlines()
        i = 1
        for x in lines:
            if x.startswith("#"):
                continue
            elif i in octmatchlinelist:
                octcatlinelist.append(int(x.split('#')[1]))
            i += 1
        f.close()
        # septcatlinelist.sort(), octcatlinelist.sort()
        f = open(septextracted[0], 'r')
        lines = f.readlines()
        for x in lines:
            if x.startswith("N"):
                continue
            elif int(x.split(' ')[0]) in septcatlinelist:
                RA.append(float(x.split(' ')[5]))
                Dec.append(float(x.split(' ')[6]))
        f.close()
        f = open(octfile, 'r')
        lines = f.readlines()
        k = 0
        for x in lines:
            if int(x.split(' ')[0]) in octcatlinelist:
                oct_total.append('%s %s %s\n' % (x.rstrip('\n'), RA[k], Dec[k]))
                k += 1
        del septmatchlinelist[:], octmatchlinelist[:], septcatlinelist[:], octcatlinelist[:], RA[:], Dec[:]
        f.close()
    f = open('oct_total.txt', 'w')
    for item in oct_total:
        f.write('%s' % item)


def twomassmatch():
    twomassRA = []
    twomassDec = []
    mag = []
    lineno = []
    octRA = []
    octDec = []
    f = open('2mass.txt', 'r')
    lines = f.readlines()
    i = 1
    for x in lines:
        if x.startswith('#') == False:
            twomassRA.append(x.split()[0])
            twomassDec.append(x.split()[1])
            mag.append("#" + x.split()[6])
            lineno.append("#" + str(i))
            i += 1

    f.close()
    ascii.write([twomassRA, twomassDec, mag, lineno], '2mass_RADec.txt', names=('#RA', 'Dec', 'Mag', 'Line number'),
                overwrite=True)
    f = open('oct_total.txt', 'r')
    lines = f.readlines()
    for x in lines:
        octRA.append(x.split(' ')[34])
        octDec.append(x.split(' ')[35])
    ascii.write([octRA, octDec], 'Oct_RADec.txt', names=('#RA', 'Dec'), overwrite=True)
    iraf.xyxymatch(input='Oct_RADec.txt', reference='2mass_RADec.txt', output="2mass_Oct_matched.txt", tolerance=30,
                   verbose="no")


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


def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
            return None
    if len(lst) %2 == 1:
            return lst[((len(lst)+1)/2)-1]
    else:
            return float(sum(lst[(len(lst)/2)-1:(len(lst)/2)+1]))/2.0


extractOct()
extractSept()
makeplot()
removeNearObjs()
findAperture()
#matchOctSept()
#twomassmatch()
