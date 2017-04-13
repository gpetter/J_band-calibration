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
import astropy.stats as st
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

x_origin, y_origin = 0, 0

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
septsectionslist =[]
stars = []

RA = []
Dec = []

# tune these values as necessary, one for each chip
class_star = float(.982)
seeing_fwhm = [0.4452, 0.4452, 0.4452, .4452]
aperture_list_asec = [0, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6,
                      4.8, 5, 5.2, 5.4, 5.6, 5.8, 6, 6.2, 6.4, 6.8]
aperture_list_pix = [i*6.285 for i in aperture_list_asec]
aperture_list_asec_radius = [float(i)/2 for i in aperture_list_asec]
aperture_list_asec_radius.pop(0)
flux_aper = []
fwhm = []
# create work paths
if not os.path.exists(sept):
    os.makedirs(sept)
if not os.path.exists(oct):
    os.makedirs(oct)

with open(oct_list_filename, "r") as fimg:
    lines = fimg.readlines()
    for x in lines:
        x = x.replace('\n', '')
        if x.endswith('.fits'):
            oct_list.append(x)

with open(sept_list_filename, "r") as fimg:
    lines = fimg.readlines()
    for x in lines:
        x = x.replace('\n', '')
        if x.endswith('.fits'):
            septlist.append(x)


def extractSept():
    sew = sewpy.SEW(
        params=["NUMBER", "X_IMAGE", "Y_IMAGE", "FLAGS", "CLASS_STAR", "ALPHA_SKY", "DELTA_SKY", "FWHM_IMAGE",
                "FLUX_APER(%s)" % (len(aperture_list_pix))], config={"MAG_ZEROPOINT": 29.1653, "PIXEL_SCALE": 0.159,
                                    "SEEING_FWHM": 0.4452, "FILTER_NAME": "gauss_3.0_7x7.conv",
                                    "CHECKIMAGE_TYPE": 'SEGMENTATION', "PHOT_APERTURES": (aperture_list_pix)})
    for num in range(0, septlist.__len__()):
        out = sew(sept + septlist[num])
        fullname = sept + septlist[num].split(".fits")[0] + "_sexout.txt"
        septextracted.append(fullname)
        ascii.write(out["table"], fullname, overwrite=True)


def extractOct():
    for num in range(0, oct_list.__len__()):
        sew = sewpy.SEW(params=["NUMBER", "X_IMAGE", "Y_IMAGE", "FLAGS", "CLASS_STAR", "FWHM_IMAGE",
                                "FLUX_APER(%s)" % (len(aperture_list_pix))],
                        config={"MAG_ZEROPOINT": 29.1653, "PIXEL_SCALE": 0.159, "SEEING_FWHM": seeing_fwhm[num],
                                "FILTER_NAME": "gauss_3.0_7x7.conv", "CHECKIMAGE_TYPE": 'SEGMENTATION',
                                "CHECKIMAGE_NAME": ('c%s_check.fits' % (num + 1)),
                                "PHOT_APERTURES": (aperture_list_pix)})
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
        temp1 = []
        temp2 = []
        for line in fileinput.input(cleanname, inplace=1):
            x = line.split(' ')
            if (not line.startswith('N')) and (int(x[3]) == 0) and (
                130 < float(x[1]) < 2250) and (130 < float(x[2]) < 2250) and (float(x[23]) > 0):
                temp1.append(x[1])
                temp2.append(x[2])
                sys.stdout.write(line)
        fileinput.close()
        copyfile(cleanname, plotname)
        with open(('%soct.xml') % (num+1), 'w') as f:
            for y in range(0, temp1.__len__()):
                f.write('circle(%s,%s,%s)' % (temp1[y], temp2[y], 30) + '\n')


def makeplot():
    for num in range(0, octplot.__len__()):
        thefile = octplot[num]
        with open(thefile, 'r') as f:
            lines = f.readlines()
            for x in lines:
                flux_aper.append((-2.5 * (math.log10(float(x.split(' ')[23])))))
                fwhm.append(x.split(' ')[5])

        plt.figure(num)
        plt.plot(flux_aper, fwhm, 'ro')
        plt.title('c%s Flux vs. FWHM' % (num + 1))
        plt.xlabel('-2.5*log_10(flux)')
        plt.ylabel('FWHM')
        axes = plt.gca()
        axes.set_xlim([-15, -2])
        axes.set_ylim([0, 20])
        del flux_aper[:], fwhm[:]
        for line in fileinput.input(thefile, inplace=1):
            if line.startswith('N') == False and ((-2.5 * (math.log10(float(line.split(' ')[23])))) < -8.5) and (
                        (-2.5 * (math.log10(float(line.split(' ')[23])))) > -12):
                flux_aper.append((-2.5 * (math.log10(float(line.split(' ')[23])))))
                fwhm.append(float(line.split(' ')[5]))
                stars.append(line.split(' ')[0])
                sys.stdout.write(line)
        fileinput.close()
        med = median(fwhm)
        sig = st.biweight_midvariance(fwhm)
        print sig
        del flux_aper[:], fwhm[:]
        for line in fileinput.input(thefile, inplace=1):
            if line.startswith('N') == False and float(line.split(' ')[5]) < (med +  sig) and \
                            float(line.split(' ')[5]) > (med - sig):
                flux_aper.append((-2.5 * (math.log10(float(line.split(' ')[23])))))
                fwhm.append(line.split(' ')[5])
                sys.stdout.write(line)
        fileinput.close()
        plt.plot(flux_aper, fwhm, 'bo')
        del flux_aper[:], fwhm[:]

        plt.savefig(('%s_flux vs. fwhm.png' % ('c' + str(num + 1))), bbox_inches='tight')
        plt.show()
        plt.clf()
        plt.close()


def removeNearObjs():
    for num in range(0, oct_extracted.__len__()):
        # get xy for extracted sources
        thefile = oct_extracted[num]
        with open(thefile, "r") as f:
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
        iraf.xyxymatch(input=full, reference=full, output=str_out, tolerance=5, separation=45, verbose="no")

        with open(str_out) as f:
            lines = f.readlines()
            for x in lines:
                x = x.strip()
                if not x or x.startswith('#'):
                    continue
                linelist.append(int(x.split()[5])-1)

        tempx = []
        tempy = []
        templineno = []
        with open(thefile, "r") as f:
            lines = f.readlines()
            i = 1
            for x in lines:
                if x.startswith("N"):
                    continue
                elif i in linelist:
                    templineno.append(x.split(' ')[0])
                i += 1
        del linelist[:]
        for line in fileinput.input(octplot[num], inplace=1):
            if line.split(' ')[0] in templineno:
                tempx.append(line.split(' ')[1])
                tempy.append(line.split(' ')[2])
                sys.stdout.write(line)
        fileinput.close()
        with open(('c%s.xml' % (num + 1)), 'w') as f:
            for y in range(0, tempx.__len__()):
                f.write('circle(%s,%s,%s)' % (tempx[y], tempy[y], 45) + '\n')

def findAperture():

    for num in range(0, octplot.__len__()):
        plt.figure(num)
        with open(octplot[num], 'r')as f:
            lines = f.readlines()
            for x in lines:
                for i in range(7, (len(aperture_list_pix)+6)):
                    flux_list.append((float(x.split(' ')[i])) / (float(x.split(' ')[len(aperture_list_asec)+5])))
                plt.plot(aperture_list_asec_radius, flux_list, 'ro')
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
    with open(septextracted[0], "r") as f:
        lines = f.readlines()
        tempx1, tempx2, tempx3, tempy1, tempy2, tempy3, templine1, templine2, templine3, RA1, RA2, RA3, RA4, Dec1, Dec2, Dec3, Dec4 = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
        for x in lines:
            line = x.split(' ')
            if (not x.startswith("N")) and int(line[3]) == 0 and float(line[13]) > 3000 and\
                        4300 > float(line[1]) > 2300 and 4300 > float(line[2]) > 2300:
                septxlist.append(line[1])
                septylist.append(line[2])
                septlineno.append('#%s' % line[0])
                RA1.append('#%s' % line[5])
                Dec1.append('#%s' % line[6])

            elif (not x.startswith("N")) and int(line[3]) == 0 and \
                        4300 > float(line[1]) > 2300 and 130 < float(line[2]) < 2300:
                tempx1.append(line[1])
                tempy1.append(line[2])
                templine1.append('#%s' % line[0])
                RA2.append('#%s' % line[5])
                Dec2.append('#%s' % line[6])
            elif (not x.startswith("N")) and int(line[3]) == 0 and \
                       130 < float(line[1]) < 2300 and 130 < float(line[2]) < 2300:
                tempx2.append(line[1])
                tempy2.append(line[2])
                templine2.append('#%s' % line[0])
                RA3.append('#%s' % line[5])
                Dec3.append('#%s' % line[6])
            elif (not x.startswith("N")) and int(line[3]) == 0 and \
                        130 < float(line[1]) < 2300 and 4300 > float(line[2]) > 2300:
                tempx3.append(line[1])
                tempy3.append(line[2])
                templine3.append('#%s' % line[0])
                RA4.append('#%s' % line[5])
                Dec4.append('#%s' % line[6])

    with open(('match1.xml'), 'w') as f:
        for y in range(0, septxlist.__len__()):
            f.write('circle(%s,%s,%s)' % (septxlist[y], septylist[y], 30) + '\n')

    for num in range(0,4):
        septsectionslist.append((sept + '%s_xysept_full.txt' % (num + 1)))
    ascii.write([septxlist, septylist, septlineno, RA1, Dec1], septsectionslist[0], names=('#x', 'y', 'line no.', 'RA', 'Dec'), overwrite=True)
    ascii.write([tempx1, tempy1, templine1, RA2, Dec2], septsectionslist[1], names=('#x', 'y', 'line no.', 'RA', 'Dec'), overwrite=True)
    ascii.write([tempx2, tempy2, templine2, RA3, Dec3], septsectionslist[2], names=('#x', 'y', 'line no.', 'RA', 'Dec'), overwrite=True)
    ascii.write([tempx3, tempy3, templine3, RA4, Dec4], septsectionslist[3], names=('#x', 'y', 'line no.', 'RA', 'Dec'), overwrite=True)

    for num in range(0, octclean.__len__()):
        octfile = octclean[num]

        with open(octfile, "r") as f:
            lines = f.readlines()
            for x in lines:
                line = x.split(' ')
                if int(line[3]) == 0:
                    octxlist.append(x.split(' ')[1])
                    octylist.append(x.split(' ')[2])
                    octlineno.append('#%s' % x.split(' ')[0])

        octxy = oct + 'c%s_xyoct_full.txt' % (num + 1)
        ascii.write([octxlist, octylist, octlineno], octxy, names=('#x', 'y', 'line no.'), overwrite=True)
        del octxlist[:], octylist[:], octlineno[:]
        str_out = 'c%s_oct_sept_match.txt' % (num + 1)

        if num == 0:
            x_origin, y_origin = 2150, 2150
        elif num == 1:
            x_origin, y_origin = 2150, 0
        elif num == 2:
            x_origin, y_origin = 0, 0
        else:
            x_origin, y_origin = 0, 2150

        iraf.xyxymatch(input=(septsectionslist[num]), reference=octxy, xin=x_origin, yin=y_origin, output=str_out, tolerance=30, verbose="no", separation=55)
        with open(str_out, "r") as f:
            lines = f.readlines()
            for x in lines:
                x = x.strip()
                if not x or x.startswith('#'):
                    continue
                septmatchlinelist.append(int(x.split()[5]) - 1)
                octmatchlinelist.append(int(x.split()[4]) - 1)

        if len(septmatchlinelist) != len(set(septmatchlinelist)) or len(octmatchlinelist) != len(set(octmatchlinelist)):
            print "Error: tweak xyxymatch parameters so that a coordinate in one set of data isn't matched to multiple coordinates in the other. First try increasing separation value"

        with open((septsectionslist[num]), 'r') as f:
            lines = f.readlines()
            length = len(lines)
            i = 1
            while i != length:
                x = lines[i]
                if x.startswith("#"):
                    continue
                elif len(septmatchlinelist) > 0 and i == septmatchlinelist[0]:
                    septcatlinelist.append(int(x.split('#')[1]))
                    RA.append(float(x.split('#')[2]))
                    Dec.append(float(x.split('#')[3]))
                    septmatchlinelist.pop(0)
                    i = 1
                    continue
                i += 1

        with open((octxy), 'r') as f:
            lines = f.readlines()
            length = len(lines)
            i = 1
            while i != length:
                x = lines[i]
                if x.startswith("#"):
                    continue
                elif len(octmatchlinelist) > 0 and i == octmatchlinelist[0]:
                    #if num == 0:
                        #print i, x
                    octcatlinelist.append(int(x.split('#')[1]))
                    octmatchlinelist.pop(0)
                    i = 1
                    continue
                i += 1
        with open(octfile, 'r') as f:
            lines = f.readlines()
            length = len(lines)
            k = 0
            i = 0
            while (i) != length:
                x = lines[i]
                if x.startswith("#"):
                    continue
                elif len(octcatlinelist) > 0 and int(x.split(' ')[0]) == octcatlinelist[0]:
                    oct_total.append('%s %s %s %s\n' % (x.rstrip('\n'), RA[k], Dec[k], (num+1)))
                    octcatlinelist.pop(0)
                    i = 0
                    k += 1
                    continue
                i += 1
        del septmatchlinelist[:], octmatchlinelist[:], septcatlinelist[:], octcatlinelist[:], RA[:], Dec[:]
    with open('oct_total.txt', 'w') as f:
        for item in oct_total:
            f.write('%s' % item)


def twomassmatch():
    twomassRA = []
    twomassDec = []
    mag = []
    lineno = []
    octRA = []
    octDec = []
    octlineno = []
    octflux = []
    chipno = []
    with open('2mass.txt', 'r') as f:
        lines = f.readlines()
        i = 1
        for x in lines:
            if x.startswith('#') == False:
                twomassRA.append(float(x.split()[0])*1000)
                twomassDec.append(float(x.split()[1])*1000)
                mag.append("#" + x.split()[6])
                lineno.append("#" + str(i))
                i += 1

    ascii.write([twomassRA, twomassDec, mag, lineno], '2mass_RADec.txt', names=('#RA', 'Dec', 'Mag', 'Line number'),
                overwrite=True)
    with open('oct_total.txt', 'r') as f:
        lines = f.readlines()
        for x in lines:
            line = x.split(' ')
            octRA.append(float(line[len(line)-3])*1000)
            octDec.append(float(line[len(line)-2])*1000)
            octlineno.append(line[0])
            octflux.append(line[len(line)-4])
            chipno.append(line[len(line)-1])

    ascii.write([octRA, octDec, octlineno, chipno, octflux], 'Oct_RADec.txt', names=('#RA', 'Dec', 'line', 'chip', 'flux'), overwrite=True)
    iraf.xyxymatch(input='Oct_RADec.txt', reference='2mass_RADec.txt', output="2mass_Oct_matched.txt", tolerance=3,
                   verbose="no", xref=0, yref=0, xin=0, yin=0)

def calculateZP():
    twommatch = []
    octmatch = []
    Mag = []
    flux = []
    octlineno = []
    chipno = []
    with open("2mass_Oct_matched.txt", 'r') as f:
        lines = f.readlines()
        for x in lines:
            x = x.strip()
            line = x.split()
            if not x or x.startswith('#'):
                continue
            twommatch.append((int(line[4])-1))
            octmatch.append((int(line[5])-1))

    if len(twommatch) != len(set(twommatch)) or len(octmatch) != len(set(octmatch)):
        print "Error: tweak xyxymatch parameters so that a coordinate in one set of data isn't matched to multiple coordinates in the other. First try increasing separation"

    with open('2mass_RADec.txt', 'r') as f:
        lines = f.readlines()
        length = len(lines)
        i = 1
        while i != length:
            x = lines[i]
            if x.startswith("#"):
                continue
            elif len(twommatch) > 0 and i == twommatch[0]:
                Mag.append(x.split('#')[1])
                twommatch.pop(0)
                i = 1
                continue
            i += 1
    with open('Oct_RADec.txt', 'r') as f:
        lines = f.readlines()
        length = len(lines)
        i = 1
        while i != length:
            x = lines[i]
            if x.startswith("#"):
                continue
            elif len(octmatch) > 0 and i == octmatch[0]:
                flux.append(x.split(' ')[4].replace('\n', ''))
                octlineno.append(x.split(' ')[2])
                chipno.append(x.split(' ')[3])
                octmatch.pop(0)
                i = 1
                continue
            i += 1
    zp = []
    for num in range(0, len(Mag)):
        if float(flux[num]) < 1 or (octlineno[num] not in stars):
            continue
        zp.append(float(Mag[num]) + float(2.5*(math.log10(float(flux[num])))))
    plt.figure(30)
    plt.hist(zp)
    axes = plt.gca()
    axes.set_xlim([25, 28])
    axes.set_ylim([0, 10])
    plt.show()



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

