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
import numpy as np
import astropy.stats as st
import fileinput
from astropy import units as u
from astropy.coordinates import SkyCoord
from itertools import combinations
from itertools import repeat

cur_dir = os.getcwd()
sept = cur_dir + "/Sept/"
oct = cur_dir + "/Oct/"
userdata = cur_dir + '/Userdata/'
sex = cur_dir + '/Sex/'
reference = cur_dir + '/Reference/'

if len(sys.argv) > 1:
    (oct_list_filename, sept_list_filename) = (oct + sys.argv[1], sept + sys.argv[2])
else:
    oct_list_filename = oct + 'oct.txt'
    sept_list_filename = sept + 'sept.txt'
oct_list_filename = oct + 'oct.txt'
sept_list_filename = sept + 'sept.txt'
bandline = 0
indexjhk = 0

# parameters: edit these
############################################################################

# band: 'j' or 'k'
band = 'j'
# mode, run or display writeup
writeup = True

# scale: arcsec/pixel, measure with imexam
scale = 0.159
# fwhm for oct and sept image: in arcsec
seeing_fwhm_oct = .4452
seeing_fwhm_sept = .7791
# lower and upper bounds of pixels to cut. use to remove false sources from the artifacts near the edge of the image
loweroct = 130
upperoct = 2250
lowersept = 130
uppersept = 4300

############################################################################

if band == 'j':
    bandline = 6
elif band == 'k':
    bandline = 14
    indexjhk = 2

# lists containing filenames of all files needed
oct_list, oct_extracted, octclean, septlist, septextracted, septclean = [], [], [], [], [], []
masteroctlines, masterseptlines = [[], [], [], []], [[], [], [], []]

# lists of xy coordinates used for xyxymatch
octplot, septsectionslist, liststars, stars, septstars, RA, Dec = [], [], [], [], [], [], []

octzp, octzperr = 0, 0

# tune these values as necessary, one for each chip
aperture_list_asec = [0, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4,
                      4.6, 4.8, 5, 5.2, 5.4, 5.6, 5.8, 6, 6.2, 6.4, 6.6]

aperture_list_pix = [i * (1 / scale) for i in aperture_list_asec]
aperture_list_asec_radius = [float(i) / 2 for i in aperture_list_asec]
aperture_list_asec_radius.pop(0)
convergeap, convergesept = 0, 0
separation = (1 / scale) * aperture_list_asec[len(aperture_list_asec) - 1] + 15
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


# run Sextractor on september data and each of the four October images.
# make copies to edit and truncate
def extractSept():
    sew = sewpy.SEW(
        params=["NUMBER", "X_IMAGE", "Y_IMAGE", "FLAGS", "CLASS_STAR", "ALPHA_SKY", "DELTA_SKY", "FWHM_IMAGE",
                "FLUX_APER(%s)" % (len(aperture_list_pix))], config={"PIXEL_SCALE": scale,
                                                                     "SEEING_FWHM": seeing_fwhm_sept,
                                                                     "FILTER_NAME": (sex + "gauss_5.0_9x9.conv"),
                                                                     "CHECKIMAGE_TYPE": 'SEGMENTATION',
                                                                     "PHOT_APERTURES": (aperture_list_pix),
                                                                     "CHECKIMAGE_NAME": (userdata + 'sept_check.fits')})
    for num in range(len(septlist)):
        out = sew(sept + septlist[num])
        sexname = sept + septlist[num].split(".fits")[0] + "_sexout.txt"
        cleanname = sept + septlist[num].split(".fits")[0] + "_clean.txt"
        septextracted.append(sexname)
        septclean.append(cleanname)
        ascii.write(out["table"], sexname, overwrite=True)
        copyfile(sexname, cleanname)
        for line in fileinput.input(cleanname, inplace=1):
            x = line.split(' ')
            if (not line.startswith('N')) and (float(x[26]) > 0) and int(x[3]) == 0 and float(x[4]) > .1:
                sys.stdout.write(line)


def extractOct():
    for num in range(len(oct_list)):
        sew = sewpy.SEW(params=["NUMBER", "X_IMAGE", "Y_IMAGE", "FLAGS", "CLASS_STAR", "FWHM_IMAGE",
                                "FLUX_APER(%s)" % (len(aperture_list_pix))],
                        config={"PIXEL_SCALE": scale, "SEEING_FWHM": seeing_fwhm_oct,
                                "FILTER_NAME": (sex + "gauss_3.0_7x7.conv"), "CHECKIMAGE_TYPE": 'SEGMENTATION',
                                "CHECKIMAGE_NAME": (userdata + 'c%s_check.fits' % (num + 1)),
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
        temp1, temp2 = [], []
        for line in fileinput.input(cleanname, inplace=1):
            x = line.split(' ')
            if (not line.startswith('N')) and (int(x[3]) == 0) and (
                            loweroct < float(x[1]) < upperoct) and (loweroct < float(x[2]) < upperoct) and (
                float(x[24]) > 0):
                temp1.append(x[1])
                temp2.append(x[2])
                sys.stdout.write(line)
        copyfile(cleanname, plotname)
        with open((userdata + '%soct.xml') % (num + 1), 'w') as f:
            for y in range(len(temp1)):
                f.write('circle(%s,%s,%s)' % (temp1[y], temp2[y], 30) + '\n')


# tune values to permit what can be counted as a star
# range in x axis which selects sources on the horizontal branch of graph
minlogflux, maxlogflux = -13.5, -7
# coeffecient of biweight midvariance sigma from the median to select sources, one for each chip
sigmamultiple = [2.5, 5, 2, 5]
# if data is bad and horizontal branch is not clearly defined, (has a significant spread in y-axis), use this.
# this is the maximum FWHM value that will be used, so that sources which fall above this line are excluded
# from median and BWMV calculation
maxfwhm = 3.4


# create plot of all sources' -2.5log(flux_4asec) vs FWHM. Distinguish those sources which lie on a horizontal
# branch of the diagram, indicating star-like
def makeplot(month):
    flux_aper, fwhm = [], []
    if month == 'oct':
        for num in range(len(octplot)):
            thefile = octplot[num]
            with open(thefile, 'r') as f:
                lines = f.readlines()
                for x in lines:
                    line = x.split(' ')
                    flux_aper.append((-2.5 * (math.log10(float(line[24])))))
                    fwhm.append(line[5])

            plt.figure(num)
            plt.plot(flux_aper, fwhm, 'ro')
            plt.title('c%s Flux vs. FWHM' % (num + 1)), plt.xlabel('-2.5*log_10(flux)'), plt.ylabel('FWHM')
            axes = plt.gca()
            axes.set_xlim([-16, -2]), axes.set_ylim([0, 20])
            del flux_aper[:], fwhm[:]
            for line in fileinput.input(thefile, inplace=1):
                x = line.split(' ')
                if (not line.startswith('N')) and (minlogflux < (-2.5 * (math.log10(float(x[24])))) < maxlogflux) and \
                                float(x[5]) < maxfwhm:
                    flux_aper.append((-2.5 * (math.log10(float(x[24])))))
                    fwhm.append(float(x[5]))
                    sys.stdout.write(line)
            med = median(fwhm)
            sig = st.biweight_midvariance(fwhm)
            del flux_aper[:], fwhm[:]
            for line in fileinput.input(thefile, inplace=1):
                x = line.split(' ')
                if (not line.startswith('N')) and ((med - sigmamultiple[num] * sig) < float(x[5]) <
                                                       (med + sigmamultiple[num] * sig)):
                    flux_aper.append((-2.5 * (math.log10(float(x[24])))))
                    fwhm.append(x[5])
                    sys.stdout.write(line)

            plt.plot(flux_aper, fwhm, 'bo')
            del flux_aper[:], fwhm[:]

            plt.savefig(('%s_flux vs. fwhm.png' % ('c' + str(num + 1))), bbox_inches='tight')
            if writeup == True:
                plt.show()
            plt.clf()
            plt.close()
    elif month == 'sept':
        thefile = septclean[0]
        with open(thefile, 'r') as f:
            lines = f.readlines()
            for x in lines:
                line = x.split(' ')
                flux_aper.append((-2.5 * (math.log10(float(line[26])))))
                fwhm.append(line[7])
                if float(line[7]) > 80:
                    print line[7]
                    print line[0]
        plt.figure(1000)
        plt.plot(flux_aper, fwhm, 'ro')
        plt.title('Sept Flux vs. FWHM'), plt.xlabel('-2.5*log_10(flux)'), plt.ylabel('FWHM')
        axes = plt.gca()
        # axes.set_xlim([-15, -2]), axes.set_ylim([0, 20])
        del flux_aper[:], fwhm[:]
        for line in fileinput.input(thefile, inplace=1):
            x = line.split(' ')
            if (not line.startswith('N')) and (-15 < (-2.5 * (math.log10(float(x[26])))) < -10):
                flux_aper.append((-2.5 * (math.log10(float(x[26])))))
                fwhm.append(float(x[7]))
                sys.stdout.write(line)
        med = median(fwhm)
        sig = st.biweight_midvariance(fwhm)
        del flux_aper[:], fwhm[:]
        for line in fileinput.input(thefile, inplace=1):
            x = line.split(' ')
            if ((med - 1.5 * sig) < float(x[7]) < (med + 1.5 * sig)):
                flux_aper.append((-2.5 * (math.log10(float(x[26])))))
                fwhm.append(x[7])
                sys.stdout.write(line)

        plt.plot(flux_aper, fwhm, 'bo')
        del flux_aper[:], fwhm[:]

        plt.savefig((sept + 'Sept flux vs. fwhm.png'), bbox_inches='tight')
        if writeup == True:
            plt.show()
        plt.clf()
        plt.close()


# calculate distance between every source. if two objects are closer than 3.3 arcsec, throw them both out
def removeNearObjs(month):
    global separation
    octxylist = []
    if month == 'oct':
        for num in range(len(oct_extracted)):
            # get xy for extracted sources
            thefile = oct_extracted[num]
            with open(thefile, "r") as f:
                lines = f.readlines()
                for x in lines:
                    line = x.split(' ')
                    if not x.startswith("N"):
                        octxylist.append((float(line[1]), float(line[2])))

            neighbors = []
            for i in range(len(octxylist) - 1):
                a = np.array([distance(*pair) for pair in zip(repeat(octxylist[i]), octxylist[(i + 1):])])
                b = np.where(a < separation)[0]
                if len(b) > 0:
                    neighbors.append(i + 1)
                    for z in b:
                        neighbors.append(z + i + 2)

            tempx, tempy = [], []
            for line in fileinput.input(octplot[num], inplace=1):
                x = line.split(' ')
                if int(x[0]) not in neighbors:
                    tempx.append(x[1])
                    tempy.append(x[2])
                    stars.append(int(x[0]))
                    sys.stdout.write(line)
            liststars.append(stars[:])

            del neighbors[:], octxylist[:], stars[:]

            with open((userdata + 'c%s.xml' % (num + 1)), 'w') as f:
                for y in range(len(tempx)):
                    f.write('circle(%s,%s,%s)' % (tempx[y], tempy[y], separation - 15) + '\n')
    elif month == 'sept':
        # get xy for extracted sources
        thefile = septextracted[0]
        septxylist = []
        with open(thefile, "r") as f:
            lines = f.readlines()
            for x in lines:
                line = x.split(' ')
                if not x.startswith("N"):
                    septxylist.append((float(line[1]), float(line[2])))

        neighbors = []
        for i in range(len(septxylist) - 1):
            a = np.array([distance(*pair) for pair in zip(repeat(septxylist[i]), septxylist[(i + 1):])])
            b = np.where(a < separation)[0]
            if len(b) > 0:
                neighbors.append(i + 1)
                for z in b:
                    neighbors.append(z + i + 2)

        tempx, tempy = [], []
        for line in fileinput.input(septclean[0], inplace=1):
            x = line.split(' ')
            if int(x[0]) not in neighbors:
                tempx.append(x[1])
                tempy.append(x[2])
                septstars.append(int(x[0]))
                sys.stdout.write(line)
        del neighbors[:], septxylist[:]

        with open((userdata + 'septstars.xml'), 'w') as f:
            for y in range(len(tempx)):
                f.write('circle(%s,%s,%s)' % (tempx[y], tempy[y], separation - 15) + '\n')


# for sources which fit stellar sequence and have no neighbors, plot curve of growth plot.
# take worst quality of the four chips as the radius aperture to use to calculate zeropoints
def findAperture(month):
    global convergeap, convergesept
    if month == 'oct':
        med, converge, flux_list = [], [], []
        for num in range(len(octplot)):
            plt.figure(num)
            with open(octplot[num], 'r') as f:
                lines = f.readlines()
                for x in lines:
                    line = x.split(' ')
                    for i in range(7, (len(aperture_list_pix) + 6)):
                        flux_list.append((float(line[i])) / (float(line[len(aperture_list_asec) + 5])))
                        # if ((float(line[i])) / (float(line[len(aperture_list_asec)+5]))) > 1:
                        # print line[0]
                    plt.plot(aperture_list_asec_radius, flux_list, 'ro')
                    del flux_list[:]
            with open(octplot[num], 'r') as f:
                lines = f.readlines()
                for z in range(7, (len(aperture_list_pix) + 6)):
                    for x in lines:
                        line = x.split(' ')
                        flux_list.append((float(line[z])) / (float(line[len(aperture_list_asec) + 5])))
                    med.append(float(median(flux_list)))
                    del flux_list[:]
            m = np.array(med)
            converge.append((np.where(m > .99)[0])[0])
            plt.plot(aperture_list_asec_radius, med, 'k')
            del med[:]
            plt.title('c%s F(<R)/F(R) vs. R' % (num + 1))
            axes = plt.gca()
            axes.set_ylim([0, 1.2])
            plt.xlabel('Aperture radius (asec)')
            plt.ylabel(
                'Flux(<%s asec)/Flux(%s asec)' % ((aperture_list_asec_radius[len(aperture_list_asec_radius) - 1]),
                                                  (aperture_list_asec_radius[len(aperture_list_asec_radius) - 1])))
            plotname = 'c%s_apertureplot.png' % (num + 1)
            if writeup:
                plt.show()
            plt.savefig(plotname, bbox_inches='tight')
            plt.clf()
            plt.close()
        convergeap = (max(converge) + 7)
        worst = aperture_list_asec_radius[max(converge)]
        print 'Use %s arcsec radius to calculate zeropoints for October. \n' % (worst)
    elif month == 'sept':
        med = []
        flux_list = []
        plt.figure(75)
        thefile = septclean[0]
        with open(thefile, 'r') as f:
            lines = f.readlines()
            for x in lines:
                line = x.split(' ')
                for i in range(9, (len(aperture_list_pix) + 8)):
                    flux_list.append((float(line[i])) / (float(line[len(aperture_list_asec) + 7])))
                plt.plot(aperture_list_asec_radius, flux_list, 'ro')
                del flux_list[:]
        with open(thefile, 'r') as f:
            lines = f.readlines()
            for z in range(9, (len(aperture_list_pix) + 8)):
                for x in lines:
                    line = x.split(' ')
                    flux_list.append((float(line[z])) / (float(line[len(aperture_list_asec) + 7])))
                med.append(float(median(flux_list)))
                del flux_list[:]
        m = np.array(med)
        convergent = (np.where(m > .99)[0])[0]
        plt.plot(aperture_list_asec_radius, med, 'k')
        del med[:]
        plt.title('F(<R)/F(R) vs. R')
        axes = plt.gca()
        axes.set_ylim([0, 1.2])
        plt.xlabel('Aperture radius (asec)')
        plt.ylabel(
            'Flux(<%s asec)/Flux(%s asec)' % ((aperture_list_asec_radius[len(aperture_list_asec_radius) - 1]),
                                              (aperture_list_asec_radius[len(aperture_list_asec_radius) - 1])))
        plotname = sept + 'Sept_apertureplot.png'
        if writeup:
            plt.show()
        plt.savefig(plotname, bbox_inches='tight')
        plt.clf()
        plt.close()
        convergesept = convergent + 9
        worst = aperture_list_asec_radius[convergent]
        print 'Use %s arcsec radius to calculate zeropoints for September. ' % (worst)


def matchOctSept():
    global liststars, masteroctlines, masterseptlines
    septxlist, septylist, septmatchlinelist, octmatchlinelist, septcatlinelist, octcatlinelist, oct_total, octlineno, \
    septlineno, octxlist, octylist = [], [], [], [], [], [], [], [], [], [], []
    # create list of acceptable sources in september
    with open(septextracted[0], "r") as f:
        lines = f.readlines()
        RAs, Decs = [], []

        for x in lines:
            line = x.split(' ')
            if (not x.startswith("N")) and int(line[3]) == 0 \
                    and uppersept > float(line[1]) > lowersept and uppersept > float(line[2]) > lowersept:
                septxlist.append(line[1])
                septylist.append(line[2])
                septlineno.append('#%s' % line[0])
                RAs.append('#%s' % line[5])
                Decs.append('#%s' % line[6])

    with open((userdata + 'match1.xml'), 'w') as f:
        for y in range(len(septxlist)):
            f.write('circle(%s,%s,%s)' % (septxlist[y], septylist[y], 30) + '\n')

    septxy = ((sept + 'xysept_full.txt'))
    ascii.write([septxlist, septylist, septlineno, RAs, Decs], septxy,
                names=('#x', 'y', 'line no.', 'RA', 'Dec'), overwrite=True)

    for num in range(len(octclean)):
        octfile = octclean[num]

        with open(octfile, "r") as f:
            lines = f.readlines()
            for x in lines:
                line = x.split(' ')
                if int(line[3]) == 0:
                    octxlist.append(line[1])
                    octylist.append(line[2])
                    octlineno.append('#%s' % line[0])

        octxy = oct + 'c%s_xyoct_full.txt' % (num + 1)
        ascii.write([octxlist, octylist, octlineno], octxy, names=('#x', 'y', 'line no.'), overwrite=True)
        del octxlist[:], octylist[:], octlineno[:]
        str_out = 'c%s_oct_sept_match.txt' % (num + 1)
        refpnt = (reference + '%sref.txt' % (num + 1))

        iraf.xyxymatch(input=(sept + 'xysept_full.txt'), reference=octxy, output=str_out, xmag=1.0, ymag=1.0,
                       refpoints=refpnt, matching="tolerance",
                       tolerance=3, verbose="no")
        with open(str_out, "r") as f:
            lines = f.readlines()
            for x in lines:
                x = x.strip()
                if not x or x.startswith('#'):
                    continue
                septmatchlinelist.append(int(x.split()[5]) - 1)
                octmatchlinelist.append(int(x.split()[4]) - 1)

        if len(septmatchlinelist) != len(set(septmatchlinelist)) or len(octmatchlinelist) != len(set(octmatchlinelist)):
            print "Error: tweak xyxymatch parameters so that a coordinate in one set of data isn't matched to " \
                  "multiple coordinates in the other"
            print num
            sys.exit()

        with open((sept + 'xysept_full.txt'), 'r') as fs:
            tempsept, tempoct = np.array(septmatchlinelist), np.array(octmatchlinelist)
            with open(octxy, 'r') as f:
                lines = fs.readlines()
                linesoct = f.readlines()
                i = 1
                for x in lines:
                    line = x.split('#')
                    if x.startswith("#"):
                        continue
                    elif i in septmatchlinelist:
                        septcatlinelist.append(int(line[1]))
                        masterseptlines[num].append(int(line[1]))
                        RA.append(float(line[2]))
                        Dec.append(float(line[3]))
                        p = linesoct[octmatchlinelist[np.where(tempsept == i)[0][0]]].split('#')
                        octcatlinelist.append(int(p[1]))
                        masteroctlines[num].append(int(p[1]))
                    i += 1

        temptotal = []
        with open(octfile, 'r') as f:
            tempoct = np.array(octcatlinelist)
            lines = f.readlines()
            for x in lines:
                z = int(x.split(' ')[0])
                if x.startswith("#"):
                    continue
                elif z in octcatlinelist:
                    a = np.where(tempoct == z)[0][0]
                    temptotal.append('%s %s %s %s\n' % (x.rstrip('\n'), RA[a], Dec[a], (num + 1)))

        for item in range(len(temptotal)):
            if int(temptotal[item].split(' ')[0]) in liststars[num]:
                oct_total.append(temptotal[item])
        del septmatchlinelist[:], octmatchlinelist[:], septcatlinelist[:], octcatlinelist[:], RA[:], Dec[:], temptotal[
                                                                                                             :]

    with open('oct_total.txt', 'w') as f:
        for item in oct_total:
            f.write('%s' % item)


def twomassmatch():
    twomassRA, twomassDec, mag, lineno, octRA, octDec, octlineno, octflux, chipno = [], [], [], [], [], [], [], [], []
    global bandline, indexjhk
    with open('2mass.txt', 'r') as f:
        lines = f.readlines()
        i = 1
        for x in lines:
            line = x.split()
            if x.startswith('#') == False and x.split()[21][indexjhk] == '0' \
                    and x.split()[23] == '0' and x.split()[24] == '0' \
                    and (x.split()[18][indexjhk] == 'A' or x.split()[18][indexjhk] == 'B'):
                twomassRA.append(float(line[0]) * 1000)
                twomassDec.append(float(line[1]) * 1000)
                mag.append("#" + line[bandline])
                lineno.append("#" + str(i))
                i += 1

    ascii.write([twomassRA, twomassDec, mag, lineno], '2mass_RADec.txt', names=('#RA', 'Dec', 'Mag', 'Line number'),
                overwrite=True)
    with open('oct_total.txt', 'r') as f:
        lines = f.readlines()
        for x in lines:
            line = x.split(' ')
            octRA.append(float(line[len(line) - 3]) * 1000)
            octDec.append(float(line[len(line) - 2]) * 1000)
            octlineno.append(line[0])
            octflux.append(line[convergeap])
            chipno.append(line[len(line) - 1])

    ascii.write([octRA, octDec, octlineno, chipno, octflux], 'Oct_RADec.txt', names=('#RA', 'Dec', 'line', 'chip',
                                                                                     'flux'), overwrite=True)
    iraf.xyxymatch(input='Oct_RADec.txt', reference='2mass_RADec.txt', output="2mass_Oct_matched.txt", tolerance=3,
                   verbose="no", xref=0, yref=0, xin=0, yin=0)


def calculateOctZPs():
    global octzp, octzperr
    twommatch, octmatch, Mag, flux, octlineno, chipno = [], [], [], [], [], []
    with open("2mass_Oct_matched.txt", 'r') as f:
        lines = f.readlines()
        for x in lines:
            x = x.strip()
            line = x.split()
            if not x or x.startswith('#'):
                continue
            twommatch.append((int(line[4]) - 1))
            octmatch.append((int(line[5]) - 1))

    if len(twommatch) != len(set(twommatch)) or len(octmatch) != len(set(octmatch)):
        print "Error: tweak xyxymatch parameters so that a coordinate in one set of data isn't matched to " \
              "multiple coordinates in the other."
        sys.exit()

    with open('2mass_RADec.txt', 'r') as f:
        temptwo_m = np.array(twommatch)
        with open('Oct_RADec.txt', 'r') as fs:
            lines = f.readlines()
            linesoct = fs.readlines()
            i = 1
            for x in lines:
                z = x.split('#')
                if x.startswith("#"):
                    continue
                elif i in twommatch:
                    b = octmatch[np.where(temptwo_m == i)[0][0]]
                    Mag.append(z[1])
                    p = linesoct[b].split(' ')
                    flux.append(p[4].replace('\n', ''))
                    octlineno.append(p[2])
                    chipno.append(p[3])
                i += 1

    zp = []
    for num in range(0, len(Mag)):
        if float(flux[num]) < 1:
            continue
        zp.append(float(Mag[num]) + float(2.5 * (math.log10(float(flux[num])))))
        # if float(Mag[num]) + float(2.5 * (math.log10(float(flux[num])))) < 26.4:
        # print flux[num]

    octzp, octzperr = median(zp), st.biweight_midvariance(zp)
    print "%s zeropoints calculated in October" % (len(zp))
    plt.figure(30), plt.hist(zp, bins=15), plt.title("October Zeropoints Histogram")
    plt.text(min(zp), len(zp) / 5, "Med: %s \n BWMV: %s" % (octzp, octzperr))
    axes = plt.gca()
    # axes.set_xlim([25, 28])
    if writeup == True:
        plt.show()
    plt.savefig('oct_zps.png', bbox_inches='tight'), plt.close()


def calculateoctmags():
    global convergeap
    for file in range(len(oct_extracted)):
        temp = []
        templine = []
        title = ''
        filename = oct_extracted[file]
        with open(filename, 'r+') as f:
            lines = f.readlines()
            for x in lines:
                line = x.split(' ')
                if x.startswith("N"):
                    x = x.strip('\n')
                    title = '%s MAG\n' % (x)
                    continue
                elif float(line[convergeap]) > 0 and int(line[3]) == 0:
                    temp.append((-2.5 * math.log10(float(line[convergeap]))) + octzp)
                    templine.append(x.replace('\n', ''))
            f.truncate()
        with open(filename, 'w') as f:
            f.write(title)
            for x in range(len(temp)):
                f.write('%s %s\n' % (templine[x], temp[x]))


def calculateSeptZPs():
    zp = []
    for num in range(len(oct_extracted)):
        thefile = oct_extracted[num]
        mags, fluxes = [], []
        something = np.array(masteroctlines[num])
        with open(thefile, 'r') as f:
            with open(septextracted[0], 'r') as fs:
                lines = f.readlines()
                linessept = fs.readlines()
                for x in lines:
                    if x.startswith("N"):
                        continue
                    elif int(x.split(' ')[0]) in masteroctlines[num]:
                        mags.append(float(x.split(' ')[len(x.split(' ')) - 1]))
                        b = masterseptlines[num][np.where(something == int(x.split(' ')[0]))[0][0]]
                        if int(linessept[b].split(' ')[0]) in septstars:
                            fluxes.append(float(linessept[b].split(' ')[convergesept]))
                        else:
                            mags.pop()
        for z in range(len(mags)):
            zp.append(mags[z] + 2.5 * math.log10(fluxes[z]))
    plt.figure(69)
    med = median(zp)
    bwmv = st.biweight_midvariance(zp)
    plt.hist(zp, bins=15, color='g'), plt.title("September Zeropoints Histogram")
    plt.text(min(zp), 7, "Med: %s\n BMWV: %s" % (med, bwmv))
    if writeup == True:
        plt.show()
    plt.savefig('septzeropoints.png'), plt.clf(), plt.close()
    pm = u'\u00B1'
    print "Zeropoint for September: %s %s %s" % (med, pm, bwmv)


def extractionComment():
    print "First, Sextractor is run on September data. The pixel scale is given as 0.159 asec/pixel, and the FWHM is " \
          "measured to be 4.9 pixels. Therefore the SEEING_FWHM is calculated to be 0.78 asec. This value " \
          "is used to configure Sextractor. The parameters retrieved are x and y " \
          "coordinates, flags, star class, RA & Dec, FWHM_Image, and fluxes from apertures of size 0.4'' to 6.6 '', " \
          "in diameter, counting by 0.2'' \n \n The same is done on October data, though RAs and Decs are excluded, " \
          "and FWHM is taken to be 2.8, or 0.4452 arcsec. Catalogs and segmentation files " \
          "are written for both October and September data.\n \n"


def plotComment():
    print "Now sources are selected from the October catalogues which throw no flags & lie on the stellar sequence. " \
          "In order to select these, for each chip, -2.5log_10(flux) is plotted against FWHM. " \
          "Take a range of flux values for which the graph flattens out, and does not show signs of saturation " \
          "Then accept sources which fall within a given number of BWMV sigmas of the median of sources within" \
          " this range. Store these unsaturated stars in a list to be referenced later.\n \n"


def apertureComment():
    print "Distances are calculated between all sources in each image. Sources which are closer than 3.3 arcsec " \
          "to each other are thrown out of the set of stars selected in the previous step. " \
          "Now a text file remains with sources which are likely stars and have " \
          "no nearby companions. These stars are used to make curve of growth plots for each chip. The median of " \
          "each spread is taken, and the aperture at which the median converges to at least 99% of the total flux is " \
          "stored. The worst image quality of the four chips is taken to be the aperture which will be used to " \
          "calculate zeropoints later on.\n \n"


def matchoctseptcomment():
    print "Now all sources in October which don't throw flags and are a reasonable distance from the edge of the " \
          "image are matched with sources in september which meet the same criteria, using xyxymatch. " \
          "In order to match, reference points must be provided due to the" \
          " extreme transformation/rotation. Once they are matched, RA's and Decs are ripped from a source in " \
          "September data and appended to the corresponding source in October. All sources in October which now" \
          "have RA's and Decs now are compiled into one text file.\n \n"


def twomassmatchcomment():
    print "xyxymatch is used to match 2mass sources which satisfy \n :ph_qual == A || B, \n cc_flg, gal_contam," \
          " mp_flg == 0: \n \n with sources in the master October list, which now have RAs and Decs. " \
          "A magnitude is grabbed from a matched 2mass star, and a flux at the previously determined convergent " \
          "aperture is grabbed from the corresponding october source, and a zeropoint is calculated." \
          " This is done for all stars in October data which have a 2mass match."


def OctZPcomment():
    print "Now a distribution of zeropoints for October data remains, and the median of that distribution is taken" \
          " to be a good estimate. This is then used to calculate magnitudes for every source in October. "


def SeptZPcomment():
    print "Now we run all the same tests on September data that we did on October data, finding probable stars, " \
          "removing neighbors to prevent contamination, and finding the convergent aperture size. Now, referring " \
          "to matched stars between Oct and Sept, and given that we now know magnitudes for Oct stars, and the" \
          " convergent aperture size for Sept, we can calculate zeropoints for September stars. A histogram of these" \
          " zeropoints is displayed, and the result is printed. The estimate for the Sept zeropoint is the median" \
          " plus/minus the biweight midvariance of the spread of zeropoints."


def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
        return None
    if len(lst) % 2 == 1:
        return lst[((len(lst) + 1) / 2) - 1]
    else:
        return float(sum(lst[(len(lst) / 2) - 1:(len(lst) / 2) + 1])) / 2.0


def distance(p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    return math.hypot(x2 - x1, y2 - y1)


def run():
    extractSept()
    extractOct()
    if writeup:
        extractionComment()
    if writeup:
        plotComment()
    makeplot('oct')
    if writeup:
        apertureComment()
    removeNearObjs('oct')
    findAperture('oct')

    matchOctSept()
    if writeup:
        matchoctseptcomment()
    if writeup:
        twomassmatchcomment()
    twomassmatch()
    if writeup:
        OctZPcomment()
    calculateOctZPs()
    calculateoctmags()
    if writeup:
        SeptZPcomment()
    makeplot('sept')
    removeNearObjs('sept')
    findAperture('sept')
    calculateSeptZPs()


run()
