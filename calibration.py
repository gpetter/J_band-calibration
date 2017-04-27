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
userdata = cur_dir + '/Userdata/'
sex = cur_dir + '/Sex/'



# list of october fits filenames
oct_list = []

# contains october Sextractor catalogues
oct_extracted = []

octclean = []
septlist = []
septextracted = []
septclean = []
masteroctlines = [[], [], [], []]
masterseptlines = [[], [], [], []]

x_origin, y_origin = 0, 0

# lists of xy coordinates used for xyxymatch
octxlist, octylist, linelist, octplot, septsectionslist, liststars, stars, RA, Dec = [], [], [], [], [], [], [], [], []

octzp, octzperr = 0, 0

# tune these values as necessary, one for each chip
seeing_fwhm = .4452
aperture_list_asec = [0, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4,
                      4.6, 4.8, 5, 5.2, 5.4, 5.6, 5.8, 6, 6.2, 6.4, 6.6]
aperture_list_pix = [i*6.28930817 for i in aperture_list_asec]
aperture_list_asec_radius = [float(i)/2 for i in aperture_list_asec]
aperture_list_asec_radius.pop(0)
convergeap, convergesept = 0, 0

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

#run Sextractor on september data and each of the four October images.
#make copies to edit and truncate
def extractSept():
    sew = sewpy.SEW(
        params=["NUMBER", "X_IMAGE", "Y_IMAGE", "FLAGS", "CLASS_STAR", "ALPHA_SKY", "DELTA_SKY", "FWHM_IMAGE",
                "FLUX_APER(%s)" % (len(aperture_list_pix))], config={"MAG_ZEROPOINT": 29.1653, "PIXEL_SCALE": 0.159,
                                    "SEEING_FWHM": .7791, "FILTER_NAME": (sex + "gauss_5.0_9x9.conv"),
                                    "CHECKIMAGE_TYPE": 'SEGMENTATION', "PHOT_APERTURES": (aperture_list_pix),
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
                        config={"MAG_ZEROPOINT": 26.6, "PIXEL_SCALE": 0.159, "SEEING_FWHM": seeing_fwhm,
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
                130 < float(x[1]) < 2250) and (130 < float(x[2]) < 2250) and (float(x[24]) > 0):
                temp1.append(x[1])
                temp2.append(x[2])
                sys.stdout.write(line)
        copyfile(cleanname, plotname)
        with open((userdata + '%soct.xml') % (num+1), 'w') as f:
            for y in range(len(temp1)):
                f.write('circle(%s,%s,%s)' % (temp1[y], temp2[y], 30) + '\n')\

#tune values to permit what can be counted as a star
sigmamultiple = [2.5, 2, 2, 3]

#create plot of all sources' -2.5log(flux_4asec) vs FWHM. Distinguish those sources which lie on a horizontal
# branch of the diagram, indicating star-like
def makeplot(month):
    flux_aper, fwhm = [], []

    global liststars
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
                if (not line.startswith('N')) and (-13.5 < (-2.5 * (math.log10(float(x[24])))) < -7) and float(x[5]) < 3.2:
                    flux_aper.append((-2.5 * (math.log10(float(x[24])))))
                    fwhm.append(float(x[5]))

                    sys.stdout.write(line)
            med = median(fwhm)
            sig = st.biweight_midvariance(fwhm)
            del flux_aper[:], fwhm[:]
            for line in fileinput.input(thefile, inplace=1):
                x = line.split(' ')
                if (not line.startswith('N')) and ((med - sigmamultiple[num]*sig) < float(x[5]) <
                                                       (med + sigmamultiple[num]*sig)):
                    flux_aper.append((-2.5 * (math.log10(float(x[24])))))
                    fwhm.append(x[5])
                    stars.append(int(x[0]))
                    sys.stdout.write(line)
            liststars.append(stars[:])

            plt.plot(flux_aper, fwhm, 'bo')
            del flux_aper[:], fwhm[:], stars[:]

            plt.savefig(('%s_flux vs. fwhm.png' % ('c' + str(num + 1))), bbox_inches='tight')
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
        #axes.set_xlim([-15, -2]), axes.set_ylim([0, 20])
        del flux_aper[:], fwhm[:]
        for line in fileinput.input(thefile, inplace=1):
            x = line.split(' ')
            if (not line.startswith('N')) and (-16 < (-2.5 * (math.log10(float(x[26])))) < -10):
                flux_aper.append((-2.5 * (math.log10(float(x[26])))))
                fwhm.append(float(x[7]))

                sys.stdout.write(line)
        med = median(fwhm)
        sig = st.biweight_midvariance(fwhm)
        del flux_aper[:], fwhm[:]
        for line in fileinput.input(thefile, inplace=1):
            x = line.split(' ')
            if ((med - 1.5*sig) < float(x[7]) < (med + 1.5*sig)):
                flux_aper.append((-2.5 * (math.log10(float(x[26])))))
                fwhm.append(x[7])
                #stars.append(int(x[0]))
                sys.stdout.write(line)
        #liststars.append(stars[:])

        plt.plot(flux_aper, fwhm, 'bo')
        del flux_aper[:], fwhm[:], stars[:]

        plt.savefig((sept + 'Sept flux vs. fwhm.png'), bbox_inches='tight')
        plt.show()
        plt.clf()
        plt.close()


def removeNearObjs(month):
    if month == 'oct':
        for num in range(len(oct_extracted)):
            # get xy for extracted sources
            thefile = oct_extracted[num]
            with open(thefile, "r") as f:
                lines = f.readlines()
                for x in lines:
                    line = x.split(' ')
                    if not x.startswith("N"):
                        octxlist.append(line[1])
                        octylist.append(line[2])

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
                x = line.split(' ')
                if x[0] in templineno:
                    tempx.append(x[1])
                    tempy.append(x[2])
                    sys.stdout.write(line)
            with open((userdata + 'c%s.xml' % (num + 1)), 'w') as f:
                for y in range(len(tempx)):
                    f.write('circle(%s,%s,%s)' % (tempx[y], tempy[y], 45) + '\n')
    elif month == 'sept':
        # get xy for extracted sources
        thefile = septextracted[0]
        with open(thefile, "r") as f:
            lines = f.readlines()
            for x in lines:
                line = x.split(' ')
                if not x.startswith("N"):
                    octxlist.append(line[1])
                    octylist.append(line[2])

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
            x = line.split(' ')
            if x[0] in templineno:
                tempx.append(x[1])
                tempy.append(x[2])
                sys.stdout.write(line)
        with open(('c%s.xml' % (num + 1)), 'w') as f:
            for y in range(len(tempx)):
                f.write('circle(%s,%s,%s)' % (tempx[y], tempy[y], 45) + '\n')



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
                    for i in range(7, (len(aperture_list_pix)+6)):
                        flux_list.append((float(line[i])) / (float(line[len(aperture_list_asec)+5])))
                        #if ((float(line[i])) / (float(line[len(aperture_list_asec)+5]))) > 1:
                            #print line[0]
                    plt.plot(aperture_list_asec_radius, flux_list, 'ro')
                    del flux_list[:]
            with open(octplot[num], 'r') as f:
                lines = f.readlines()
                for z in range(7, (len(aperture_list_pix)+6)):
                    for x in lines:
                        line = x.split(' ')
                        flux_list.append((float(line[z])) / (float(line[len(aperture_list_asec)+5])))
                    med.append(float(median(flux_list)))
                    del flux_list[:]
            m = np.array(med)
            converge.append((np.where(m > .99)[0])[0])
            plt.plot(aperture_list_asec_radius, med, 'k')
            del med[:]
            plt.title('F(<R)/F(R) vs. R')
            axes = plt.gca()
            axes.set_ylim([0, 1.2])
            plt.xlabel('Aperture radius (asec)')
            plt.ylabel('Flux(<%s asec)/Flux(%s asec)' % ((aperture_list_asec_radius[len(aperture_list_asec_radius)-1]),
                                                         (aperture_list_asec_radius[len(aperture_list_asec_radius)-1])))
            plotname = 'c%s_apertureplot.png' % (num + 1)
            plt.savefig(plotname, bbox_inches='tight')
            plt.clf()
            plt.close()
        convergeap = (max(converge)+7)
        worst = aperture_list_asec_radius[max(converge)]
        print 'Use %s arcsec radius to calculate zeropoints for October. ' % (worst)
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
        plt.show()
        plt.savefig(plotname, bbox_inches='tight')
        plt.clf()
        plt.close()
        convergesept = convergent+9
        worst = aperture_list_asec_radius[convergent]
        print 'Use %s arcsec radius to calculate zeropoints for September. ' % (worst)








def matchOctSept():
    global liststars, masteroctlines, masterseptlines
    septxlist, septylist, septmatchlinelist, octmatchlinelist, septcatlinelist, octcatlinelist, oct_total, octlineno,\
    septlineno  = [], [], [], [], [], [], [], [], []
    #create list of acceptable sources in september
    with open(septextracted[0], "r") as f:
        lines = f.readlines()
        RAs, Decs = [], []

        for x in lines:
            line = x.split(' ')
            if (not x.startswith("N")) and int(line[3]) == 0 and float(line[13]) > 4000 \
                    and 4300 > float(line[1]) > 130 and 4300 > float(line[2]) > 130:

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



        iraf.xyxymatch(input=(sept + 'xysept_full.txt'), reference=octxy, output=str_out, xmag=1.0, ymag=1.0,
                       refpoints=('%sref.txt' % (num+1)), matching="tolerance",
                       tolerance=7, verbose="no")
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
                  "multiple coordinates in the other. First try increasing separation value"
            print num
            sys.exit()

        with open((sept+'xysept_full.txt'), 'r') as f:
            lines = f.readlines()
            length = len(lines)
            i = 1
            while i != length:
                x = lines[i]
                line = x.split('#')
                if x.startswith("#"):
                    continue
                elif len(septmatchlinelist) > 0 and i == septmatchlinelist[0]:
                    septcatlinelist.append(int(line[1]))
                    masterseptlines[num].append(int(line[1]))
                    RA.append(float(line[2]))
                    Dec.append(float(line[3]))
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
                    masteroctlines[num].append(int(x.split('#')[1]))
                    octmatchlinelist.pop(0)
                    i = 1
                    continue
                i += 1
        temptotal = []
        with open(octfile, 'r') as f:
            lines = f.readlines()
            length = len(lines)
            k = 0
            i = 0
            while (i) != length:
                x = lines[i]

                if x.startswith("#"):
                    continue
                elif ((len(octcatlinelist) > 0) and (int(x.split(' ')[0]) == octcatlinelist[0])):
                    temptotal.append('%s %s %s %s\n' % (x.rstrip('\n'), RA[k], Dec[k], (num + 1)))
                    octcatlinelist.pop(0)
                    i = 0
                    k += 1
                    continue
                i += 1

        for item in range(len(temptotal)):
            if int(temptotal[item].split(' ')[0]) in liststars[num]:
                oct_total.append(temptotal[item])
        del septmatchlinelist[:], octmatchlinelist[:], septcatlinelist[:], octcatlinelist[:], RA[:], Dec[:], temptotal[:]

    with open('oct_total.txt', 'w') as f:
        for item in oct_total:
            f.write('%s' % item)


def twomassmatch():
    twomassRA, twomassDec, mag, lineno, octRA, octDec, octlineno, octflux, chipno = [], [], [], [], [], [], [], [], []
    with open('2mass.txt', 'r') as f:
        lines = f.readlines()
        i = 1
        for x in lines:
            line = x.split()
            if x.startswith('#') == False:
                twomassRA.append(float(line[0])*1000)
                twomassDec.append(float(line[1])*1000)
                mag.append("#" + line[6])
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
            octflux.append(line[convergeap])
            chipno.append(line[len(line)-1])

    ascii.write([octRA, octDec, octlineno, chipno, octflux], 'Oct_RADec.txt', names=('#RA', 'Dec', 'line', 'chip',
                                                                                     'flux'), overwrite=True)
    iraf.xyxymatch(input='Oct_RADec.txt', reference='2mass_RADec.txt', output="2mass_Oct_matched.txt", tolerance=3,
                   verbose="no", xref=0, yref=0, xin=0, yin=0)

def calculateZP():
    global octzp, octzperr
    twommatch, octmatch, Mag, flux, octlineno, chipno = [], [], [], [], [], []
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
        print "Error: tweak xyxymatch parameters so that a coordinate in one set of data isn't matched to " \
              "multiple coordinates in the other."
        sys.exit()

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
            line = x.split(' ')
            if x.startswith("#"):
                continue
            elif len(octmatch) > 0 and i == octmatch[0]:
                flux.append(line[4].replace('\n', ''))
                octlineno.append(line[2])
                chipno.append(line[3])
                octmatch.pop(0)
                i = 1
                continue
            i += 1
    zp = []
    for num in range(0, len(Mag)):
        if float(flux[num]) < 1:
            continue
        zp.append(float(Mag[num]) + float(2.5*(math.log10(float(flux[num])))))
    plt.figure(30)
    octzp = median(zp)
    octzperr = st.biweight_midvariance(zp)
    print "%s zeropoints calculated in October" % (len(zp))
    plt.hist(zp, bins=20)
    plt.title("October Zeropoints Histogram")
    axes = plt.gca()
    #axes.set_xlim([25, 28])
    plt.show()
    plt.savefig('oct_zps.png', bbox_inches='tight')
    plt.close()

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
                    x=x.strip('\n')
                    title = '%s MAG\n' % (x)
                    continue
                elif float(line[convergeap]) > 0 and int(line[3]) == 0:
                    temp.append((-2.5*math.log10(float(line[convergeap])))+octzp)
                    templine.append(x.replace('\n', ''))
            f.truncate()
        with open(filename, 'w') as f:
            f.write(title)
            for x in range(len(temp)):
                f.write('%s %s\n' % (templine[x], temp[x]))

def calculateSeptMags():

    for num in range(len(oct_extracted)):
        thefile = oct_extracted[num]
        zp, mags, fluxes= [], [], []
        something = np.array(masteroctlines[num])
        with open(thefile, 'r') as f:
            with open(septextracted[0], 'r') as fs:
                lines = f.readlines()
                linessept = fs.readlines()
                for x in lines:
                    if x.startswith("N"):
                        continue
                    elif int(x.split(' ')[0]) in masteroctlines[num]:
                        mags.append(float(x.split(' ')[len(x.split(' '))-1]))
                        a = np.where(something == int(x.split(' ')[0]))[0][0]
                        b = masterseptlines[num][a]
                        #print a, b
                        fluxes.append(float(linessept[b].split(' ')[convergesept]))


        #print mags, fluxes

        for z in range(len(mags)):
            zp.append(mags[z]+2.5*math.log10(fluxes[z]))
    plt.figure(69)
    plt.clf()
    plt.hist(zp, bins=20, color='g')
    plt.title("September Zeropoints Histogram")
    plt.show()
    plt.savefig('septzeropoints.png')
    plt.clf(), plt.close()






def extractionComment():
    print "First, Sextractor is run on September data. The pixel scale is given as 0.159 asec/pixel, and the FWHM is " \
          "measured to be 4.9 pixels. Therefore the SEEING_FWHM is calculated to be 0.78 asec. This value along " \
          "with the zeropoint estimate is used to configure Sextractor. The parameters retrieved are x and y" \
          "coordinates, flags, star class, RA & Dec, FWHM_Image, and fluxes from apertures of size 0.4'' to 6.6 '', " \
          "counting by 0.2'' \n \n The same is done on October data, though RAs and Decs are excluded, and" \
          "FWHM is taken to be 2.8, or 0.4452 arcsec. Catalogs " \
          "are written to files for both October and September data.\n"


def plotComment():
    print "Now sources are selected from the October catalogues which throw no flags, and have a CLASS_STAR  " \
          "above a certain value specified. The flux of that source through the 4'' aperture and the FWHM are " \
          "stored. Then for each chip, -2.5log_10(flux) is plotted against FWHM. Looking at this plot, observe " \
          "whether or not your CLASS_STAR value selected objects on the stellar sequence. If not, take the" \
          "median of FHWMs for sources within a range, and accept sources within a number of sigmas to be stars."

def apertureComment():
    print "xyxymatch is run on a catalog against itself with a separation of 6.3 arcsec, to remove sources with" \
          "companions. Sources which are left, are then cut with the criteria established in the previous step, " \
          "fitting well on the stellar sequence. Now a text file remains with sources which are likely stars and have" \
          "no nearby companions. These stars are used to make curve of growth plots for each chip. The median of " \
          "each spread is taken, and the aperture at which the median converges to at least .99 of the total flux is " \
          "stored. The worst image quality of the four chips is taken to be the aperture which will be used to " \
          "calculate zeropoints later on."

def matchoctseptcomment():
    print "Now all sources in october which don't throw flags and are a reasonable distance from the edge of the image" \
          " are matched with sources in september which meet the same criteria, and which are substantively bright" \
          "(given that Sept data is much more deep). In order to match, reference points must be provided due to the" \
          "severe transformation/rotation. Once they are matched, RA's and Decs are ripped from a source in " \
          "September data and appended to the corresponding source in October. All sources in October which now" \
          "have RA's and Decs now are compiled into one text file."
def twomassmatchcomment():
    print "xyxymatch is used to match (downloaded )2mass sources with sources in the master October list." \
          "A magnitude is grabbed from a matched 2mass star, and a flux at the previously determined convergent" \
          "aperture is grabbed from the corresponding october source, and a zeropoint is calculated."
def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
            return None
    if len(lst) %2 == 1:
            return lst[((len(lst)+1)/2)-1]
    else:
            return float(sum(lst[(len(lst)/2)-1:(len(lst)/2)+1]))/2.0

#extractSept()
#extractOct()
#makeplot('oct')
#removeNearObjs('oct')
#findAperture('oct')
#matchOctSept()
#twomassmatch()
#calculateZP()
#calculateoctmags()
#makeplot('sept')
#findAperture('sept')
#calculateSeptMags()