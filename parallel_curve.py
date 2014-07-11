# -*- coding: utf-8 -*-
"""
Created on Fri May 16 12:58:49 2014

@author: Parke
"""
import pointers
import curvetools as gcurve
import dbasetools as gdb
from multiprocessing import Pool
from optparse import OptionParser
from astropy.table import Table

parser = OptionParser()
parser.add_option("-c", "--cat", "--catalog", action="store", type="string", 
                  dest="cat_file")
parser.add_option('-t', '--tstep', '--dt', action='store', type='float', 
                  dest='tstep')
parser.add_option('-b', '--band', action='store', type='string', 
                  dest='band')
parser.add_option('-p', '--processes', action='store', type='int', dest='Np')
parser.add_option("-a", "--aperture", action="store", type="float", dest="ap_radius")
parser.add_option("-i", "--inner", action="store", type="float", dest="annulus1")
parser.add_option("-o", "--outer", action="store", type="float", dest="annulus2")
(options, args) = parser.parse_args()

cat = Table.read(options.cat_file, format='ascii')
bk_annulus = [options.annulus1,options.annulus2]

def generate_curve(star):
    ra, dec, kid = star
    ofile = '{}{}.{}.csv'.format(pointers.curve_folder, kid, options.band)
    tranges = gdb.fGetTimeRanges(options.band, [ra,dec])
    gcurve.write_curve(options.band, [ra,dec], tranges, options.ap_radius, 
                       outfile=ofile, annulus=bk_annulus, stepsz=options.tstep, 
                       userr=True, usehrbg=True, 
                       calpath=pointers.calibration_folder)

if __name__ == '__main__':
    pool = Pool(processes=options.Np)
    stars = [[row['RA'],row['Dec'],row['Kepler ID']] for row in cat]
    files = pool.map(generate_curve, stars)
#    files = map(generate_curve, cat)
    print files