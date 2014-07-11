# -*- coding: utf-8 -*-
"""
Created on Fri May 16 12:58:49 2014

@author: Parke
"""
import pointers, utils
from multiprocessing import Pool
from optparse import OptionParser
from astropy.table import Table
import variability as vb
import pdfutils as pu
import pdb

parser = OptionParser()
parser.add_option("-c", "--cat", "--catalog", action="store", type="string", 
                  dest="cat_file")
parser.add_option('-b', '--band', action='store', type='string', dest='band')
parser.add_option('-p', '--processes', action='store', type='int', dest='Np')
(options, args) = parser.parse_args()

cat = Table.read(options.cat_file, format='ascii')

def compute_xnoise(kid):
    curvefile = '{}{}.{}.csv'.format(pointers.curve_folder, 
                                     kid, options.band)
    curve, isBad = utils.read_curve(curvefile,cull_bad_data=True)
    limitflag, x_noise, err0, err1 = '', None, None, None
    if not isBad and len(curve) > 3:
        try:
            #get prob dist func of sigma/mu
            pdf, x_noise = vb.excess_noise_PDF(curve['cps'], curve['cps_err'])
        
            width_appx = 1.0/pdf(x_noise)
            limit683 = pu.upper_limit(pdf, 0.683, normalized=True, x0=0.0, 
                                      xpeak=x_noise, x1guess=x_noise + 0.5*width_appx)
            if pdf(limit683) > pdf(0.0): #if there is a well-defined 68.3% interval
                xlo, xhi = pu.confidence_interval(pdf, x_noise, normalized=True)
                err0 = x_noise - xlo
                err1 = xhi - x_noise
            else:
                limitflag = '<'
                x_noise = pu.upper_limit(pdf, normalized=True, x0=0.0, xpeak=x_noise,
                                         x1guess=x_noise + 2*width_appx)
        except:
            pass
    row = [kid,limitflag,x_noise,err0,err1]
    return row
    

if __name__ == '__main__':
    kids = cat['Kepler ID']
    pool = Pool(processes=options.Np)
    rows = pool.map(compute_xnoise, kids)
#    compute_xnoise(kids[153])
#    rows = []
#    for i, kid in enumerate(kids):
#        print i
#        rows.append(compute_xnoise(kid))
#   rows = map(compute_xnoise, kids)
#    pdb.set_trace()
    xcat = Table(names=['Kepler ID','flag','x_noise','-err','+err'],
                 masked=True, dtype=[float,str,float,float,float])
    for row in rows: xcat.add_row(row, [r == None for r in row])
    xcat.sort(['Kepler ID'])
    utils.table_to_csv(xcat, pointers.xnoisecat_file)