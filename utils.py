# -*- coding: utf-8 -*-
"""
Utilities for interacting with the gPhoton pipeline and data.

Created on Thu May 15 14:16:29 2014

@author: Parke
"""

from astropy.table import Table, Column
import astropy.io.fits as fits
import numpy as np
from scipy.stats import mode
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pdb, pointers

def read_curve(filename,cull_bad_data=False):
    
    names = ['t0', 't1', 'ap_radius', 'exptime', 'cps', 'cps_err', 'flux', 
             'flux_err', 'mag', 'mag_err', 'bkgnd_radius0', 'bkgnd_radius1', 
             'bkgnd_cps', 'response', 'counts', 'apcorrect1', 'apcorrect2']
    curve_tbl = Table.read(filename,format='ascii',names=names)
    if not curve_tbl.masked: curve_tbl = Table(curve_tbl, masked=True)
    
    gaps = curve_tbl['t0'][1:] - curve_tbl['t0'][:-1]
    dt = guess_tstep(curve_tbl)
    
    exposure_start = np.concatenate([np.array([False]),(gaps > dt)])
    exposure = exposure_start.cumsum()
    exp_column = Column(data=exposure, name='exposure', dtype=int)
    curve_tbl.add_column(exp_column, index=0)
    
#    pdb.set_trace()
    if cull_bad_data:
        remove_line = [curve_tbl['cps'].mask]
        remove_line.append(np.logical_not(curve_tbl['cps'] >= 0.0))
        for exp in np.arange(exposure[-1]+1):
            total_counts = sum(curve_tbl['counts'][exposure == exp])
            if total_counts == 0:
                remove_line.append(exposure == exp)
        remove_line = reduce(np.logical_or, remove_line)
        keep_line = np.logical_not(remove_line)
        curve_tbl = curve_tbl[keep_line]
    
    #augment errors of pts with zero counts
    zero_counts = (curve_tbl['counts'] == 0)
    curve_tbl['cps'][zero_counts] = 1.0/curve_tbl['exptime'][zero_counts]
    
    isBad = (len(curve_tbl) == 0) or (sum(curve_tbl['cps'].mask) > 0)
    return curve_tbl, isBad
    
def table_to_csv(table, filename):
    table.write(filename, format='ascii', delimiter=',', fill_values=('--',''))
    
def exposure_endpts(curve_tbl):
    exposure = curve_tbl['exposure']
    endpts_boolean = exposure[1:] - exposure[:-1]
    endpts_indices = np.nonzero(endpts_boolean)
    return endpts_indices

def exposure_times(curve_tbl):
    exp = np.unique(curve_tbl['exposure'])
    steps = curve_tbl['t1'] - curve_tbl['t0']
    def sumgroup(i):
        in_exp = (curve_tbl['exposure'] == i)
        return sum(steps[in_exp])
    return map(sumgroup,exp)
    
def guess_tstep(curve_tbl):
    gaps = curve_tbl['t0'][1:] - curve_tbl['t0'][:-1]
    return mode(gaps)[0][0]
    
def collapse_tvec(tvec,maxgap,newgap=None):
    if not newgap: newgap = maxgap
    tvec = tvec
    gaps = np.concatenate(([0.0],tvec[1:] - tvec[:-1]))
    gapindex = np.nonzero(gaps > maxgap)[0]
    gaplen = gaps[gaps > maxgap]
    gaps[gaps > maxgap] = newgap
    tvec_new = gaps.cumsum()
    midgap = (tvec_new[gapindex-1] + tvec_new[gapindex])/2.0
    return tvec_new, midgap, gaplen
    
def quickplot(curve,yvals='cps',spaced=False):
    if curve.__class__ == str:
        curve = read_curve(curve,cull_bad_data=True)
    
    if spaced:
        t = curve['t0'].data.data
        t = t - t[0]
        xlbl = 'Time from Data Start [s]'
    else:
        tstep = guess_tstep(curve)
        t = collapse_tvec(curve['t0'].data.data, 3*tstep) + tstep/2.
        xlbl = 'Arbitrary Time (Exposure Gaps Shortened) [s]'
    
    fig = plt.figure()
    plt.errorbar(x=t, y=curve[yvals], 
             yerr=curve['{}_err'.format(yvals)], fmt='.')
    plt.ylabel(yvals)
    plt.xlabel(xlbl)
    
    if not spaced:    
        end_indices = exposure_endpts(curve)
        end_times = t[end_indices] + tstep/2
        draw_line = lambda end_time: plt.axvline(end_time,ls=':',color='k')
        map(draw_line,end_times)
    
    plt.show()
    return fig
    
def thin(cat):
    good_kids = np.load(pointers.good_kids_file)
    keep = np.in1d(cat['Kepler ID'],good_kids)
    return cat[keep]
    
def Fbol_SDSSg(Teff_vec, gmag_vec):
    """Computes the bolometric flux in somewhat arbitrary, magnitude based
    units for comparison with other fluxes. 
    
    Flux is computed as 10**(-2.5*mag), such that it can be compared with
    other magnitudes.
    """
    Teff_vec = np.float64(Teff_vec)
    
    hdulist = fits.open('../Scratchwork/SDSS_filter_curves.fits')
    data = hdulist[2].data
    wvln = data['wavelength'] #ang
    resp = data['respt']
    wvln = wvln/1e4 #microns
    minw, maxw = wvln[0], wvln[-1]
    thru = interp1d(wvln, resp, bounds_error=False, fill_value=0.0)
    
    hc_k = 14388.8 #microns*K
    hc4pi = 4.99 #10^18 J*micron
    rad = 7.57e-16 #10^18 J micron-3 K-4
    def g2bol_ratio(Teff):
        plnk = lambda w: hc4pi / w**5 / (np.exp(hc_k/Teff/w) - 1.0) #J micron-4
        total = rad*Teff**4 #J micron-3
        part = quad(lambda w: plnk(w)*thru(w), minw, maxw)[0]
        return part/total
        #I checked and integrating plnk from 0 to inf gives the same as total
    
    g2bol_vec = np.array(map(g2bol_ratio, Teff_vec))
    FgF0 = 10**(-2.5*gmag_vec) #ratio of flux to a reference flux of mag 0
    FbolF0 = FgF0/g2bol_vec
    return FbolF0
    
def find_loners(radec, radec_all, radius):
    """Returns the indices of only the stars with no neighbors within radius.
    
    radec contains the coordinates of the stars of interest as an Nx2 numpy 
    array (ra, dec), and radec_all of every source that might contaminate the
    stellar flux.
    """
    
    loners = np.ones(len(radec))
    for i,(ra,dec) in enumerate(radec):
        dra = abs(radec_all[:,0] - ra)
        ddec = abs(radec_all[:,1] - dec)
        keep = np.logical_and(dra < radius, ddec < radius)
        r = np.sqrt((dra[keep]**2 + ddec[keep]**2))
        r = r[r != 0]
        if any(r < radius):
            loners[i] = False
            
    return loners