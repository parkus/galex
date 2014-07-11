# -*- coding: utf-8 -*-
"""
Contains functions to generate various plots related to the GALEX SPI search.

Created on Sun Jun 29 16:03:57 2014

@author: Parke
"""

import pdb, utils, pointers, strutils
import plotutils as pu
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from astropy.table import Table

fnt = mpl.rcParams['font.size']

if pointers.mini:
    starcat = utils.thin(Table.read(pointers.starcat_file, format='ascii'))
    propcat = utils.thin(Table.read(pointers.propcat_file, format='ascii'))
    noisecat = utils.thin(Table.read(pointers.xnoisecat_file, format='ascii'))
    koicat = utils.thin(Table.read(pointers.koicat_file, format='ascii'))
    detections = reduce(np.logical_and, [noisecat['flag'].mask, 
                        noisecat['x_noise'] > 0.0, noisecat['x_noise'] < 2.0])
else:
    starcat = Table.read(pointers.starcat_file, format='ascii')
    propcat = Table.read(pointers.propcat_file, format='ascii')
    koicat = Table.read(pointers.koicat_file, format='ascii')

def filter_throughput(ax):
    #prepare
    ax.set_xlabel('Wavelength [$\AA$]')
    ax.set_ylabel('Log(Flux), Throughput')
    ax.set_yticklabels('')
    
    #plot example spectrum
    h = fits.open('../Media/hd20630 UV spectrum.fits')
    wave, flux = h[1].data['wave'][0], h[1].data['flux'][0]
    Nsmooth = 300
    def binup(x):
        end = len(x)//Nsmooth*Nsmooth
        x = x[:end]
        arr = np.reshape(x,[len(x)/Nsmooth, Nsmooth])
        return np.mean(arr, 1)
    wave, flux = map(binup, [wave, flux])
    flux = np.log10(flux)
    flux -= np.nanmin(flux)
    flux = flux/np.nanmax(flux)
    ax.plot(wave, flux, 'k')
    ax.legend(['Ex. Spectrum: HD 20630'], loc='lower right', fontsize=0.9*fnt)
    h.close()
    
    #plot throughput curves
    for band, color in zip(['FUV','NUV'], ['b','r']):
        curve = Table.read('../Media/'+band+' throughput curve.txt', 
                               format='ascii', delimiter=' ', 
                               data_start=1, names=['wav','thru'])
        curve['thru'] = curve['thru']/max(curve['thru'])
        ax.fill(curve['wav'], curve['thru'], color=color, alpha=0.4)
        ax.plot(curve['wav'], curve['thru'], color=color)
        mn = (sum(curve['wav']**2*curve['thru'])/
              sum(curve['wav']*curve['thru']))
        ax.text(mn, 0.55, band, ha='center', color=color)
        
    pu.tight_axis_limits(ax, xory='x')
    plt.draw()

def UVvRa(ax, band):
    #prep
    ax.set_xlabel('$R_{planet}/a$ (R$_\oplus$/AU)')
    ax.set_ylabel('{} - g'.format(band))
    
    #keep only the planets with good radii
    keep = np.logical_and(koicat['r_p']/koicat['r_p_err1'] > 3,
                          koicat['r_p'] < 50.0)
    xcat = koicat[keep]
    
    #compute x
    x = xcat['r_p']/xcat['a']
    
    #compute y
    def makey(kid):
        star = starcat['Kepler ID'] == kid
        if starcat.mask[band][star] or propcat.mask['g mag'][star]:
            return 0.0
        else:
            return starcat[band][star].data.data - propcat['g mag'][star].data.data
    y = np.array(map(makey, xcat['Kepler ID']))
    
    #keep only good y
    keep = (y != 0.0)
    x, y = x[keep], y[keep]
    
    ax.semilogx(x, y, 'k.')
    plt.draw()

def exlightcurve(ax, kid=9468475,tag='9468475',band='FUV',yticks=None):
    curve = utils.read_curve(pointers.curvefile(kid,band), cull_bad_data=True)[0]
    t = (curve['t0'] + curve['t1']).data.data/2.0
    t, gapi, gaplen = utils.collapse_tvec(t, 800.0)
    gapt = (t[gapi] + t[gapi-1])/2.0
    y = curve['cps']
    yerr = curve['cps_err']
    
    ax.errorbar(t, y, yerr, fmt='ko', capsize=0, ms=3)
    ax.set_xlabel('Time (scale only) (s)')
    ax.set_ylabel('Counts/s')
    ylim = ax.get_ylim()
    if yticks: ax.set_yticks(yticks)
    
    ytxt = ylim[0] + 0.8*(ylim[1] - ylim[0])
    for gt,glen in zip(gapt,gaplen):
        lenstr = strutils.pretty_deltatime(glen)
        ax.text(gt, ytxt, lenstr, ha='center', va='top', rotation='vertical', 
             color='gray')
    
    pu.tight_axis_limits(ax, xory='x')
    plt.draw()
    
def VarvUV(ax, band):
    xdata = starcat[band] - propcat['g mag']
    ax.set_ylabel('{} Excess Noise'.format(band))
    goodx = np.logical_not(xdata.mask)
    good = np.logical_and(detections, goodx)
    ax.plot(xdata[good], noisecat['x_noise'][good], 'k.')
    ylim = ax.get_ylim()
    ax.set_ylim([0,ylim[1]])
    ax.set_xlabel('{} - g'.format(band))
    plt.draw()
        
def VarvRa(ax, band):
    #prep
    ax.set_xlabel('$R_{planet}/a$ (R$_\oplus$/AU)')
    ax.set_ylabel('{} Excess Noise'.format(band))
    
    #keep only the planets with good radii and excess noise
    goodx = np.logical_and(koicat['r_p']/koicat['r_p_err1'] > 3,
                          koicat['r_p'] < 50.0)
    keep = np.logical_and(detections, goodx)
    xcat = koicat[keep]
    
    #compute x
    x = xcat['r_p']/xcat['a']
    
    ax.semilogx(x, noisecat['x_noise'][keep], 'k.')
    plt.draw()