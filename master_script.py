# -*- coding: utf-8 -*-
"""The set of code that contains all the instructions necessary to completely
run the GALEX time-series analysis script and produce the variability catalog
and the numbers, tables, and figures used in the associated paper.

Created on Tue May 06 15:22:06 2014

@author: Parke
"""
import aplpy, utils, os, multiprocessing, pointers, time, imp, pdb
import numpy as np
import my_numpy as mynp
import dbasetools as gdb
import catalogIO as catIO
import variability as vb
import scipy.optimize as opt
import pdfutils as pu
from astropy.table import Table, Column

redo = False #toggle switch for where to restart the analysis

#GLOBAL VARIABLES
ap_radius = 0.0025 # aperture radius for photometry, deg
mag_cut = 18
tstep = 60.0
Teff_sun = 5777 #K
bk_annulus = np.array([1.5,2])*ap_radius
band = 'FUV'
#sample_criteria = ['coordformat=dec','ra_galex=!\\null','FUV<18']
sample_criteria = ['coordformat=dec','ra_galex=!\\null']

#%% GENERATE THE CATALOG OF SOURCES TO BE ANALYZED

if redo:
    
    #get the kep ids for all sources in the kepler catalog with galex coordinates
    #all of these appear in the "gold match" kepler-galex cross match catalog
    #there are issues with the GALEX IDs (too long?), so for now I'm ignoring them
    #there do not seem to be duplicates in the HTML tables
    
    #now get the ra and dec from the MCAT database gphoton uses
    cols = ['kic_kepler_id','ra_GALEX','dec_GALEX','NUV','FUV']
    starcat = catIO.fetch_MASTcatalog_data('kepler/kepler_fov',cols,
                                           sample_criteria)
                                           
    #clean up column labels
    old_cols = ['Kepler_ID','col1','col2','col3','col4']
    new_cols = ['Kepler ID','RA','Dec','NUV','FUV']
    map(starcat.rename_column, old_cols, new_cols)
    starcat.sort('Kepler ID')

    #add columns for gPhoton exptime
    starcat.add_column(Column(name='FUV gtime', length=len(starcat), 
                              dtype=float, unit='s'))
    starcat.add_column(Column(name='NUV gtime', length=len(starcat), 
                              dtype=float, unit='s'))
    
    #save as csv
    utils.table_to_csv(starcat, pointers.starcat_file)
else:
    starcat = Table.read(pointers.starcat_file, format='ascii')
  
#%% FETCH AND COMPUTE STELLAR PROPERTIES OF INTEREST
redo = True
if redo:
    cols = ['kic_kepler_id','kic_radius','kic_logg','kic_teff','kic_feh','g']
    propcat = catIO.fetch_MASTcatalog_data('kepler/kepler_fov',cols,
                                           sample_criteria)    
    
    #clean up column labels
    old_cols = ['Kepler_ID','Radius (solar=1.0)','Log_G',
                'Metallicity (solar=0.0)','g']
    new_cols = ['Kepler ID','R','logg','Fe/H','g mag']
    map(propcat.rename_column, old_cols, new_cols)
    propcat.sort('Kepler ID')
    
    #create a luminosity column
    L = propcat['R']**2 * (propcat['Teff']/Teff_sun)**4
    L.name = 'L'
    propcat.add_column(L)
    
    #save as csv
    utils.table_to_csv(propcat, pointers.propcat_file)
#else:
#    propcat = Table.read(pointers.propcat_file, format='ascii')
    
#%% FETCH PLANET PROPERTIES OF INTEREST
    #searching just the FUV < 18 stars yielded 2 candidates and 3 nds
redo = False
if redo:      
    catcols = ['Kepler ID','a', 'a_err1', 'a_err2', 'r_p', 'r_p_err1', 
               'r_p_err2', 'NExSci Disposition']
    dtypes = [int] + [float]*6 + [str]
    koicat = Table(masked=True, names=catcols, dtype=dtypes)
    kids = starcat['Kepler ID'].data.data
    for chunk in mynp.chunks(kids, 500):
        criteria=['kepid='+','.join(map(str, chunk)), 
                  'koi_disposition=!FALSE POSITIVE']
        cols = ['kepid','koi_sma','koi_sma_err1','koi_sma_err2','koi_prad',
                'koi_prad_err1','koi_prad_err2','koi_disposition']
        cat = catIO.fetch_MASTcatalog_data('kepler/koi', cols, criteria)
        if len(cat) > 0:
            if cat.masked:
                for row, rowmask in zip(cat, cat.mask):
                    koicat.add_row(row.data, rowmask)
            else:
                for row in cat:
                    koicat.add_row(row.data)
    
    utils.table_to_csv(koicat, pointers.koicat_file)
else:
    koicat = Table.read(pointers.koicat_file, format='ascii')
    
#%% DETERMINE THE AVAILABILITY OF G DATA
redo = False
if redo:
    for b in ['FUV','NUV']:
        for i,row in enumerate(starcat):
            ra, dec = row['RA'], row['Dec']
            tranges = gdb.fGetTimeRanges(b,[ra,dec])
            tottime = sum(tranges[:,1] - tranges[:,0])
            colname = '{} gtime'.format(b)
            starcat[colname][i] = tottime

#%% KEEP ONLY GOOD DATA
#TODO: CULL STARS THAT ARE...
    #CONTAMINATED BY FLUX FROM ANOTHER SOURCE
    #have insufficient data
if redo:
    keep_list = [propcat['logg'] > 2.5,
                 propcat['logg'] < 4.5,
                 np.logical_not(propcat['logg'].mask)]
    keep = reduce(np.logical_and,keep_list)
    good_kids = propcat['Kepler ID'][keep].data.data
    np.save(pointers.good_kids_file, good_kids)
    
#%% GENERATE AN IMAGE TO CHECK IF EXTRACTION REGIONS MATCH SOURCES

if redo:
    import imagetools as gimg
    import dbasetools as gdb
    
    ra = 292.3875 #deg
    dec = 42.3763 #deg
    s = 0.1 #deg
    
    #find the start and end of the longest exposure covering the area
    tranges = gdb.fGetTimeRanges('NUV',[ra,dec])
    imax = np.argmax(tranges[:,1] - tranges[:,0])
    gimg.write_images('NUV',[ra,dec],[tranges[imax,:]],[2*s,2*s],
                      write_cnt=pointers.miniimg_file,
                      calpath=pointers.calibration_folder)
                      
if redo:
    ra_in = np.array(minicat['RA'])
    ra_in = ra_in[np.logical_and(ra_in < (ra+s), ra_in > (ra-s))]
    dec_in = np.array(minicat['Dec'])
    dec_in = dec_in[np.logical_and(dec_in < (dec+s), dec_in > (dec-s))]
    miniimg = aplpy.FITSFigure(pointers.miniimg_file)
    miniimg.show_grayscale(stretch='arcsinh')
    miniimg.show_circles(ra_in, dec_in, radius=ap_radius, edgecolor='red', 
                         layer='apertures')

#%% GENERATE LIGHTCURVES FOR THE OBJECTS

if redo:
    if pointers.compname == 'PARKENOTEBOOK': Np = 10
    elif pointers.compname == 'COMMUNAL-PC': Np = 14
    else: Np = multiprocessing.cpu_count()*4
    call = 'python parallel_curve.py -c {} -t {} -b {} -p {} -a {} -i {} -o {}'.format(
            pointers.minicat_file, tstep, band, Np, ap_radius, bk_annulus[0], 
            bk_annulus[1])
    os.system(call)
    
#%% COMPUTE EXCESS NOISE VALUES

if redo:
    call = 'python parallel_xnoise.py -c {} -b {} -p 2'.format(
    pointers.minicat_file, band)
    os.system(call)
        
        
        
        