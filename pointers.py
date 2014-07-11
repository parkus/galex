# -*- coding: utf-8 -*-
"""
Created on Fri May 16 14:11:40 2014

@author: Parke
"""
import sys
from os import environ

compname = environ['COMPUTERNAME']
gPhoton_folder = '../gPhoton/source/'
calibration_folder = '../gPhoton/cal/'

mini = True
cat_folder = '../Scratchwork/minitest/' if mini else '../Catalogs/'

starcat_file = cat_folder + 'catalog_galex.csv' 
propcat_file = cat_folder + 'catalog_stellar_properties.csv'
planetcat_file = cat_folder + 'catalog_planet_properties.csv'
koicat_file = cat_folder + 'catalog_koi_properties.csv'
good_kids_file = cat_folder + 'list_good_stars.npy'
miniimg_file = cat_folder + 'Scratchwork/image_count.fits'
xnoisecat_file = cat_folder + 'catalog_xnoise.csv'
    

if compname == 'PARKENOTEBOOK':
    curve_folder = 'c:\\Users\\Parke\\Documents\\Grad School\\GALEX\\Lightcurves\\'
    global_code_folder = 'c:\\Users\\Parke\\Google Drive\\Python'
elif compname == 'COMMUNAL-PC':
    curve_folder = 'e:\\GALEX Data\\'
    global_code_folder = 'e:\\Google Drive\\Python'
	
else:
    curve_folder = '../Lightcurves/'
    print 'Not sure where the folder for lightcurves is on this system. Using',
    print curve_folder,'by default.'

if sys.path.count(gPhoton_folder) == 0: sys.path.append(gPhoton_folder)
if sys.path.count(global_code_folder) == 0: sys.path.append(global_code_folder)
    
def curvefile(kid, band):
    return '{}{}.{}.csv'.format(curve_folder, kid, band)