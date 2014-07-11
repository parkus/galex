# -*- coding: utf-8 -*-
"""
Created on Tue May 06 12:13:34 2014

@author: Parke
"""
#import pdb

def fetch_MASTcatalog_data(catalog_url,columns,selection_criteria=[],
                           maxrecords=5e5):
    """Fetches Kepler-GALEX cross-match data for gold-standard (confidently
    matched) stars and returns a pandas database.
    
    The catalog URL should NOT have '/' on either end. Ex. 'kepler/kgmatch' 
    See http://archive.stsci.edu/vo/mast_services.html#MISSION for a list of
    the catalogs that can be searched in this way.
    
    You can use > and < in the criteria. They will be converted to the proper
    HTML characters. 
    
    Be careful of escape characters (or just use string literals). E.g., use 
    ['parameter=!\\null'] or [r'parameter=!\null'] instead of [r'parameter=\null'].
    """
    
    from urllib import urlopen
    from astropy.table import Table
    
    selection_criteria = [s.replace('>','=%3E') for s in selection_criteria]
    selection_criteria = [s.replace('<','=%3C') for s in selection_criteria]
    
    #create root URL
    url =  'http://archive.stsci.edu/' + catalog_url + '/search.php?'
    #first the criteria for what entries to return
    url += '&'.join(selection_criteria) if selection_criteria else ''
    #now what data to return for those entries
    url += '&selectedColumnsCsv=' + ','.join(columns)
    #return it as a CSV
    url += '&outputformat=CSV'
    #get only so many records (couldn't figure out syntax for all)
    url += '&max_records=' + str(int(maxrecords))
    #slap on the necessary action=search bit
    url += '&action=Search'  
    
#    pdb.set_trace()
    f = urlopen(url)
    data = Table.read(f,format='ascii',delimiter=',',guess=False,data_start=2)
    return data