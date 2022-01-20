import subprocess
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky

"""
Map psr coordinate to gain factor via actual MGPS psf and decide if survey actually capable of discovering the pulsar 
"""

def convert_equatorial_coordinate_to_pixel(equatorial_coordinates, bore_sight, utc_time):
    """
    Credits: Weiwei's Mosaic
    https://docs.astropy.org/en/stable/wcs/index.html#using-astropy-wcs
    """
    step = 1/10000000000.

    wcs_properties = wcs.WCS(naxis=2)
    wcs_properties.wcs.crpix = [0, 0]
    wcs_properties.wcs.cdelt = [step, step]
    wcs_properties.wcs.crval = [bore_sight.ra.deg,bore_sight.dec.deg]
    wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    eq_coordinates = np.array([[equatorial_coordinates.ra.deg,equatorial_coordinates.dec.deg]])    
    scaled_pixel_coordinates = wcs_properties.wcs_world2pix(eq_coordinates, 0)
    pixel_coordinates = scaled_pixel_coordinates*step 
     
    return pixel_coordinates


def get_coherent_beam_coords(meta):
    """
    Return beam coordinates for all coherent beams
    """     
    with open(meta,'r') as f:
        data = f.read()
    all_info = literal_eval(data)  
 
    # Get all key value pairs for beams and sort them based on beam number
    vals = list(all_info['beams'].values())
    keys = list(all_info['beams'].keys())
    vals = [x for _, x in sorted(zip(keys,vals))]
    keys = sorted(keys)
   
    
    # Convert beam coordinates into Astropy dataframe
    coherent_beams_ra = [] 
    coherent_beams_dec = [] 

    for i in range(len(vals)):
        if 'unset' in vals[i]:
            continue
        coherent_beams_ra.append(vals[i].split(',')[-2])
        coherent_beams_dec.append(vals[i].split(',')[-1])
    
    # Convert equatorial beam coordinates to pixel coordinates
    beam_coords = SkyCoord(frame='icrs', ra= coherent_beams_ra, dec= coherent_beams_dec, unit=(u.hour, u.deg))

    return beam_coords


ALL_META_PATH = '/MGPS/survey_stats/all_metafiles'
SURVEY_BEAM_RADIUS = 1.1*0.2488 # 10 percent tolerance over the actual value


cnt=0
# Read MGPS pointing coordinates
pointings = np.loadtxt('MGPS-L_all_pointing_coordinates.txt', dtype=str)[:,0]
ras = np.loadtxt('MGPS-L_all_pointing_coordinates.txt', dtype=str)[:,3]
decs = np.loadtxt('MGPS-L_all_pointing_coordinates.txt', dtype=str)[:,4]
pointing_gls = np.loadtxt('MGPS-L_all_pointing_coordinates.txt', dtype=str)[:,1]
pointing_gbs = np.loadtxt('MGPS-L_all_pointing_coordinates.txt', dtype=str)[:,2]
pointing_coords = SkyCoord(ra=ras, dec=decs, unit=(u.hourangle, u.deg), frame='icrs')
pointing_ras=pointing_coords.ra.deg
pointing_decs=pointing_coords.dec.deg


# Read the discovery file
discoveries = np.loadtxt('soi.disc')
gls = discoveries[:,2]
gbs = discoveries[:,3]
psr_discovery_coords = SkyCoord(l=gls, b=gbs, frame='galactic', unit=(u.deg, u.deg)).transform_to('icrs') 
psr_discovery_ras = psr_discovery_coords.ra.deg
psr_discovery_decs = psr_discovery_coords.dec.deg


# Catalog match 
idx, d2d, d3d = psr_discovery_coords.match_to_catalog_3d(pointing_coords)

psr_remaining=[]

# Filter pulsars based on whether they make it within the survey beam radius
for i in range(len(psr_discovery_coords)):
    if d2d.deg[i] < SURVEY_BEAM_RADIUS:
        psr_remaining.append([psr_discovery_ras[i], psr_discovery_decs[i], pointings[idx[i]], ras[idx[i]], decs[idx[i]], d2d.deg[i]])
        cnt+=1
print "Number of remaining discoveries based on survey beam radius cutoff: %d"%(cnt)


# Map pointing to meta file and get closest beam position + separation
for psr in psr_remaining:
    pointing = psr[2]
 
    try: 
        meta_path = subprocess.check_output("grep -rnw -l \"{}\" {}/beegfs/DATA/TRAPUM/SCI-20200703-MK-01".format(pointing, ALL_META_PATH), shell=True)    
        print "Pointing has been observed" 
    except Exception as error:
        print "Pointing has not yet been observed"
        continue
    
    if len(meta_path.decode('utf-8').split('\n')) > 2:
        print "Observed multiple times"
        first_meta_file = meta_path.decode('utf-8').split('\n')[0]      
         second_meta_file = meta_path.decode('utf-8').split('\n')[1]
        print (first_meta_file, second_meta_file)
        with open(first_meta_file, 'r') as f:
            data = f.read()
        all_info = literal_eval(data)
        first = all_info['output_dir']
        with open(second_meta_file, 'r') as f:
            data = f.read()
        all_info = literal_eval(data)
        second = all_info['output_dir']


