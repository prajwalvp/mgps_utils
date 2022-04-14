import sys
import math
import subprocess
import numpy as np
from astropy import wcs
from astropy.time import Time
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
from ast import literal_eval
#import matplotlib.pyplot as plt
#from matplotlib.patches import Ellipse

"""
Map psr coordinate to gain factor via actual MGPS psf and decide if survey actually capable of discovering the pulsar 
"""


#ALL_META_PATH = '/MGPS/survey_stats/all_metafiles'
ALL_META_PATH = '/work/prajwal.voraganti/simulations'
SURVEY_BEAM_RADIUS = 1.1*0.2488 # 10 percent tolerance over the actual value


def convert_equatorial_coordinate_to_pixel(equatorial_coordinates, bore_sight):
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


def get_pixel_coherent_beam_coordinates(beam_coords, boresight_coords):
    """
    Convert coherent beam equatorial coordinates to pixel coordinates
    """
    pixel_beam_ras=[]
    pixel_beam_decs=[]
    # Convert equatorial beam coordinates to pixel coordinates
    for beam_coord in beam_coords:
        pixel_coordinates = convert_equatorial_coordinate_to_pixel(beam_coord, boresight_coords)
        pixel_beam_ras.append(boresight_coords.ra.deg + pixel_coordinates[0][0])
        pixel_beam_decs.append(boresight_coords.dec.deg + pixel_coordinates[0][1])

    return pixel_beam_ras, pixel_beam_decs





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
    
    beam_coords = SkyCoord(frame='icrs', ra= coherent_beams_ra, dec= coherent_beams_dec, unit=(u.hour, u.deg))

    return beam_coords


def pointInEllipse(x,y,xp,yp,d,D,angle):
    """
    tests if a point[xp,yp] is within
    boundaries defined by the ellipse
    of center[x,y], diameter d D, and tilted at angle
    """

    cosa=math.cos(angle)
    sina=math.sin(angle)
    dd=d/2*d/2
    DD=D/2*D/2

    a = math.pow(cosa*(xp-x)+sina*(yp-y),2)
    b = math.pow(sina*(xp-x)-cosa*(yp-y),2)
    ellipse=(a/dd)+(b/DD)
    if ellipse <= 1:
        return True
    else:
        return False


def check_discovery(meta, psr):

    with open(meta, 'r') as f:
        data = f.read()
    all_info = literal_eval(data)

    
    print "Meta file of interest:%s"%(meta)
    # Initialise plot
    #plt.clf()
    #fig, ax = plt.subplots()
    #ax.set_aspect(1)

    # utc time for sign convention issue
    utc_time = all_info['utc_start'].replace(" ","T").replace("/","-")
    print "UTC time: %s"%utc_time

    # Pointing name
    pointing_name = all_info['boresight'].split(',')[0]
    print "Pointing of interest: %s"%pointing_name

    # Boresight coordinates
    boresight_ra = all_info['boresight'].split(',')[-2]
    boresight_dec = all_info['boresight'].split(',')[-1]
    boresight_coords = SkyCoord(frame='icrs', ra=boresight_ra, dec=boresight_dec, unit=(u.hour, u.deg))
    boresight_ra_deg = boresight_coords.ra.value
    boresight_dec_deg = boresight_coords.dec.value

    # Get all coherent beam coords
    coherent_beam_coords = get_coherent_beam_coords(meta)

    # Beam shapes
    beam_width = 2.0*float(literal_eval(all_info['beamshape'])['x'])
    beam_height = 2.0*float(literal_eval(all_info['beamshape'])['y'])
    if utc_time < Time('2021-10-05T00:00:00'):
        beam_angle =  float(literal_eval(all_info['beamshape'])['angle'])
    else:
        beam_angle = -1.0*float(literal_eval(all_info['beamshape'])['angle']) # Sign change in FBFUSE tiling rollout   


    # Plot psr in pixel coordinates
    psr_coords = SkyCoord(frame='icrs', ra=psr[0], dec=psr[1], unit=(u.deg, u.deg))
    psr_pixel_coordinates = convert_equatorial_coordinate_to_pixel(psr_coords, boresight_coords)
    psr_pixel_ra = boresight_ra_deg + psr_pixel_coordinates[0][0]
    psr_pixel_dec = boresight_dec_deg + psr_pixel_coordinates[0][1]
    psr_pixel_coords = SkyCoord(frame='icrs', ra=psr_pixel_ra, dec = psr_pixel_dec, unit=(u.deg, u.deg))
    #ax.plot(psr_pixel_ra, psr_pixel_dec,'*',label='psr',markersize=7.5)


    #Plot all ellipses
    


    # Add best ellipse
    pixel_beam_ras, pixel_beam_decs = get_pixel_coherent_beam_coordinates(coherent_beam_coords, boresight_coords)
    psr_idx, psr_d2d, psr_d3d = psr_coords.match_to_catalog_sky(coherent_beam_coords)
    best_ellipse = Ellipse(xy=(pixel_beam_ras[psr_idx], pixel_beam_decs[psr_idx]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc='grey', lw=1.5)
    #ax.add_patch(best_ellipse)

    best_beam = 'beam_%d'%(psr_idx)
    print "Best beam: %s"%(best_beam)
    print "Separation: %f deg."%(psr_d2d.deg)
    # Return closest 5 beams
    all_seps = psr_coords.separation(coherent_beam_coords)
    all_beams_sorted = np.argsort(all_seps)
    ellipse_2 = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[1]], pixel_beam_decs[all_beams_sorted[1]]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc ='none', lw=1.5) 
    ellipse_3 = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[2]], pixel_beam_decs[all_beams_sorted[2]]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc='none', lw=1.5) 
    ellipse_4 = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[3]], pixel_beam_decs[all_beams_sorted[3]]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc ='none', lw=1.5) 
    ellipse_5 = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[4]], pixel_beam_decs[all_beams_sorted[4]]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc='none', lw=1.5) 
    ellipse_6 = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[5]], pixel_beam_decs[all_beams_sorted[5]]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc='none', lw=1.5) 
    #ax.add_patch(ellipse_2)
    #ax.add_patch(ellipse_3)
    #ax.add_patch(ellipse_4)
    #ax.add_patch(ellipse_5)
    #ax.add_patch(ellipse_6)

    #ax.set_xlabel('Right Ascension (Degrees)')
    #ax.set_ylabel('Declination (Degrees)')

    #ax.set_xlim(psr_pixel_ra - 3*beam_width, psr_pixel_ra + 3*beam_width)
    #ax.set_ylim(psr_pixel_dec - 3*beam_width, psr_pixel_dec + 3*beam_width)


    # Return true if the pulsar is in any of the 3 closest beams
       
    if pointInEllipse(pixel_beam_ras[psr_idx], pixel_beam_decs[psr_idx], psr_pixel_ra, psr_pixel_dec, beam_width, beam_height, beam_angle):
       #plt.savefig('%s_%s_%s_detected.png'%(utc_time, best_beam, sys.argv[1]), dpi=300)
       return True
    if pointInEllipse(pixel_beam_ras[all_beams_sorted[1]], pixel_beam_decs[all_beams_sorted[1]], psr_pixel_ra, psr_pixel_dec, beam_width, beam_height, beam_angle):
       #plt.savefig('%s_%s_%s_detected.png'%(utc_time, best_beam, sys.argv[1]), dpi=300)
       return True
    if pointInEllipse(pixel_beam_ras[all_beams_sorted[2]], pixel_beam_decs[all_beams_sorted[2]], psr_pixel_ra, psr_pixel_dec, beam_width, beam_height, beam_angle):
       #plt.savefig('%s_%s_%s_detected.png'%(utc_time, best_beam, sys.argv[1]), dpi=300)
       return True
    
    #plt.savefig('%s_%s_%s_missed.png'%(utc_time, best_beam, sys.argv[1]), dpi=300)
    return False    

def RotatedGaussian2DPDF(x, y, xMean, yMean, xSigma, ySigma, angle):
        angle = -(angle - np.pi)
        a = np.power(np.cos(angle), 2)/(2*xSigma**2) + np.power(np.sin(angle), 2)/(2*ySigma**2)
        b = - np.sin(2*angle)/(4*xSigma**2) + np.sin(2*angle)/(4*ySigma**2)
        c = np.power(np.sin(angle), 2)/(2*xSigma**2) + np.power(np.cos(angle), 2)/(2*ySigma**2)

        return np.exp(-(a*np.power(x-xMean, 2) + 2*b*(x-xMean)*(y-yMean) + c*np.power(y-yMean, 2)))


def evaulate_discovery(meta, psr):
    """
    Model the synthesised beam as a 2D Gaussian and calculate the new offset gain factor and check if the psr is still above S/N threshold.
    """
    with open(meta,'r') as f:
        data = f.read()
    all_info = literal_eval(data) 

    # Initialise plot
    #plt.clf()
    #fig, ax = plt.subplots()
    #ax.set_aspect(1)

    # utc time for sign convention issue
    utc_time = all_info['utc_start'].replace(" ","T").replace("/","-")
    print "UTC time: %s"%utc_time

    # Get boresight coords and survey beam radius
    boresight_ra = all_info['boresight'].split(',')[-2]     
    boresight_dec = all_info['boresight'].split(',')[-1]
    boresight_coords = SkyCoord(frame='icrs', ra=boresight_ra, dec=boresight_dec, unit=(u.hour, u.deg))
    boresight_ra_deg = boresight_coords.ra.deg
    boresight_dec_deg = boresight_coords.dec.deg

    # Get all coherent beam coords
    coherent_beam_coords = get_coherent_beam_coords(meta)


    # Get beam dimensions
    beam_width = 2.0*float(literal_eval(all_info['beamshape'])['x'])
    beam_height = 2.0*float(literal_eval(all_info['beamshape'])['y'])
    if utc_time < Time('2021-10-05T00:00:00'):
        beam_angle =  float(literal_eval(all_info['beamshape'])['angle'])
    else:
        beam_angle = -1.0*float(literal_eval(all_info['beamshape'])['angle']) # Sign change in FBFUSE tiling rollout   

    # sigma = FWHM / 2.355; 2 coordinates for 2d gaussian 
    x_sigma_pixel = literal_eval(all_info['beamshape'])['x']/(0.5*2.355)
    y_sigma_pixel = literal_eval(all_info['beamshape'])['y']/(0.5*2.355)
    rad_angle = np.deg2rad(beam_angle)


    # Plot psr in pixel coordinates
    psr_coords = SkyCoord(frame='icrs', ra=psr[0], dec=psr[1], unit=(u.deg, u.deg))
    psr_pixel_coordinates = convert_equatorial_coordinate_to_pixel(psr_coords, boresight_coords)
    psr_pixel_ra = boresight_ra_deg + psr_pixel_coordinates[0][0]
    psr_pixel_dec = boresight_dec_deg + psr_pixel_coordinates[0][1]
    psr_pixel_coords = SkyCoord(frame='icrs', ra=psr_pixel_ra, dec = psr_pixel_dec, unit=(u.deg, u.deg))
    #ax.plot(psr_pixel_ra, psr_pixel_dec,'*',label='psr',markersize=7.5)



    # Add best ellipse
    pixel_beam_ras, pixel_beam_decs = get_pixel_coherent_beam_coordinates(coherent_beam_coords, boresight_coords)
    psr_idx, psr_d2d, psr_d3d = psr_coords.match_to_catalog_sky(coherent_beam_coords)
    #best_ellipse = Ellipse(xy=(pixel_beam_ras[psr_idx], pixel_beam_decs[psr_idx]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc='grey', lw=1.5)
    #ax.add_patch(best_ellipse)

    pixel_beam_ra = pixel_beam_ras[psr_idx]
    pixel_beam_dec = pixel_beam_decs[psr_idx]


    best_beam = 'beam_%d'%(psr_idx)
    print "Best beam: %s"%(best_beam)
    print "Separation: %f deg."%(psr_d2d.deg)
    # Return closest 5 beams
    all_seps = psr_coords.separation(coherent_beam_coords)
    all_beams_sorted = np.argsort(all_seps)
    #ellipse_2 = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[1]], pixel_beam_decs[all_beams_sorted[1]]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc ='none', lw=1.5) 
    #ellipse_3 = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[2]], pixel_beam_decs[all_beams_sorted[2]]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc='none', lw=1.5) 
    #ellipse_4 = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[3]], pixel_beam_decs[all_beams_sorted[3]]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc ='none', lw=1.5) 
    #ellipse_5 = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[4]], pixel_beam_decs[all_beams_sorted[4]]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc='none', lw=1.5) 
    #ellipse_6 = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[5]], pixel_beam_decs[all_beams_sorted[5]]), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc='none', lw=1.5) 
    #ax.add_patch(ellipse_2)
    #ax.add_patch(ellipse_3)
        
    #ax.set_xlabel('Right Ascension (Degrees)')
    #ax.set_ylabel('Declination (Degrees)')

    #ax.set_xlim(psr_pixel_ra - 3*beam_width, psr_pixel_ra + 3*beam_width)
    #ax.set_ylim(psr_pixel_dec - 3*beam_width, psr_pixel_dec + 3*beam_width)




    # Fractional gain
    fg = RotatedGaussian2DPDF(psr_pixel_ra, psr_pixel_dec, pixel_beam_ra, pixel_beam_dec, x_sigma_pixel, y_sigma_pixel, rad_angle)

    # New SNR
    print psr[6], psr[7], fg 
    new_snr = (float(psr[6])/float(psr[7]))*fg

    # Decide if psr S/N still above threshold
    if new_snr > 9.0:
        #plt.savefig('%s_%s_%s_detected.png'%(utc_time, best_beam, sys.argv[1]))
        return True
    else:
        #plt.savefig('%s_%s_%s_missed.png'%(utc_time, best_beam, sys.argv[1]))
        return False

cnt=0
# Read MGPS pointing coordinates
pointing_coordinates_list = 'MGPS-L_all_pointing_coordinates.txt'
pointings = np.loadtxt(pointing_coordinates_list, dtype=str)[:,0]
ras = np.loadtxt(pointing_coordinates_list, dtype=str)[:,3]
decs = np.loadtxt(pointing_coordinates_list, dtype=str)[:,4]
pointing_gls = np.loadtxt(pointing_coordinates_list, dtype=str)[:,1]
pointing_gbs = np.loadtxt(pointing_coordinates_list, dtype=str)[:,2]
pointing_coords = SkyCoord(ra=ras, dec=decs, unit=(u.hourangle, u.deg), frame='icrs')
pointing_ras=pointing_coords.ra.deg
pointing_decs=pointing_coords.dec.deg


# Read the discovery file
discoveries = np.loadtxt(sys.argv[1])
gls = discoveries[:,2]
gbs = discoveries[:,3]
snrs = discoveries[:,4]
og_fgs = discoveries[:,5]

psr_discovery_coords = SkyCoord(l=gls, b=gbs, frame='galactic', unit=(u.deg, u.deg)).transform_to('icrs') 
psr_discovery_ras = psr_discovery_coords.ra.deg
psr_discovery_decs = psr_discovery_coords.dec.deg


# Catalog match 
idx, d2d, d3d = psr_discovery_coords.match_to_catalog_3d(pointing_coords)

psr_remaining=[]

# Filter pulsars based on whether they make it within the survey beam radius
for i in range(len(psr_discovery_coords)):
    if d2d.deg[i] < SURVEY_BEAM_RADIUS:
        psr_remaining.append([psr_discovery_ras[i], psr_discovery_decs[i], pointings[idx[i]], ras[idx[i]], decs[idx[i]], d2d.deg[i], snrs[i], og_fgs[i]])
        cnt+=1
print "Number of remaining discoveries based on survey beam radius cutoff: %d"%(cnt)

true_discoveries=[]

# Map pointing to meta file and get closest beam position + orientation to decide if psr within beam
for psr in psr_remaining:
    pointing = psr[2]
 
    try: 
        meta_path = subprocess.check_output("grep -rnw -l \"{}\" {}/beegfs/DATA/TRAPUM/SCI-20200703-MK-01".format(pointing, ALL_META_PATH), shell=True)    
        print "Pointing has been observed" 
    except Exception as error:
        print "Pointing has not yet been observed"
        continue
    
    if len(meta_path.decode('utf-8').split('\n')) > 2:
        print "Observed multiple times..."
        first_meta_file = meta_path.decode('utf-8').split('\n')[0]      
        second_meta_file = meta_path.decode('utf-8').split('\n')[1]
        print first_meta_file, second_meta_file

        # Check if pulsars falls within any of the beams from either epochs
        #if check_discovery(first_meta_file, psr): 
        if evaluate_discovery(first_meta_file, psr):
           print "Disovery is found by synthesis beams and considered true. Will be retained"  
           true_discoveries.append(psr)
                
        #if check_discovery(second_meta_file, psr): 
        if evaluate_discovery(second_meta_file, psr): 
           print "Disovery is found by synthesis beams and considered true. Will be retained"  
           true_discoveries.append(psr)  

    else: 
        print "Observed once..."    
        meta_file = meta_path.decode('utf-8').split('\n')[0]        
        #if check_discovery(meta_file, psr):
        if evaulate_discovery(meta_file, psr):
           print "Disovery is found by synthesis beams and considered true. Will be retained"  
           true_discoveries.append(psr)  
        else:
           print "Discovery outside the coherent beam region. Will be ignored" 

with open('true_discoveries_%s.txt'%(sys.argv[1]), 'w') as f:
    for item in true_discoveries:
        f.write("%s\n" % item)
    f.close()
    
