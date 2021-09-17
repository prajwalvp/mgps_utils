#Modified version of Prajwalp's pointing_pattern.py that includes both the survey beam and the inoherent beam.
#In addition of searching for PSRCAT hits, it also looks for FermiLAT candidates.
#It requires that the direcotry with table 4FGL_DR2_Ppsr.fits is specifyied in line 92.

import re
import glob
import sys
import numpy as np
from astropy import wcs
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import ICRS, Galactic, SkyCoord, EarthLocation, AltAz
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Circle
from ast import literal_eval
from psrqpy import QueryATNF
#from astropy.utils.iers import conf
#conf.auto_max_age = None

# MeerKAT location
meerkat = EarthLocation(lat=-30.713*u.deg, lon=21.4*u.deg)

#Constants
survey_beam_radius = 0.2488 # FWHM/sqrt(5)
parkes_beam_radius = 0.5*0.23333333 # beam width is 14 arcmin at L-Band for Parkes
survey_beam_area = np.pi*survey_beam_radius*survey_beam_radius
parkes_beam_area = np.pi*parkes_beam_radius*parkes_beam_radius

# From Weiwei's Mosaic
def convert_equatorial_coordinate_to_pixel(equatorial_coordinates, bore_sight):
    """
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
    pixel_coordinates = scaled_pixel_coordinates * step
 
    #print (bore_sight.ra.deg, bore_sight.dec.deg)
    #print (scaled_pixel_coordinates)
    #print (step)
     
    return pixel_coordinates


#Boresight coordinates
survey_beam_radius = 0.2488360131953144 # From Ewan's simulations
incoherent_beam_radius = 2.5*0.2488360131953144 # From Ewan's simulations
plt.rc('axes', labelsize=12)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
plt.rc('ytick', labelsize=8) 
plt.rcParams['figure.figsize'] = (10.0, 10.0)
plt.rcParams['font.family'] = "serif"
plt.rcParams['figure.dpi'] = 200



def get_tiling_plot(meta):
    with open(meta,'r') as f:
        data = f.read()
    all_info = literal_eval(data)   

    plt.clf()
    fig, ax = plt.subplots()
    ax.set_aspect(1)    
     
    # Boresight
    boresight_ra = all_info['boresight'].split(',')[-2]     
    boresight_dec = all_info['boresight'].split(',')[-1]
    pointing_name = all_info['boresight'].split(',')[0]
    boresight_coords = SkyCoord(frame='icrs', ra=boresight_ra, dec=boresight_dec, unit=(u.hour, u.deg))
    boresight_ra_deg = boresight_coords.ra.value
    boresight_dec_deg = boresight_coords.dec.value


    # Beam shapes    
    beam_width = 2.0*float(literal_eval(all_info['beamshape'])['x'])
    beam_height = 2.0*float(literal_eval(all_info['beamshape'])['y'])
    beam_angle = literal_eval(all_info['beamshape'])['angle']
    
    

    vals = list(all_info['beams'].values())
    keys = list(all_info['beams'].keys())

    # Add known pulsar
    q = QueryATNF(params=['JNAME','RAJ','DECJ'],circular_boundary=(boresight_ra,boresight_dec,incoherent_beam_radius))
    pulsar_list = np.array((q.table['JNAME','RAJ','DECJ']))


    # Add ellipses

    for i in range(len(vals)):
        beam_ra = vals[i].split(',')[-2]
        beam_dec = vals[i].split(',')[-1]
        beam_no = keys[i][-3:]
    
        # Convert equatorial beam coordinates to pixel coordinates
        beam_coords = SkyCoord(frame='icrs', ra=beam_ra, dec=beam_dec, unit=(u.hour, u.deg))
        pixel_coordinates = convert_equatorial_coordinate_to_pixel(beam_coords,boresight_coords)
        pixel_beam_ra = boresight_ra_deg + pixel_coordinates[0][0]
        pixel_beam_dec = boresight_dec_deg + pixel_coordinates[0][1]
    
    
        ellipse = Ellipse(xy=(pixel_beam_ra, pixel_beam_dec), width=beam_width, height=beam_height,angle = beam_angle, edgecolor='blue', fc='none', lw=1.5)
        ax.add_patch(ellipse)
        # Add beam numbers
        ax.annotate(beam_no, (pixel_beam_ra,pixel_beam_dec),fontsize=3)  
    
    
    # Convert known pulsar from equatorial to pixel too
    for psr in pulsar_list:
        psr_coords = SkyCoord(frame='icrs', ra=psr[1], dec=psr[2], unit=(u.hour, u.deg))
        psr_pixel_coordinates = convert_equatorial_coordinate_to_pixel(psr_coords,boresight_coords)
        psr_pixel_ra = boresight_ra_deg + psr_pixel_coordinates[0][0]
        psr_pixel_dec = boresight_dec_deg + psr_pixel_coordinates[0][1]
        ax.plot(psr_pixel_ra, psr_pixel_dec,'*',label=psr[0],markersize=10.0)    

    # Convert good FermiLAT candidates to pixel too
    #for cand in gamma_list:
    #    cand_coords = cand[2]
    #    cand_pixel_coordinates = convert_equatorial_coordinate_to_pixel(cand_coords,boresight_coords)
    #    cand_pixel_ra = boresight_ra_deg + cand_pixel_coordinates[0][0]
    #    cand_pixel_dec = boresight_dec_deg + cand_pixel_coordinates[0][1]
    #    ax.plot(cand_pixel_ra, cand_pixel_dec,'*',label=cand[0]+", P="+str(cand[1]),markersize=10.0)  

    #Adding Vishnu's pulsars.
#    pulsar_vishnu=[["J1632 (cand)","16:32:27.47","-45:20:02.8"],["J1632 (follow-up)","16:32:27.47","-45:02:02.8"]]
#    for psr in pulsar_vishnu:
#        psr_coords = SkyCoord(frame='icrs', ra=psr[1], dec=psr[2], unit=(u.hour, u.deg))
#        psr_pixel_coordinates = convert_equatorial_coordinate_to_pixel(psr_coords,boresight_coords)
#        psr_pixel_ra = boresight_ra_deg + psr_pixel_coordinates[0][0]
#        psr_pixel_dec = boresight_dec_deg + psr_pixel_coordinates[0][1]
#        ax.plot(psr_pixel_ra, psr_pixel_dec,'*',label=psr[0],markersize=10.0)   

    #Forcing the plot by adding a fake pulsar at the center.
    psr= [boresight_ra,boresight_dec]
    psr_coords = SkyCoord(frame='icrs', ra=psr[0], dec=psr[1], unit=(u.hour, u.deg))
    psr_pixel_coordinates = convert_equatorial_coordinate_to_pixel(psr_coords,boresight_coords)
    psr_pixel_ra = boresight_ra_deg + psr_pixel_coordinates[0][0]
    psr_pixel_dec = boresight_dec_deg + psr_pixel_coordinates[0][1]
    ax.plot(psr_pixel_ra, psr_pixel_dec,'*',markersize=0.0) 
         
    # Survey beam radius    
    circle1 = Circle((boresight_ra_deg,boresight_dec_deg), survey_beam_radius, color='red',linestyle='--',linewidth=2.5,fill=False,label='Survey beam')
    ax.add_patch(circle1)

    # Incoherent beam radius
    circle1 = Circle((boresight_ra_deg,boresight_dec_deg), incoherent_beam_radius, color='green',linestyle='--',linewidth=2.5,fill=False,label='Incoherent beam')
    ax.add_patch(circle1)


    # Parkes beam diameter    
    #circle_parkes = Circle((psr_pixel_ra,psr_pixel_dec), parkes_beam_radius, color='purple',linestyle='--',linewidth=1.5,fill=False,label='Parkes beam')
    #ax.add_patch(circle_parkes)

    ### Get elevation 
    
    # Location on sky
    pointing_coords = SkyCoord(frame='icrs', ra=boresight_ra, dec=boresight_dec, unit=(u.hour, u.deg))

    # UTC Time
    #pattern = "/metafiles/(.*?).meta"
    utc_time = all_info['utc_start'].replace(" ","T").replace("/","-")
    #utc_time = re.search(pattern, meta).group(1)
    time = Time(utc_time)
    
    # Calculate alt/elevation
    coord_altaz = pointing_coords.transform_to(AltAz(obstime=time,location=meerkat))
    elv_value = coord_altaz.alt.deg
    
    ### Survey beam filling ratio
    no_of_beams = len(vals) - 1 # Remove incoherent beam
    total_area = np.pi*0.25*beam_width*beam_height*no_of_beams
    survey_beam_fill_factor = (total_area/survey_beam_area)/0.91
    

    #plotting ornaments
    ax.set_xlabel('Right Ascension (Degrees)')
    ax.set_ylabel('Declination (Degrees)')
    ax.set_title("Pointing: %s, Elevation: %f deg., SBCF=%f "%(pointing_name, elv_value, survey_beam_fill_factor))
    #ax.set_xlim(134.3,134.4)
    #ax.set_ylim(-44.35,-44.25)
    plt.legend(prop={"size":6})
    plt.savefig('{}.png'.format(pointing_name),dpi=300)
    #plt.sh()



#files=glob.glob(sys.argv[1])
#for file in files:
#	print("")
#	print(file)
get_tiling_plot(sys.argv[1])
