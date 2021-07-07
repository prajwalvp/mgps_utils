#Modified version of Prajwalp's pointing_pattern.py that includes both the survey beam and the inoherent beam.
#In addition of searching for PSRCAT hits, it also looks for FermiLAT candidates.
#It requires that the direcotry with table 4FGL_DR2_Ppsr.fits is specifyied in line 92.
import re
import sys
import glob
import optparse
import numpy as np
import astropy.units as u
from astropy import wcs
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import ICRS, Galactic, SkyCoord, EarthLocation, AltAz
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Circle
from ast import literal_eval
from psrqpy import QueryATNF


# MeerKAT location
meerkat = EarthLocation(lat=-30.713*u.deg, lon=21.4*u.deg)

# Survey beam and incoherent beam at L-Band 
survey_beam_radius = 0.2488360131953144 # From Ewan's simulations
survey_beam_area = np.pi*survey_beam_radius*survey_beam_radius
incoherent_beam_radius = 2.5*0.2488360131953144 # From Ewan's simulations
primary_beam_area = np.pi*incoherent_beam_radius*incoherent_beam_radius 

# Plotting ornaments
plt.rc('axes', labelsize=12)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
plt.rc('ytick', labelsize=8) 
plt.rcParams['figure.figsize'] = (10.0, 10.0)
plt.rcParams['font.family'] = "serif"
plt.rcParams['figure.dpi'] = 200

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

    return pixel_coordinates

def check_Fermi_association(coords):
    """
    Check if given coords has a Fermi association
    """
    # Add FermiLAT sources with a high chance of being a pulsar.

    fgl4_fits = fits.open('/home/miquel/Desktop/MPIfR/TRAPUM/gamma_pulsars/4FGL_DR2_Ppsr.fits')[1].data
    fgl4_pos = SkyCoord(fgl4_fits.field('RAJ2000'),
                        fgl4_fits.field('DEJ2000'),
                        unit=u.deg)
    fgl4_name = fgl4_fits.field('Source_Name')
    fgl4_val = fgl4_fits.field('P(psr)')
    gamma_list = []
    i=0
    for position in fgl4_pos:
        if position.separation(boresight_coords).deg<=survey_beam_radius*1.1 and fgl4_val[i]>=0.5:
            gamma_list.append([fgl4_name[i].split()[1],fgl4_val[i],position])
            print(fgl4_name[i].split()[1]+", with P="+str(fgl4_val[i])+" found at "+str(position))
        i=i+1

    # Convert good FermiLAT candidates to pixel too
    for cand in gamma_list:
        cand_coords = cand[2]
        cand_pixel_coordinates = convert_equatorial_coordinate_to_pixel(cand_coords,boresight_coords)
        cand_pixel_ra = boresight_ra_deg + cand_pixel_coordinates[0][0]
        cand_pixel_dec = boresight_dec_deg + cand_pixel_coordinates[0][1]
        ax.plot(cand_pixel_ra, cand_pixel_dec,'*',label=cand[0]+", P="+str(cand[1]),markersize=10.0)  
    


def query_known_pulsar_ATNF(boresight_ra, boresight_dec, incoherent_beam_radius):
    """
    Query known pulsars from ATNF database alone using psrqpy
    """
    q = QueryATNF(params=['JNAME','RAJ','DECJ'],circular_boundary=(boresight_ra,boresight_dec,incoherent_beam_radius))
    pulsar_list = np.array((q.table['JNAME','RAJ','DECJ']))
    return pulsar_list


def query_known_pulsar_scraper(boresight_ra, boresight_dec, incoherent_beam_radius):
    """
    Query for known pulsars in field via a webparser with the pulsar survey scraper 
    """     
    #subprocess.check_call("python ")


def get_tiling_plot(opts):
    """
    Generate tiling plot for a meta file
    """
    plt.clf()
    fig, ax = plt.subplots()
    ax.set_aspect(1)
    
    with open(opts.meta,'r') as f:
        data = f.read()
    all_info = literal_eval(data)   
    
    vals = list(all_info['beams'].values())
    keys = list(all_info['beams'].keys())
     
    # Boresight
    boresight_ra = all_info['boresight'].split(',')[-2]     
    boresight_dec = all_info['boresight'].split(',')[-1]
    pointing_name = all_info['boresight'].split(',')[0]
    boresight_coords = SkyCoord(frame='icrs', ra=boresight_ra, dec=boresight_dec, unit=(u.hour, u.deg))
    boresight_ra_deg = boresight_coords.ra.value
    boresight_dec_deg = boresight_coords.dec.value

    # Location on sky
    pointing_coords = SkyCoord(frame='icrs', ra=boresight_ra, dec=boresight_dec, unit=(u.hour, u.deg))
    
    # Beam shapes    
    beam_width = 2.0*float(literal_eval(all_info['beamshape'])['x'])
    beam_height = 2.0*float(literal_eval(all_info['beamshape'])['y'])
    beam_angle = literal_eval(all_info['beamshape'])['angle']
    
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
    
        ellipse = Ellipse(xy=(pixel_beam_ra, pixel_beam_dec), 
                          width=beam_width, 
                          height=beam_height,
                          angle = beam_angle, 
                          edgecolor='blue', 
                          fc='none', 
                          lw=1.5)
        ax.add_patch(ellipse)
        
        # Add beam numbers
        ax.annotate(beam_no, (pixel_beam_ra,pixel_beam_dec),fontsize=3)  
    
   
    # Known pulsar plotting
    if opts.query_atnf:
        pulsar_list = query_known_pulsar_ATNF(boresight_ra, boresight_dec, incoherent_beam_radius)
 
    elif opts.query_pss:
        pulsar_list = query_known_pulsar_pss(boresight_ra, boresight_dec, incoherent_beam_radius)

    # Convert known pulsar from equatorial to pixel too
    if pulsar_list in locals():
        print("Plotting the pulsars in incoherent beam field")    
        for psr in pulsar_list:
            psr_coords = SkyCoord(frame='icrs', ra=psr[1], dec=psr[2], unit=(u.hour, u.deg))
            psr_pixel_coordinates = convert_equatorial_coordinate_to_pixel(psr_coords,boresight_coords)
            psr_pixel_ra = boresight_ra_deg + psr_pixel_coordinates[0][0]
            psr_pixel_dec = boresight_dec_deg + psr_pixel_coordinates[0][1]
            ax.plot(psr_pixel_ra, psr_pixel_dec,'*',label=psr[0],markersize=10.0)    


    #Forcing the plot by adding a fake pulsar at the center.
    psr= [boresight_ra,boresight_dec]
    psr_coords = SkyCoord(frame='icrs', ra=psr[0], dec=psr[1], unit=(u.hour, u.deg))
    psr_pixel_coordinates = convert_equatorial_coordinate_to_pixel(psr_coords,boresight_coords)
    psr_pixel_ra = boresight_ra_deg + psr_pixel_coordinates[0][0]
    psr_pixel_dec = boresight_dec_deg + psr_pixel_coordinates[0][1]
    ax.plot(psr_pixel_ra, psr_pixel_dec,'*',markersize=0.0) 
         
    # Survey beam radius    
    circle1 = Circle((boresight_ra_deg,boresight_dec_deg), 
                      survey_beam_radius, 
                      color='red',
                      linestyle='--',
                      linewidth=2.5,
                      fill=False,
                      label='Survey beam')
    ax.add_patch(circle1)

    # Incoherent beam radius
    circle1 = Circle((boresight_ra_deg,boresight_dec_deg), 
                      incoherent_beam_radius, 
                      color='green',
                      linestyle='--',
                      linewidth=2.5,
                      fill=False,
                      label='Incoherent beam')
    ax.add_patch(circle1)


    # Add extra specified radius
    if opts.user_radius != None:
        user_circle = Circle((psr_pixel_ra,psr_pixel_dec),
                              opts.user_radius,
                              color='black',
                              linestyle='--',
                              linewidth=1.5,
                              fill=False,
                              label=opts.user_radius_name) 

    # UTC Time
    utc_time = all_info['utc_start'].replace(" ","T").replace("/","-")
    time = Time(utc_time)
    
    # Calculate telescope alt/elevation
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
    plt.legend(prop={"size":6})
    plt.savefig('{}.png'.format(pointing_name),dpi=200)


if __name__ =="__main__":
    # Set options
    parser = optparse.OptionParser()
    parser.add_option('--meta', type=str, help='Path to meta file for pointing', dest='meta')
    parser.add_option('--fermi', type=int, help='Flag to find Fermi sources in field', dest='fermi', default=0)
    parser.add_option('--fermi_fits', type=str, help='Path to FITS file of Fermi sources', dest='fermi_fits', default=None)
    parser.add_option('--query_atnf', type=int, help='Flag to query ATNF pulsars in the field', dest='query_atnf', default=0)
    parser.add_option('--query_pss', type=int, help='Flag to query all pulsars in the field', dest='query_pss', default=1)
    parser.add_option('--user_radius', type=float, help='Specify a beam radius to compare with tiling', dest='user_radius', default=None)
    parser.add_option('--user_radius_name', type=str, help='Name of beam radius (e.g. Parkes, GBT)', dest='user_radius_name', default=None)
    opts, args = parser.parse_args()

    # Generate Pointing pattern plot
    get_tiling_plot(opts)
    
   

         
