from astropy import wcs
from matplotlib.patches import Ellipse,Circle
import re
import astropy.units as u
import re
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import ICRS, Galactic
from ast import literal_eval
from psrqpy import QueryATNF

# MeerKAT location
meerkat = EarthLocation(lat=-30.713*u.deg, lon=21.4*u.deg)
survey_beam_radius = 0.2488 # FWHM/sqrt(5)
primary_beam_area = np.pi*survey_beam_radius*survey_beam_radius

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


#Boresight coordinates
survey_beam_radius = 0.2488360131953144 # From simulations
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14) 
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
    q = QueryATNF(params=['JNAME','RAJ','DECJ'],circular_boundary=(boresight_ra,boresight_dec,survey_beam_radius*np.sqrt(5)))
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
        ax.plot(psr_pixel_ra, psr_pixel_dec,'*',label=psr[0],color='k',markersize=10.0)    
         
    # Survey beam radius    
    circle1 = Circle((boresight_ra_deg,boresight_dec_deg), survey_beam_radius, color='red',linestyle='--',linewidth=2.5,fill=False,label='Survey beam')
    ax.add_patch(circle1)

    ### Get elevation 
    
    # Location on sky
    pointing_coords = SkyCoord(frame='icrs', ra=boresight_ra, dec=boresight_dec, unit=(u.hour, u.deg))

    # UTC Time
    pattern = "metafiles/(.*?).meta"
    utc_time = re.search(pattern, meta).group(1)
    time = Time(utc_time)
    
    # Calculate alt/elevation
    coord_altaz = pointing_coords.transform_to(AltAz(obstime=time,location=meerkat))
    elv_value = coord_altaz.alt.deg
    
    ### Survey beam filling ratio
    no_of_beams = len(vals) - 1 # Remove incoherent beam
    total_area = np.pi*0.25*beam_width*beam_height*no_of_beams
    survey_beam_fill_factor = (total_area/primary_beam_area)/0.91
    

    #plotting ornaments
    ax.set_xlabel('Right Ascension (Degrees)')
    ax.set_ylabel('Declination (Degrees)')
    ax.set_title("Pointing: %s, Elevation: %f deg., SBCF=%f "%(pointing_name, elv_value, survey_beam_fill_factor))
    plt.legend()
    plt.show()
