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

# Supply a meta file from a MeerKAT pointing and this script gives the following information:
# 1. PNG file showing the beam tiling with known pulsars within the incoherent beam
# 2. PNG file showing the a zoomed version of the known pulsar within a tiling and surrounding beams
# 3. Check for associated Fermi sources in DR2 release and report them in a csv file
# 4. Gives csv files of all pulsars in field, separation from boresight, best beam to find the pulsar
# 5. Query HESS, Gaia via webparser (?) 

############# Some fixed parameters ###################

#MeerKAT location
meerkat = EarthLocation(lat=-30.713*u.deg, lon=21.4*u.deg)

# Survey beam and incoherent beam at L-Band
survey_beam_radius = 0.2488360131953144 # From Ewan's simulations
survey_beam_area = np.pi*survey_beam_radius*survey_beam_radius
incoherent_beam_radius = 2.5*0.2488360131953144 # From Ewan's simulations
primary_beam_area = np.pi*incoherent_beam_radius*incoherent_beam_radius


#########################################################

def convert_equatorial_coordinate_to_pixel(equatorial_coordinates, bore_sight):
    """
    Credits : Weiwei's mosaic (https://github.com/wchenastro/Mosaic)
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
    Credits: Colin Clark 
    Check if given coords has a Fermi association
    """

def query_known_pulsar_ATNF(boresight_ra, boresight_dec, incoherent_beam_radius):
    """
    Credits: Matthew Pitkin (https://github.com/mattpitkin/psrqpy)
    Query known pulsars from ATNF database alone using psrqpy
    """
    q = QueryATNF(params=['JNAME','RAJ','DECJ','F0','DM'],circular_boundary=(boresight_ra,boresight_dec,incoherent_beam_radius))
    pulsar_list = np.array((q.table['JNAME','RAJ','DECJ','F0','DM']))
    return pulsar_list    
    
def query_known_pulsar_scraper(boresight_ra, boresight_dec, incoherent_beam_radius):
    """
    Query for known pulsars in field via a webparser with the pulsar survey scraper
    """
    #subprocess.check_call("python ")

def get_beam_tiling_plot(opts):
    """
    Generate beam tiling plot for a meta file 
    """


