#!/usr/bin/python
import os
import sys
import math
import glob
import json
import getpass
import logging
import subprocess
import optparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle
import pandas as pd
from psrqpy import QueryATNF
import astropy.units as u
from astropy import wcs
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import ICRS, Galactic, SkyCoord, EarthLocation, AltAz, match_coordinates_sky


"""
This script reads in a given apsuse.meta file for a TRAPUM/MeerKAT observation and does the following:

a) Generates a plot of the coherent beam tiling marked with the survey beam and incoherent beam and known pulsars within the incoherent beam
b) Looks for possible associations with Fermi DR2, reports them and plots r95 regions on the beam tiling plots
c) Takes in a user defined coordinate e.g. 07:37 -37:39 and plots it with the tiling
d) Takes in a user defined beam radius to help compare the tiling with other telescope beams. 
e) Cross matches with ATNF or Pulsar scraper to get possible known pulsars expected
f) Generates a beam list to retain based on known pulsars expected in the vicinity: (localisation should be good)
g) Reads in customised google spreadsheets for unpublished known pulsar cross matching (e.g. HTRU reprocessing discoveries) 

More features and suggestions are welcome!
"""

############  CONSTANTS ############################
# MeerKAT location
meerkat = EarthLocation(lat=-30.713*u.deg, lon=21.4*u.deg)

# Beam dimensions
# beam width is 14 arcmin at L-Band for Parkes
parkes_beam_radius = 0.5*0.23333333
survey_beam_radius = 0.2488360131953144  # FWHM/sqrt(5)
# survey_beam_radius = 0.42536 # FWHM/sqrt(5)
# effelsberg_beam_radius
survey_beam_area = np.pi*survey_beam_radius*survey_beam_radius
parkes_beam_area = np.pi*parkes_beam_radius*parkes_beam_radius
incoherent_beam_radius = math.sqrt(5)*survey_beam_radius


########### Matplotlib settings ##############################
plt.rc('axes', labelsize=12)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
plt.rc('ytick', labelsize=8)
plt.rcParams['figure.figsize'] = (10.0, 10.0)
plt.rcParams['font.family'] = "serif"
plt.rcParams['figure.dpi'] = 400


###############################################################

# Logging format
log = logging.getLogger('meta_extractor')
FORMAT = "[%(levelname)s - %(asctime)s - %(filename)s:%(lineno)s] %(message)s"
logging.basicConfig(format=FORMAT)
log.setLevel('INFO')


def convert_equatorial_coordinate_to_pixel(equatorial_coordinates, bore_sight, utc_time):
    """
    Credits: Weiwei's Mosaic
    https://docs.astropy.org/en/stable/wcs/index.html#using-astropy-wcs
    """
    step = 1/10000000000.

    wcs_properties = wcs.WCS(naxis=2)
    wcs_properties.wcs.crpix = [0, 0]
    wcs_properties.wcs.cdelt = [step, step]
    wcs_properties.wcs.crval = [bore_sight.ra.deg, bore_sight.dec.deg]
    wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    eq_coordinates = np.array(
        [[equatorial_coordinates.ra.deg, equatorial_coordinates.dec.deg]])
    scaled_pixel_coordinates = wcs_properties.wcs_world2pix(eq_coordinates, 0)
    pixel_coordinates = scaled_pixel_coordinates*step

    return pixel_coordinates


def nbeams(data):
    """
    Get number of coherent beams recorded from meta file
    """
    n = 0
    for name, target in data["beams"].items():
        if name.startswith("cfbf") and "unset" not in target:
            n += 1
    return n


def get_coherent_beam_coords(meta):
    """
    Return beam coordinates for all coherent beams
    """
    with open(meta, 'r') as f:
        all_info = json.load(f)

    # Get all key value pairs for beams and sort them based on beam number
    vals = list(all_info['beams'].values())
    keys = list(all_info['beams'].keys())
    vals = [x for _, x in sorted(zip(keys, vals))]
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
    beam_coords = SkyCoord(frame='icrs', ra=coherent_beams_ra,
                           dec=coherent_beams_dec, unit=(u.hour, u.deg))

    return beam_coords


def get_pixel_coherent_beam_coordinates(beam_coords, boresight_coords, utc_time):
    """
    Convert coherent beam equatorial coordinates to pixel coordinates
    """
    pixel_beam_ras = []
    pixel_beam_decs = []
    # Convert equatorial beam coordinates to pixel coordinates
    for beam_coord in beam_coords:
        pixel_coordinates = convert_equatorial_coordinate_to_pixel(
            beam_coord, boresight_coords, utc_time)
        pixel_beam_ras.append(boresight_coords.ra.deg +
                              pixel_coordinates[0][0])
        pixel_beam_decs.append(
            boresight_coords.dec.deg + pixel_coordinates[0][1])

    return pixel_beam_ras, pixel_beam_decs


def calculate_ellipse_coefficients(angle, x_sigma, y_sigma):
    """
    Calculate the coefficients (a, b, c) for the quadratic form of an ellipse.

    Args:
        angle (float): The orientation of the ellipse in radians.
        x_sigma (float): The semi-major axis of the ellipse (std. dev).
        y_sigma (float): The semi-minor axis of the ellipse (std. dev).

    Returns:
        tuple: Coefficients (a, b, c) for the quadratic form.
    """
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)
    x_sigma_sq = x_sigma * x_sigma
    y_sigma_sq = y_sigma * y_sigma

    a = (cos_angle**2) / (2 * x_sigma_sq) + (sin_angle**2) / (2 * y_sigma_sq)
    b = (-math.sin(2 * angle)) / (4 * x_sigma_sq) + \
        (math.sin(2 * angle)) / (4 * y_sigma_sq)
    c = (sin_angle**2) / (2 * x_sigma_sq) + (cos_angle**2) / (2 * y_sigma_sq)

    return a, b, c


def gaussian_ellipse_value(x, y, x_c, y_c, angle, x_sigma, y_sigma, G0=1.0):
    """
    Compute the Gaussian intensity at (x, y) with respect to an elliptical 2D Gaussian beam.

    Args:
        x, y (float): Coordinates where the value is calculated.
        x_c, y_c (float): Center of the ellipse.
        angle (float): Orientation of the ellipse in radians.
        x_sigma (float): Semi-major axis (std. dev).
        y_sigma (float): Semi-minor axis (std. dev).
        G0 (float): Peak Gaussian value at the center (default: 1.0).

    Returns:
        float: Gaussian intensity at (x, y).
    """
    # Compute ellipse coefficients
    a, b, c = calculate_ellipse_coefficients(angle, x_sigma, y_sigma)

    # Compute the quadratic form value D(x, y)
    dx = x - x_c
    dy = y - y_c
    D = a * dx**2 + 2 * b * dx * dy + c * dy**2

    # Compute Gaussian value at (x, y)
    G_xy = G0 * np.exp(-0.5 * D)

    return G_xy


def get_Fermi_rq_pulsars(opts, boresight_coords, pointing_name, utc_time, meta_output_path):
    """
    Write out Fermi radio quiet pulsars in field
    """
    try:
        rq_df = pd.read_csv(
            '{}/Fermi_radio_quiet_pulsars_public.csv'.format(os.getcwd()))
    except Exception as error:
        return None

    fermi_rq_pos = SkyCoord(rq_df['RA (deg)'], rq_df['DEC (deg)'], unit=u.deg)
    columns = ['Name', 'RA (deg)', 'DEC (deg)', 'P(ms)', 'Edot',
               'Separation(deg)', 'pointing_name', 'utc_obs', 'output_path']
    fermi_rq_df = pd.DataFrame(columns=columns)

    fermi_rq_cnt = 0
    for i, pos in enumerate(fermi_rq_pos):
        if pos.separation(boresight_coords).deg <= survey_beam_radius*1.05:
            fermi_rq_df.loc[fermi_rq_cnt] = [rq_df['Name'][i].strip('PSR '),  pos.ra.deg, pos.dec.deg, rq_df['P (ms)'][i], rq_df['Edot'][i], pos.separation(
                boresight_coords).deg, pointing_name, utc_time, output_path]
            fermi_rq_cnt += 1

    return fermi_rq_df


def get_Fermi_association(opts, boresight_coords, pointing_name, utc_time, meta_output_path):
    """
    Write out the possible Fermi associations
    """
    # Read necessary info from FITS file
    if isinstance(opts.fits_file, type(None)):
        fgl4_fits = fits.open(
            '{}/4FGL_DR2_Ppsr.fits'.format(os.getcwd()))[1].data
    else:
        fgl4_fits = fits.open(opts.fits_file)[1].data

    fgl4_pos = SkyCoord(fgl4_fits.field('RAJ2000'),
                        fgl4_fits.field('DEJ2000'), unit=u.deg)
    fgl4_name = fgl4_fits.field('Source_Name')
    fgl4_val = fgl4_fits.field('P(psr)')
    fgl4_r95_semi_major = fgl4_fits.field('Conf_95_SemiMajor')
    fgl4_r95_semi_minor = fgl4_fits.field('Conf_95_SemiMinor')

    columns = ['Fermi Name', 'RA(deg)', 'DEC (deg)', 'P(psr)', 'r95_semi_major (deg)',
               'r95_semi_minor (deg)', 'Separation (deg)', 'pointing_name', 'utc_obs', 'output_path']
    fermi_source_df = pd.DataFrame(columns=columns)

    fermi_cnt = 0
    for i, pos in enumerate(fgl4_pos):
        if pos.separation(boresight_coords).deg <= survey_beam_radius*1.05:
            fermi_source_df.loc[fermi_cnt] = [fgl4_name[i], pos.ra.deg, pos.dec.deg, fgl4_val[i], fgl4_r95_semi_major[i],
                                              fgl4_r95_semi_minor[i], pos.separation(boresight_coords).deg, pointing_name, utc_time, meta_output_path]
            fermi_cnt += 1

    return fermi_source_df


def pointInEllipse(x, y, xp, yp, d, D, angle):
    """
    tests if a point[xp,yp] is within
    boundaries defined by the ellipse
    of center[x,y], diameter d D, and tilted at angle
    """

    cosa = math.cos(angle)
    sina = math.sin(angle)
    dd = d/2*d/2
    DD = D/2*D/2

    a = math.pow(cosa*(xp-x)+sina*(yp-y), 2)
    b = math.pow(sina*(xp-x)-cosa*(yp-y), 2)
    ellipse = (a/dd)+(b/DD)
    if ellipse <= 1:
        return True
    else:
        return False


def write_keep_beams_csv(opts, best_beams, best_psrs):
    """
    Write out a file with paths to filterbanks to be retained where known pulsars are expected 
    """
    with open(opts.meta, 'r') as f:
        all_info = json.load(f)

    filterbank_path = all_info['output_dir']
    columns = ['filterbank_path', 'username', 'reason']
    keep_beams_df = pd.DataFrame(columns=columns)

    for i, psr in enumerate(best_psrs):
        keep_beams_df.loc[i] = [filterbank_path+'/'+best_beams[i],
                                opts.username, opts.survey_name + ' KNOWN PSR ({})'.format(psr)]

    if keep_beams_df.shape[0] > 0:
        csv_full_path = "{}/{}_keep_beams_suggested.csv".format(
            opts.output_path, all_info['boresight'].split(',')[0])
        if os.path.isfile(csv_full_path):
            keep_beams_df.to_csv(csv_full_path, mode='a',
                                 header=False,  index=False)
        else:
            keep_beams_df.to_csv(csv_full_path, index=False)
    else:
        log.info("No beams suggested to be retained")


def generate_info_from_meta(opts):
    """
    Main function where all the info is generated from a meta file
    """
    with open(opts.meta, 'r') as f:
        all_info = json.load(f)

    # Initialise the plot
    plt.clf()
    fig, ax = plt.subplots()
    ax.set_aspect(1)

    # UTC Time
    utc_time = all_info['utc_start'].replace(" ", "T").replace("/", "-")
    time = Time(utc_time)

    # Pointing name
    pointing_name = all_info['boresight'].split(',')[0]

    # Output path
    meta_output_path = all_info['output_dir']

    # Beamshape parameters
    beam_shape = json.loads(all_info['beamshape'])
    # Convert angle to radians
    angle = np.radians(float(beam_shape['angle'])) - np.pi
    x_sigma = beam_shape['x']
    y_sigma = beam_shape['y']

    # Check if pointing name files already exist and skip them
    if os.path.isfile('{}/{}_{}.meta.png'.format(opts.output_path, pointing_name, utc_time)):
        log.info("Info about {} already exists in output path".format(pointing_name))
        return None

    log.info("Pointing Name: {}".format(pointing_name))
    log.info("Observed UTC Time: {}".format(utc_time))

    # Boresight
    boresight_ra = all_info['boresight'].split(',')[-2]
    boresight_dec = all_info['boresight'].split(',')[-1]
    boresight_coords = SkyCoord(
        frame='icrs', ra=boresight_ra, dec=boresight_dec, unit=(u.hour, u.deg))
    boresight_ra_deg = boresight_coords.ra.value
    boresight_dec_deg = boresight_coords.dec.value

    # Forcing the plot to use boresight as the centre
    bore = [boresight_ra, boresight_dec]
    bore_coords = SkyCoord(
        frame='icrs', ra=bore[0], dec=bore[1], unit=(u.hour, u.deg))
    bore_pixel_coordinates = convert_equatorial_coordinate_to_pixel(
        bore_coords, boresight_coords, time)
    bore_pixel_ra = boresight_ra_deg + bore_pixel_coordinates[0][0]
    bore_pixel_dec = boresight_dec_deg + bore_pixel_coordinates[0][1]
    ax.plot(bore_pixel_ra, bore_pixel_dec, '*', markersize=0.0)

    # Get all coherent beam coords
    coherent_beam_coords = get_coherent_beam_coords(opts.meta)

    # Beam shapes
    beam_width = 2.0*float(json.loads(all_info['beamshape'])['x'])
    beam_height = 2.0*float(json.loads(all_info['beamshape'])['y'])
    if utc_time < Time('2021-10-05T00:00:00'):
        beam_angle = float(json.loads(all_info['beamshape'])['angle'])
    else:
        # Sign change in FBFUSE tiling rollout
        beam_angle = -1.0*float(json.loads(all_info['beamshape'])['angle'])

    # Get all key value pairs for beams and sort them based on beam number
    vals = list(all_info['beams'].values())
    keys = list(all_info['beams'].keys())
    vals = [x for _, x in sorted(zip(keys, vals))]
    keys = sorted(keys)

    # Add known pulsar list based on ATNF or Pulsar survey scraper

    if opts.kp_catalogue == 'ATNF':
        log.info("Querying the ATNF catalogue and retrieving all known pulsars")
        columns = ['JNAME', 'RA(deg)', 'DEC (deg)', 'P0 (s)', 'DM', 'Closest beam(expected)',
                   'pointing_name', 'utc_obs', 'output_path', 'Within beam?', 'Neighbour beams']
        kp_df = pd.DataFrame(columns=columns)

        q = QueryATNF(params=['JNAME', 'RAJ', 'DECJ', 'P0', 'DM'], circular_boundary=(
            boresight_ra, boresight_dec, 2.0*incoherent_beam_radius))
        pulsar_list = np.array((q.table['JNAME', 'RAJ', 'DECJ', 'P0', 'DM']))

        log.info("{} known pulsars found within the incoherent beam".format(
            len(pulsar_list)))

        best_beams = []
        best_psrs = []
        for i, psr in enumerate(pulsar_list):
            psr_coords = SkyCoord(
                frame='icrs', ra=psr[1], dec=psr[2], unit=(u.hour, u.deg))
            psr_pixel_coordinates = convert_equatorial_coordinate_to_pixel(
                psr_coords, boresight_coords, time)
            psr_pixel_ra = boresight_ra_deg + psr_pixel_coordinates[0][0]
            psr_pixel_dec = boresight_dec_deg + psr_pixel_coordinates[0][1]
            psr_pixel_coords = SkyCoord(
                frame='icrs', ra=psr_pixel_ra, dec=psr_pixel_dec, unit=(u.hour, u.deg))
            ax.plot(psr_pixel_ra, psr_pixel_dec, '*',
                    label=psr[0]+' (ATNF)', markersize=7.5)

            # Check if pulsar within survey beam, gets closest beam + checks if psr within closest beam region + 3 closest beams
            if psr_coords.separation(boresight_coords).deg <= survey_beam_radius*1.05:
                log.info(
                    "{} expected within the survey beam region".format(psr[0]))
                pixel_beam_ras, pixel_beam_decs = get_pixel_coherent_beam_coordinates(
                    coherent_beam_coords, boresight_coords, time)

                # Sort beams by sensitivity fraction by modelling the beam as a 2D Gaussian
                pixel_beam_ras, pixel_beam_decs = get_pixel_coherent_beam_coordinates(
                    coherent_beam_coords, boresight_coords, time)
                G_values = []
                for i, pixel_ra in enumerate(pixel_beam_ras):
                    G_value = gaussian_ellipse_value(
                        psr_pixel_ra, psr_pixel_dec, pixel_ra, pixel_beam_decs[i], angle, x_sigma, y_sigma)
                    G_values.append(G_value)
                G_vals = np.asarray(G_vals, dtype=float)
                all_beams_sorted = np.argsort(-G_vals)

                best_ellipse = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[0]], pixel_beam_decs[all_beams_sorted[0]]),
                                       width=beam_width, height=beam_height, angle=beam_angle, edgecolor='blue', fc='grey', lw=1.5)
                ax.add_patch(best_ellipse)

                best_beam = 'cfbf00{:03d}'.format(all_beams_sorted[0])
                best_beams.append(best_beam)
                best_psrs.append(psr[0]+'; ATNF')

                # Use geometric formula directly since contains point function is failing.
                if pointInEllipse(pixel_beam_ras[all_beams_sorted[0]], pixel_beam_decs[all_beams_sorted[0]], psr_pixel_ra, psr_pixel_dec, beam_width, beam_height, beam_angle):
                    log.info("{} is within the {} beam region".format(
                        psr[0], best_beam))
                    within_flag = 'Y'
                else:
                    within_flag = 'N'

                neighbour_beam_list = ";".join(map(str, list(
                    all_beams_sorted[1:min(7, len(all_beams_sorted))])))  # Sort by just position


            else:
                best_beam = 'Outside survey beam'
                within_flag = 'N'
                neighbour_beam_list = "None"

            kp_df.loc[i] = [psr[0], psr_coords.ra.deg, psr_coords.dec.deg, psr[3], psr[4], best_beam,
                            pointing_name, utc_time, meta_output_path,  within_flag, neighbour_beam_list]

        kp_df.to_csv(
            '{}/{}_known_psrs.csv'.format(opts.output_path, pointing_name), index=False)

        if opts.keep_beams:
            log.info("Writing out a file with beams to keep")
            write_keep_beams_csv(opts, best_beams, best_psrs)

    elif opts.kp_catalogue == 'PSS':
        log.info(
            "Using the Pulsar survey scraper to retrieve known pulsars within incoherent beam")
        columns = ['JNAME', 'RA(deg)', 'DEC (deg)', 'P0 (s)', 'DM', 'Survey', 'Boresight separation (deg)',
                   'Closest beam(expected)', 'pointing_name', 'utc_obs', 'output_path', 'Within beam?', 'Neighbour beams']
        kp_df = pd.DataFrame(columns=columns)
        command = "python get_psrs_in_field.py --tag {}  --search_coordinates \"{} {}\" --search_radius {}".format(
            pointing_name, boresight_ra, boresight_dec, 2.0*incoherent_beam_radius)
        kp_out = subprocess.check_output(command, shell=True)
        if 'No pulsars' in str(kp_out):
            log.info("No pulsars from PSS in the incoherent beam region")
        else:
            best_beams = []
            best_psrs = []
            pss_data = pd.read_csv(
                '{}/{}_known_psrs.csv'.format(os.getcwd(), pointing_name))
            os.remove("{}/{}_known_psrs.csv".format(os.getcwd(), pointing_name))
            log.info("{} known pulsars found within the incoherent beam".format(
                len(pss_data.index)))
            for index, psr in pss_data.iterrows():
                psr_coords = SkyCoord(
                    frame='icrs', ra=psr['RA (deg)'], dec=psr['DEC (deg)'], unit=(u.deg, u.deg))
                psr_pixel_coordinates = convert_equatorial_coordinate_to_pixel(
                    psr_coords, boresight_coords, time)
                psr_pixel_ra = boresight_ra_deg + psr_pixel_coordinates[0][0]
                psr_pixel_dec = boresight_dec_deg + psr_pixel_coordinates[0][1]
                ax.plot(psr_pixel_ra, psr_pixel_dec, '*',
                        label=psr['PSR']+' '+'('+psr['Survey']+')', markersize=7.5)
                if 'HTRU-S' in psr['Survey']:
                    parkes_telescope_beam = Circle((psr_pixel_ra, psr_pixel_dec), 0.1166666, linestyle='--',
                                                   linewidth=2.5, fill=False, label='HTRU-S Low-latitude Parkes beam')
                    ax.add_patch(parkes_telescope_beam)

                if psr_coords.separation(boresight_coords).deg <= survey_beam_radius*1.05:
                    log.info(
                        "{} expected within the survey beam region".format(psr['PSR']))
                    pixel_beam_ras, pixel_beam_decs = get_pixel_coherent_beam_coordinates(
                        coherent_beam_coords, boresight_coords, time)

                    G_values = []
                    for i, pixel_ra in enumerate(pixel_beam_ras):
                        G_value = gaussian_ellipse_value(
                            psr_pixel_ra, psr_pixel_dec, pixel_ra, pixel_beam_decs[i], angle, x_sigma, y_sigma)
                        G_values.append(G_value)
                    G_vals = np.asarray(G_values, dtype=float)
                    all_beams_sorted = np.argsort(-G_vals)

                    best_ellipse = Ellipse(xy=(pixel_beam_ras[all_beams_sorted[0]], pixel_beam_decs[all_beams_sorted[0]]),
                                           width=beam_width, height=beam_height, angle=beam_angle, edgecolor='blue', fc='grey', lw=1.5)
                    ax.add_patch(best_ellipse)

                    best_beam = 'cfbf00{:03d}'.format(all_beams_sorted[0])
                    best_beams.append(best_beam)
                    best_psrs.append(psr['PSR']+'; '+psr['Survey'])

                    if pointInEllipse(pixel_beam_ras[all_beams_sorted[0]], pixel_beam_decs[all_beams_sorted[0]], psr_pixel_ra, psr_pixel_dec, beam_width, beam_height, beam_angle):
                        log.info("{} is within the {} beam region".format(
                            psr['PSR'], best_beam))
                        within_flag = 'Y'
                    else:
                        within_flag = 'N'

                    neighbour_beam_list = ";".join(map(str, list(
                        all_beams_sorted[1:min(7, len(all_beams_sorted))])))  # Sort by just position

                else:
                    best_beam = 'Outside survey beam'
                    within_flag = 'N'
                    neighbour_beam_list = "None"

                kp_df.loc[index] = [psr['PSR'], psr_coords.ra.deg, psr_coords.dec.deg,
                                    psr['P (ms)'], psr['DM (pc cm^-3)'], psr['Survey'], psr['Separation (deg)'], best_beam, pointing_name, utc_time, meta_output_path,  within_flag, neighbour_beam_list]

            kp_df.to_csv(
                '{}/{}_known_psrs.csv'.format(opts.output_path, pointing_name), index=False)

            if opts.keep_beams:
                log.info("Writing out a file with beams to keep")
                write_keep_beams_csv(opts, best_beams, best_psrs)

    # Output list of pulsars based on separate unpublished spreadsheets
    if opts.unpublished_flag:
        log.info("Checking unpublished spreadsheets for known pulsars..")
        if isinstance(opts.sheet_id, type(None)):
            raise Exception("Unique Sheet ID for spreadsheet not specified!")
        Columns = ['PSR', 'P(ms)', 'DM', 'Separation (deg)',
                   'pointing_name', 'utc_obs', 'output_path']
        unpublished_list = pd.DataFrame(columns=Columns)
        unpublished_df = pd.read_csv(
            'https://docs.google.com/spreadsheets/d/{}/gviz/tq?tqx=out:csv&sheet={}'.format(opts.sheet_id, opts.sheet_name))
        gls = unpublished_df['gl (deg) ']
        gbs = unpublished_df['gb (deg) ']
        unpublished_psr_coords = SkyCoord(
            gls*u.deg, gbs*u.deg, frame='galactic').transform_to('icrs')
        max_radius = survey_beam_radius + 0.1166  # 7 arcmin for L-band at Parkes
        unpublished_cnt = 0

        for i, unpublished_psr_coord in enumerate(unpublished_psr_coords):
            if boresight_coords.separation(unpublished_psr_coord).deg < max_radius:
                unpublished_list.loc[unpublished_cnt] = [unpublished_df['PSR Name '][i], unpublished_df['P(ms) '][i], unpublished_df['DM '][i], boresight_coords.separation(
                    unpublished_psr_coord).deg, pointing_name, utc_time, meta_output_path]
                log.info("{} potentially within survey beam".format(
                    unpublished_df['PSR Name '][i]))
                unpublished_psr_pixel_coordinates = convert_equatorial_coordinate_to_pixel(
                    unpublished_psr_coord, boresight_coords, time)
                unpublished_psr_pixel_ra = boresight_ra_deg + \
                    unpublished_psr_pixel_coordinates[0][0]
                unpublished_psr_pixel_dec = boresight_dec_deg + \
                    unpublished_psr_pixel_coordinates[0][1]
                ax.plot(unpublished_psr_pixel_ra, unpublished_psr_pixel_dec, '*',
                        label=unpublished_df['PSR Name '][i] + ' (HTRU unpublished)', markersize=7.5)
                telescope_beam = Circle((unpublished_psr_pixel_ra, unpublished_psr_pixel_dec),
                                        0.1166666, linestyle='--', linewidth=2.5, fill=False, label='HTRU Parkes beam')
                ax.add_patch(telescope_beam)
                unpublished_cnt += 1
        if unpublished_list.empty:
            log.info("No unpublished pulsars within the survey beam")
        else:
            unpublished_list.to_csv(
                '{}/{}_unpublished_known_psrs.csv'.format(opts.output_path, pointing_name), index=False)
            log.info("{} Unpublished psrs found and written to {}_unpublished_known_psrs.csv".format(
                len(unpublished_list.index), pointing_name))

    # Get Fermi sources in region and plot with r95 ellipse
    if opts.fermi_flag:
        log.info("checking for Fermi associations...")
        fermi_source_df = get_Fermi_association(
            opts, boresight_coords, pointing_name, utc_time, meta_output_path)
        if fermi_source_df.empty:
            log.info("No Fermi associations within survey beam region")
        else:
            fermi_source_df.to_csv(
                '{}/{}_Fermi_associations.csv'.format(opts.output_path, pointing_name), index=False)
            log.info("{} Fermi associations found and written to {}_Fermi_associations.csv".format(
                len(fermi_source_df.index), pointing_name))
        for index, row in fermi_source_df.iterrows():
            fermi_coords = SkyCoord(
                frame='icrs', ra=row[1], dec=row[2], unit=(u.deg, u.deg))
            pixel_fermi_coordinates = convert_equatorial_coordinate_to_pixel(
                fermi_coords, boresight_coords, time)
            pixel_fermi_ra = boresight_ra_deg + pixel_fermi_coordinates[0][0]
            pixel_fermi_dec = boresight_dec_deg + pixel_fermi_coordinates[0][1]
            ax.plot(pixel_fermi_ra, pixel_fermi_dec, '*',
                    label=row[0] + ' ({})'.format(row[3]), markersize=7.5)
            if not math.isnan(row[4]):  # If source is extended, skip r95 ellipse plot
                ellipse = Ellipse(xy=(pixel_fermi_ra, pixel_fermi_dec), width=2.0 *
                                  row[4], height=2.0*row[5], edgecolor='k', fc='none', lw=1.5, linestyle='--', label='Fermi r95 region')
                ax.add_patch(ellipse)
            else:
                log.info(
                    "{} is an extended Fermi source. R95 region is invalid".format(row[0]))

        # Get Fermi radio quiet sources in field.
        log.info("Checking for Fermi radio quiet pulsars...")
        fermi_rq_df = get_Fermi_rq_pulsars(
            opts, boresight_coords, pointing_name, utc_time, meta_output_path)
        if isinstance(fermi_rq_df, type(None)):
            log.info("No Radio quiet csv file found. Ignoring checking")
        elif fermi_rq_df.empty:
            log.info("No Fermi radio quiet pulsars within survey beam region")
        else:
            fermi_rq_df.to_csv(
                '{}/{}_Fermi_radio_quiet_pulsars.csv'.format(opts.output_path, pointing_name), index=False)
            log.info("{} Fermi radio quiet pulsars found and written to {}_Fermi_radio_quiet_pulsars.csv".format(
                len(fermi_rq_df.index), pointing_name))
            pixel_beam_ras, pixel_beam_decs = get_pixel_coherent_beam_coordinates(
                coherent_beam_coords, boresight_coords, time)
            best_psrs = []
            best_beams = []
            for index, row in fermi_rq_df.iterrows():
                fermi_coords = SkyCoord(
                    frame='icrs', ra=row[1], dec=row[2], unit=(u.deg, u.deg))
                pixel_fermi_coordinates = convert_equatorial_coordinate_to_pixel(
                    fermi_coords, boresight_coords, time)
                pixel_fermi_ra = boresight_ra_deg + \
                    pixel_fermi_coordinates[0][0]
                pixel_fermi_dec = boresight_dec_deg + \
                    pixel_fermi_coordinates[0][1]
                ax.plot(pixel_fermi_ra, pixel_fermi_dec, '*',
                        label=row[0]+'(q)', markersize=7.5)
                psr_idx, psr_d2d, psr_d3d = fermi_coords.match_to_catalog_sky(
                    coherent_beam_coords)
                best_ellipse = Ellipse(xy=(pixel_beam_ras[psr_idx], pixel_beam_decs[psr_idx]), width=beam_width,
                                       height=beam_height, angle=beam_angle, edgecolor='blue', fc='grey', lw=1.5)
                ax.add_patch(best_ellipse)
                best_beam = 'cfbf00{:03d}'.format(psr_idx)
                best_beams.append(best_beam)
                best_psrs.append(row[0]+'; Fermi q')

            if opts.keep_beams:
                log.info("Appending/writing to a file with beams to keep")
                write_keep_beams_csv(opts, best_beams, best_psrs)

    # Add ellipses
    for i in range(len(vals)):
        if 'unset' in vals[i]:
            continue
        beam_ra = vals[i].split(',')[-2]
        beam_dec = vals[i].split(',')[-1]
        beam_no = keys[i][-3:]

        # Convert equatorial beam coordinates to pixel coordinates
        beam_coords = SkyCoord(frame='icrs', ra=beam_ra,
                               dec=beam_dec, unit=(u.hour, u.deg))
        pixel_coordinates = convert_equatorial_coordinate_to_pixel(
            beam_coords, boresight_coords, time)
        pixel_beam_ra = boresight_ra_deg + pixel_coordinates[0][0]
        pixel_beam_dec = boresight_dec_deg + pixel_coordinates[0][1]

        ellipse = Ellipse(xy=(pixel_beam_ra, pixel_beam_dec), width=beam_width,
                          height=beam_height, angle=beam_angle, edgecolor='blue', fc='none', lw=1.5)
        ax.add_patch(ellipse)
        # Add beam numbers
        ax.annotate(beam_no, (pixel_beam_ra, pixel_beam_dec), fontsize=2)

    # Add user specified coordinates
    if not isinstance(opts.user_coords, type(None)):
        user_coords = SkyCoord(frame='icrs', ra=opts.user_coords.split(
            ' ')[0], dec=opts.user_coords.split(' ')[1], unit=(u.hour, u.deg))
        user_pixel_coordinates = convert_equatorial_coordinate_to_pixel(
            user_coords, boresight_coords, time)
        user_pixel_ra = boresight_ra_deg + user_pixel_coordinates[0][0]
        user_pixel_dec = boresight_dec_deg + user_pixel_coordinates[0][1]
        ax.plot(user_pixel_ra, user_pixel_dec, '*',
                label='User specified', markersize=7.5)

    # Add user beam radius
    if opts.beam_radius == survey_beam_radius and opts.plot_beams == 1:
        user_circle = Circle((boresight_ra_deg, boresight_dec_deg), opts.beam_radius,
                             color='red', linestyle='--', linewidth=2.5, fill=False, label='Survey beam')
        ax.add_patch(user_circle)
    if opts.beam_radius != survey_beam_radius:
        user_circle = Circle((boresight_ra_deg, boresight_dec_deg), opts.beam_radius,
                             color='red', linestyle='--', linewidth=2.5, fill=False, label='Telescope beam')
        ax.add_patch(user_circle)

    # Incoherent beam radius
    if opts.plot_beams:
        incoherent_circle = Circle((boresight_ra_deg, boresight_dec_deg), incoherent_beam_radius,
                                   color='green', linestyle='--', linewidth=2.5, fill=False, label='Incoherent beam')
        ax.add_patch(incoherent_circle)

    # Get elevation

    # Location on sky
    pointing_coords = SkyCoord(
        frame='icrs', ra=boresight_ra, dec=boresight_dec, unit=(u.hour, u.deg))

    # Calculate alt/elevation
    coord_altaz = pointing_coords.transform_to(
        AltAz(obstime=time, location=meerkat))
    elv_value = coord_altaz.alt.deg

    # Survey beam filling ratio - Accounts for min. of half PSR power
    no_of_beams = nbeams(all_info)
    log.info("{} coherent beams tile the survey beam".format(no_of_beams))
    total_area = np.pi*0.25*beam_width*beam_height*no_of_beams
    survey_beam_fill_factor = (total_area/survey_beam_area)/0.91

    # plotting ornaments
    ax.set_xlabel('Right Ascension (Degrees)')
    ax.set_ylabel('Declination (Degrees)')
    ax.set_title("Pointing: {}, UTC: {}, Elevation: {}".format(
        pointing_name, utc_time, elv_value))
    if opts.plot_beams:
        # ax.set_xlim(boresight_coords.ra.deg - incoherent_beam_radius, boresight_coords.ra.deg + incoherent_beam_radius)
        ax.set_xlim(boresight_coords.ra.deg - survey_beam_radius,
                    boresight_coords.ra.deg + survey_beam_radius)
        # ax.set_ylim(boresight_coords.dec.deg - incoherent_beam_radius, boresight_coords.dec.deg + incoherent_beam_radius)
        ax.set_ylim(boresight_coords.dec.deg - survey_beam_radius,
                    boresight_coords.dec.deg + survey_beam_radius)
    plt.legend(prop={"size": 6})
    plt.savefig('{}/{}_{}.meta.png'.format(opts.output_path,
                pointing_name, utc_time), dpi=400)

    log.info("Output path: {}".format(opts.output_path))


if __name__ == "__main__":
    # Select options
    parser = optparse.OptionParser()
    parser.add_option('--meta_path', type=str,
                      help='Path to meta file (Full path)', dest='meta')
    parser.add_option('--check_unpublished', type=int,
                      help='Flag for checking unpublished sources (Default: 0)', dest='unpublished_flag', default=0)
    parser.add_option('--sheet_id', type=str,
                      help='Unique spreadheet ID', dest='sheet_id', default=None)
    parser.add_option('--sheet_name', type=str, help='Name of sheet of unpublished pulsars (Default: reprocessing_discoveries)',
                      dest='sheet_name', default='reprocessing_discoveries')
    parser.add_option('--check_fermi', type=int,
                      help='Check for possible Fermi associations (Default: 1)', dest='fermi_flag', default=1)
    parser.add_option('--fits_file', type=str,
                      help='Full path for Fermi fits file (Defaults to Fermi fits file in current working directory)', dest='fits_file', default=None)
    parser.add_option('--known_pulsar', type=str,
                      help='Cross match with known pulsar list as specified. Options: ATNF, PSS', dest='kp_catalogue', default='ATNF')
    parser.add_option('--user_coordinate', type=str,
                      help=' Plot user specified coordinates  e.g. 12:08 -59:36', dest='user_coords', default=None)
    parser.add_option('--user_beam_radius', type=float, help='Plot user specified beam radius in degrees (Defaults to survey beam radius for MGPS-L)',
                      dest='beam_radius', default=survey_beam_radius)
    parser.add_option('--username', type=str, help='Username of person running (Defaults to username of local machine)',
                      dest='username', default=getpass.getuser())
    parser.add_option('--survey_name', type=str,
                      help='Survey name (e.g. TRAPUM-Fermi, MGPS-L) ; Default is MGPS-L', dest='survey_name', default='MGPS-L')
    parser.add_option('--output_path', type=str, help='Path to store output files (Defaults to current working directory)',
                      dest='output_path', default=os.getcwd())
    parser.add_option('--keep_beams', type=int,
                      help='Flag to write out list of beams to keep based on expected known pulsars (Defaults to 1)', dest='keep_beams', default=1)
    parser.add_option('--plot_beams', type=int,
                      help='Flag to indicate if survey and incoherent beams should be plotted (Defaults to 1)', dest='plot_beams', default=1)
    opts, args = parser.parse_args()

    # Expand user for all path arguments
    opts.output_path = os.path.expanduser(opts.output_path)
    opts.meta = os.path.expanduser(opts.meta)
    if not isinstance(opts.fits_file, type(None)):
        opts.meta = os.path.expanduser(opts.meta)

    # Generate all info for pointing
    generate_info_from_meta(opts)
