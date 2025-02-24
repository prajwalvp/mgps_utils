import os
import sys
import glob
import math
import shutil
import subprocess
import tarfile
import optparse
import logging
import numpy as np
import pandas as pd
from ast import literal_eval
from astropy.coordinates import SkyCoord
import astropy.units as u 
import getpass

"""
Updated version of candy_scraper.py which writes out beam lists with columns: filterbank_path, username, reason
"""

log = logging.getLogger('candy_scraper2')
FORMAT = "[%(levelname)s - %(asctime)s - %(filename)s:%(lineno)s] %(message)s"
logging.basicConfig(format=FORMAT)
log.setLevel('INFO')


def make_tarfile(output_path, input_path, name):
    with tarfile.open(output_path + '/' + name, "w:gz") as tar:
        tar.add(input_path, arcname=os.path.basename(input_path))

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


def calculate_ellipse_coefficients(angle, x_sigma, y_sigma):
    """
    Calculate the coefficients (a, b, c) for the quadratic form of an ellipse.

    Args:
        angle (float): The orientation of the ellipse in radians.
        x_sigma (float): The semi-major axis of the ellipse in radians.
        y_sigma (float): The semi-minor axis of the ellipse in radians.

    Returns:
        tuple: Coefficients (a, b, c) for the quadratic form.
    """
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)
    x_sigma_sq = x_sigma * x_sigma
    y_sigma_sq = y_sigma * y_sigma

    a = (cos_angle * cos_angle) / (2 * x_sigma_sq) + (sin_angle * sin_angle) / (2 * y_sigma_sq)
    b = (-math.sin(2 * angle)) / (4 * x_sigma_sq) + (math.sin(2 * angle)) / (4 * y_sigma_sq)
    c = (sin_angle * sin_angle) / (2 * x_sigma_sq) + (cos_angle * cos_angle) / (2 * y_sigma_sq)

    return a, b, c


def get_beam_intersection(b1_coords, b2_coords, angle, x_sigma, y_sigma):
    """
    Calculate the intersection value between two beams (b1 and b2).

    Args:
        b1_coords (SkyCoord): Coordinates of the first beam in radians.
        b2_coords (SkyCoord): Coordinates of the second beam in radians.
        angle (float): The orientation of the ellipse in radians.
        x_sigma (float): The semi-major axis of the ellipse in radians.
        y_sigma (float): The semi-minor axis of the ellipse in radians.

    Returns:
        float: The intersection value between the two beams.
    """
    # Convert coordinates to radians
    x_mean = b1_coords.ra.radian
    y_mean = b1_coords.dec.radian
    x = b2_coords.ra.radian
    y = b2_coords.dec.radian

    # Calculate coefficients
    a, b, c = calculate_ellipse_coefficients(angle, x_sigma, y_sigma)

    # Compute the quadratic form
    dx = x - x_mean
    dy = y - y_mean
    d = a * dx * dx + 2 * b * dx * dy + c * dy * dy

    # Safeguard against overflow
    if d > 700:  # exp(700) is approximately 1e304, close to the limit for 64-bit floats
        return 0.0  # Return 0 for very large d (beams are effectively non-overlapping)
    else:
        return math.exp(d)


def get_neighbour_beams(opts, df):
    """
    Select up to 6 neighboring beams w.r.t. reference beam based on the intersection value.
    """ 
    columns = ['filterbank_path', 'username', 'reason']
    neighbour_set = pd.DataFrame(columns=columns)

    for i, row in df.iterrows():
        meta_file = opts.main_dir + '/' + row['metafile_path']
        coherent_beam_coords = get_coherent_beam_coords(meta_file)
        ref_coords = coherent_beam_coords[int(row['beam_name'].strip('cfbf'))]

        # Extract ellipse configuration (angle, x_sigma, y_sigma) from the meta file
        with open(meta_file, 'r') as f:
            meta_data = literal_eval(f.read())
        beam_shape = literal_eval(meta_data['beamshape'])
        angle = np.radians(beam_shape['angle']) - np.pi  # Convert angle to radians
        x_sigma = beam_shape['x']   
        y_sigma = beam_shape['y']

        # Step 1: Track original indices
        original_indices = list(range(len(coherent_beam_coords)))

        # Step 2: Exclude the reference beam
        filtered_beam_coords = []
        filtered_indices = []
        for i, beam_coord in enumerate(coherent_beam_coords):
            if beam_coord != ref_coords:
                filtered_beam_coords.append(beam_coord)
                filtered_indices.append(original_indices[i])

        # Step 3: Calculate intersection values
        intersection_values = []
        for beam_coord in filtered_beam_coords:
            intersection_value = get_beam_intersection(ref_coords, beam_coord, angle, x_sigma, y_sigma)
            intersection_values.append(intersection_value)

        # Step 4: Sort beams by intersection value and select the closest ones
        sorted_indices = np.argsort(intersection_values)
        neighbour_beam_list = [filtered_indices[i] for i in sorted_indices[:6]]  # Select the 6 closest beams


        # Add neighboring beams to the DataFrame
        for beam_num in neighbour_beam_list:
            filterbank_path = os.path.dirname(row['filterbank_path']) + '/' + 'cfbf00{:03d}'.format(int(beam_num))
            neighbour_set.loc[len(neighbour_set)] = [filterbank_path, opts.username, opts.survey_name + ' T1_CAND n']

    return neighbour_set

def write_t1_t2_beams(opts):
    """
    Write out T1, T2 beams 
    """ 
    columns = ['filterbank_path','username','reason']
    t1_t2_all = pd.DataFrame(columns=columns)

    #classified_files =  [c for c in glob.glob("{}/**/*{}*.csv".format(opts.main_dir, opts.tag)) if 'autosave' not in c]
    classified_files =  [c for c in glob.glob("{}/*{}*.csv".format(opts.main_dir, opts.tag)) if 'autosave' not in c if 'keep' not in c if 'full' not in c]
    #all_candidate_csvs = glob.glob("{}/**/candidates.csv".format(opts.main_dir))
    all_candidate_csvs = glob.glob("{}/candidates.csv".format(opts.main_dir))

    if len(classified_files) != len(all_candidate_csvs):
        raise Exception("The total number of files with input tag does not match the total candidate.csv files! Please check for multiple redundant files")

    # If this expection is not raised, I assume that there is just one candidates.csv file and a corresponding classified file
    classified_file = classified_files[0]
    candidate_csv = all_candidate_csvs[0]


    # Check if all candidates have been classified
    if (pd.read_csv(classified_file).shape[0] > pd.read_csv(candidate_csv).shape[0]):
        log.warning("The number of entries in the classification file is not equal to the number of entries in the candidates.csv. You have classified the candidates multiple times!")
    elif (pd.read_csv(classified_file).shape[0] < pd.read_csv(candidate_csv).shape[0]):
        log.warning("The number of entries in the classification file is not equal to the number of entries in the candidates.csv. All candidates were not classified!")
    
    
    #for i,classified_file in enumerate(classified_files):
        #print (i, classified_file)
    log.info('Checking {}'.format(os.path.dirname(classified_file)))
    df = pd.read_csv(classified_file)
    t1_t2_df =  (df[(df['classification']=='T1_CAND') | (df['classification']=='T2_CAND')])
        
    if t1_t2_df.empty:
        log.info('No T1/T2 files in {}. Returning empty T1/T2 list'.format(os.path.dirname(classified_file)))         
        return t1_t2_all

    else:
        candidate_meta_df = pd.read_csv(candidate_csv)
        t1_t2_beams_df = candidate_meta_df[candidate_meta_df['png_path'].isin(t1_t2_df['png'])][['filterbank_path', 'png_path']].reset_index(drop=True)
        log.info("Number of T1/T2 cands found: {}".format(t1_t2_beams_df.drop_duplicates().shape[0]))
        t1_t2_beams_df['filterbank_path'] = t1_t2_beams_df['filterbank_path'].apply(lambda x:x.replace(x,os.path.dirname(x)))             
        t1_t2_beams_df = t1_t2_beams_df.assign(username = opts.username)    
        t1_t2_beams_df_renamed = t1_t2_beams_df.rename({'png_path':'png'}, axis=1)
        merged_df = pd.merge(t1_t2_beams_df_renamed, t1_t2_df, on='png')
        merged_df = merged_df.assign(reason = opts.survey_name + ' ' + merged_df['classification'])
        t1_t2_all = pd.concat([t1_t2_all, merged_df[['filterbank_path', 'username', 'reason']]], ignore_index=True) 



        # Add neighbour beams by default if candidate is T1
        if opts.keep_neighbours == 1: 
            if t1_t2_df['classification'].str.contains('T1_CAND').sum():
                log.info("T1 candidate(s) found. Neighbouring beams will be added by default")     
                t1_df = df[df['classification']=='T1_CAND'] 
                t1_beams_df = candidate_meta_df[candidate_meta_df['png_path'].isin(t1_df['png'])][['filterbank_path','beam_name', 'metafile_path']]  
                t1_beams_df['filterbank_path'] = t1_beams_df['filterbank_path'].apply(lambda x:x.replace(x,os.path.dirname(x)))
                t1_beams_df.to_csv('t1_beams')
                t1_neighbour_set_df = get_neighbour_beams(opts, t1_beams_df)
                t1_t2_all = pd.concat([t1_t2_all, t1_neighbour_set_df], ignore_index=True) 
            else:
                log.info("No T1 candidates to search neighbours")


    return t1_t2_all 


def write_known_pulsar_beams(opts):
    """
    Write out known pulsar beams 
    """ 
    columns = ['filterbank_path','username','reason']
    kp_all = pd.DataFrame(columns=columns)

    #classified_files =  [c for c in glob.glob("{}/**/*{}*.csv".format(opts.main_dir, opts.tag)) if 'autosave' not in c]
    classified_files =  [c for c in glob.glob("{}/*{}*.csv".format(opts.main_dir, opts.tag)) if 'autosave' not in c if 'keep' not in c if 'full' not in c]
    #all_candidate_csvs = glob.glob("{}/**/candidates.csv".format(opts.main_dir))
    all_candidate_csvs = glob.glob("{}/candidates.csv".format(opts.main_dir))

    if len(classified_files) != len(all_candidate_csvs):
        raise Exception("The total number of classification files does not match the total candidate.csv files! Please check for multiple redundant files")
    
    for i,classified_file in enumerate(classified_files):
        log.info('Checking {}'.format(os.path.dirname(classified_file)))
        df = pd.read_csv(classified_file)
        kp_df =  df[df['classification']=='KNOWN_PSR'].drop_duplicates()
        
        if kp_df.empty:
            log.info('No known pulsar files in {}'.format(os.path.dirname(classified_file)))
            continue         
        else:
            candidate_meta_df = pd.read_csv(all_candidate_csvs[i])
            kp_beams_df = candidate_meta_df[candidate_meta_df['png_path'].isin(kp_df['png'])][['filterbank_path']]
            kp_beams_df['filterbank_path'] = kp_beams_df['filterbank_path'].apply(lambda x:x.replace(x,os.path.dirname(x)))

            kp_beams_df = kp_beams_df.assign(username=opts.username)
            
            kp_beams_df['reason'] = kp_df[kp_df['png'].isin(candidate_meta_df['png_path'])]['classification'].to_list()
            kp_beams_df['reason'] = opts.survey_name + ' ' + kp_beams_df['reason'].astype(str) 

            log.info("Number of known pulsar cands found: {}".format(kp_beams_df.drop_duplicates().shape[0]))
            #kp_all = kp_all.append(kp_beams_df, ignore_index=True)
            kp_all = pd.concat([kp_all, kp_beams_df], ignore_index=True)


    return kp_all 


def write_out_second_revision_tar_file(opts):
    """
    Make and tar up a CandyJar compatible directory made of T1/T2 cands for multiple viewings 
    """
    columns1= ['pointing_id','beam_id','beam_name','source_name',
              'ra','dec','gl','gb',
              'mjd_start','utc_start',
              'f0_user','f0_opt','f0_opt_err','f1_user','f1_opt','f1_opt_err',
              'acc_user', 'acc_opt', 'acc_opt_err', 
              'dm_user', 'dm_opt', 'dm_opt_err',
              'sn_fft', 'sn_fold',
              'pepoch',
              'maxdm_ymw16', 'dist_ymw16',
              'pics_trapum_ter5', 'pics_palfa',
              'png_path', 'metafile_path', 'filterbank_path', 'candidate_tarball_path']

    columns2 = ['utc','png','classification']

    t1_t2_all_meta_df = pd.DataFrame(columns=columns1)
    t1_t2_all_classified_df = pd.DataFrame(columns=columns2)

    # Make all necessary directories
    second_round_path = opts.output_dir + '/{}_{}_round1'.format(os.path.basename(os.path.normpath(opts.main_dir)),opts.tag)
    try:
        os.makedirs(second_round_path+'/plots')
        os.makedirs(second_round_path+'/metafiles')
        
    except FileExistsError:       
        shutil.rmtree(second_round_path)
        os.makedirs(second_round_path+'/plots') 
        os.makedirs(second_round_path+'/metafiles') 

    # Transfer all necessary  T1/T2 plots + metafiles
    #classified_files =  [c for c in glob.glob("{}/**/*{}*.csv".format(opts.main_dir, opts.tag)) if 'autosave' not in c]
    classified_files =  [c for c in glob.glob("{}/*{}*.csv".format(opts.main_dir, opts.tag)) if 'autosave' not in c if 'keep' not in c if 'full' not in c]
    #all_candidate_csvs = glob.glob("{}/**/candidates.csv".format(opts.main_dir)) 
    all_candidate_csvs = glob.glob("{}/candidates.csv".format(opts.main_dir)) 

    if len(classified_files) != len(all_candidate_csvs):
        raise Exception("The total number of files with input tag does not match the total candidate.csv files! Please check for multiple redundant files")

    for i,classified_file in enumerate(classified_files):
        log.info('Checking {}'.format(os.path.dirname(classified_file)))
        df = pd.read_csv(classified_file)
        t1_t2_df =  (df[(df['classification']=='T1_CAND') | (df['classification']=='T2_CAND')])
        if t1_t2_df.empty:
            continue
        else:
            candidate_meta_df = pd.read_csv(all_candidate_csvs[i])
    
            t1_t2_meta_df = candidate_meta_df[candidate_meta_df['png_path'].isin(t1_t2_df['png'])]       
            t1_t2_classified_df = t1_t2_df[t1_t2_df['png'].isin(candidate_meta_df['png_path'])]
            # Copy over metafile and png files to respective directories
            log.info("Copying over meta and png files to {}".format(second_round_path))
            for val in t1_t2_meta_df['metafile_path'].values:            
                shutil.copyfile(os.path.dirname(classified_file)+'/'+ val, 
                                            opts.output_dir+'/{}_{}_round1/metafiles/{}'.format(os.path.basename(os.path.normpath(opts.main_dir)), opts.tag, os.path.basename(val)))       
            for val in t1_t2_meta_df['png_path'].values:
                shutil.copyfile(os.path.dirname(classified_file)+'/'+val,                                
                                opts.output_dir+'/{}_{}_round1/plots/{}'.format(os.path.basename(os.path.normpath(opts.main_dir)), opts.tag,                                                 os.path.basename(val)))
            
            #t1_t2_all_classified_df = t1_t2_all_classified_df.append(t1_t2_classified_df, ignore_index=True)            
            t1_t2_all_classified_df = pd.concat([t1_t2_all_classified_df, t1_t2_classified_df], ignore_index=True)            
            #t1_t2_all_meta_df = t1_t2_all_meta_df.append(t1_t2_meta_df, ignore_index=True)
            t1_t2_all_meta_df = pd.concat([t1_t2_all_meta_df, t1_t2_meta_df], ignore_index=True)
          
    # Write out candidates.csv and first classification csv file
    log.info("Writing out candidates.csv file in {}".format(second_round_path))
    t1_t2_all_meta_df.to_csv(second_round_path+'/candidates.csv', index=False)
    log.info("Writing out first classification file in {}".format(second_round_path))
    t1_t2_all_classified_df.to_csv(second_round_path+'/{}_classified1.csv'.format(opts.tag), index=False)
    
    # Tar up the directory
    tar_name = "{}_{}_round1_T1_T2.tar.gz".format(os.path.basename(os.path.normpath(opts.main_dir)), opts.tag)
    make_tarfile(opts.output_dir, second_round_path, tar_name)    

    
    
if __name__ == "__main__":
    # Select options
    parser = optparse.OptionParser()
    parser.add_option('--name_tag', type=str, help = 'Name tag for searching classified files', dest='tag', default='autosave')
    parser.add_option('--username', type=str, help = 'MGPS/TRAPUM login username; Default is local machine username', dest='username', default=getpass.getuser())
    parser.add_option('--survey_name', type=str, help = 'Survey name (e.g. MMGPS-S, MMGPS-U) ; Default is MMGPS-S', dest='survey_name',default='MMGPS-S')
    parser.add_option('--main_dir', type=str, help='Root directory where all pointing subdirectories are stored', dest='main_dir')
    parser.add_option('--output_dir',type=str, help='Output directory where csvs/tar files will be written out (Default is same as root directory specified)', dest='output_dir')
    parser.add_option('--separate_csvs',type=int, help='Flag for writing out separate csvs for T1/T2 and pulsars (Default:0)', dest='sep_csv',default=0)
    parser.add_option('--second_revision',type=int, help='Flag for writing out a tarball of T1/T2 plots for second revision (Default:1)', dest='second_revision',default=1)
    parser.add_option('--keep_neighbours',type=int, help='Flag for keeping neighbouring beams w.r.t T1 candidates (Default:0)', dest='keep_neighbours',default=0)
    #parser.add_option('--known_psr_filter',type=str, help='Filter the non boresight beams to ensure best SNR/position beam retained', dest='filter_psr)
    opts, args = parser.parse_args()



    if opts.main_dir is None:
        log.info("Pointing directory not specified!")
        sys.exit(0) 
    if opts.tag is None:
        log.info("Tag to find classified files not specified!")
        sys.exit(0) 

    if opts.output_dir is None:
        log.debug("Output directory not specified")
        log.info("Files will be written in main specified directory")
        opts.output_dir = os.path.normpath(opts.main_dir)

    # Check for tilda in main dir and output dir
    if '~' in opts.output_dir:
        opts.output_dir = os.path.expanduser(opts.output_dir) 
    if '~' in opts.main_dir:
        opts.main_dir = os.path.expanduser(opts.main_dir)
 
    t1_t2_df = write_t1_t2_beams(opts)
    if opts.sep_csv:
        # Write out T1,T2 cand beam path+names
        log.info("Writing out beams for {} T1/T2 cands in total found in {}".format(t1_t2_df.shape[0], opts.main_dir))
        t1_t2_df.to_csv('{}/{}_{}_keep_T1_T2_cands.csv'.format(opts.output_dir, os.path.basename(os.path.normpath(opts.main_dir)), opts.tag), index = False)
    
    
    kp_df = write_known_pulsar_beams(opts) 
    if opts.sep_csv:
        # Write out known pulsar beam path+names
        log.info("Writing out beams for {} known pulsar cands in total found in {}".format(kp_df.drop_duplicates().shape[0], opts.main_dir))
        kp_df.drop_duplicates().to_csv('{}/{}_{}_keep_pulsar_cands.csv'.format(opts.output_dir, os.path.basename(os.path.normpath(opts.main_dir)), opts.tag), index = False)

    # Write out a csv file for beams to be retained
    #all_beams_df = t1_t2_df.append(kp_df, ignore_index=True)
    all_beams_df = pd.concat([t1_t2_df,kp_df], ignore_index=True)
    log.info("Writing out beams for all known pulsar and T1/T2 cands (Total: {}) found in {}".format(all_beams_df.drop_duplicates().shape[0], opts.main_dir))
    all_beams_df.drop_duplicates().to_csv('{}/{}_{}_keep_beams.csv'.format(opts.output_dir, os.path.basename(os.path.normpath(opts.main_dir)), opts.tag), index = False)
  
    # Write out a second revision file
    if opts.second_revision:
        log.info("Writing out a tarred file containing T1/T2 candidates for second revision")
        write_out_second_revision_tar_file(opts) 
 
