import os
import path
import sys
import glob
import shutil
import tarfile
import optparse
import logging
import pandas as pd


log = logging.getLogger('viewed_candidate_organiser')
FORMAT = "[%(levelname)s - %(asctime)s - %(filename)s:%(lineno)s] %(message)s"
logging.basicConfig(format=FORMAT)
log.setLevel('INFO')



def make_tarfile(output_path, input_path, name):
    with tarfile.open(output_path + '/' + name, "w:gz") as tar:
        tar.add(input_path, arcname=os.path.basename(input_path))


def write_t1_t2_beams(opts):
    """
    Write out T1, T2 beams 
    """ 
    columns = ['filterbank_path','beam_name']
    t1_t2_all = pd.DataFrame(columns=columns)

    classified_files =  [c for c in glob.glob("{}/**/*{}*.csv".format(opts.main_dir, opts.tag)) if 'autosave' not in c]
    all_candidate_csvs = glob.glob("{}/**/candidates.csv".format(opts.main_dir))

    if len(classified_files) != len(all_candidate_csvs):
        raise Exception("The total number of classification files does not match the total candidate.csv files! Please check for multiple redundant files")
    
    for i,classified_file in enumerate(classified_files):
        log.info('Checking {}'.format(os.path.dirname(classified_file)))
        df = pd.read_csv(classified_file)
        t1_t2_df =  (df[(df['classification']=='T1_CAND') | (df['classification']=='T2_CAND')])
        
        if t1_t2_df.empty:
            log.info('No T1/T2 files in {}'.format(os.path.dirname(classified_file)))
            continue         
        else:
            candidate_meta_df = pd.read_csv(all_candidate_csvs[i])
            t1_t2_beams_df = candidate_meta_df[candidate_meta_df['png_path'].isin(t1_t2_df['png'])][['filterbank_path','beam_name']]
            t1_t2_beams_df['filterbank_path'] = t1_t2_beams_df['filterbank_path'].apply(lambda x:x.replace(os.path.join(os.path.basename(os.path.dirname(x)), os.path.basename(x)), '%'))
            log.info("Number of T1/T2 cands found: {}".format(t1_t2_beams_df.shape[0]))
            t1_t2_all = t1_t2_all.append(t1_t2_beams_df, ignore_index=True)


    return t1_t2_all 


def write_known_pulsar_beams(opts):
    """
    Write out known pulsar beams 
    """ 
    columns = ['filterbank_path','beam_name']
    kp_all = pd.DataFrame(columns=columns)

    classified_files =  [c for c in glob.glob("{}/**/*{}*.csv".format(opts.main_dir, opts.tag)) if 'autosave' not in c]
    all_candidate_csvs = glob.glob("{}/**/candidates.csv".format(opts.main_dir))

    if len(classified_files) != len(all_candidate_csvs):
        raise Exception("The total number of classification files does not match the total candidate.csv files! Please check for multiple redundant files")
    
    for i,classified_file in enumerate(classified_files):
        log.info('Checking {}'.format(os.path.dirname(classified_file)))
        df = pd.read_csv(classified_file)
        kp_df =  df[df['classification']=='KNOWN_PSR']
        
        if kp_df.empty:
            log.info('No known pulsar files in {}'.format(os.path.dirname(classified_file)))
            continue         
        else:
            candidate_meta_df = pd.read_csv(all_candidate_csvs[i])
            kp_beams_df = candidate_meta_df[candidate_meta_df['png_path'].isin(kp_df['png'])][['filterbank_path','beam_name']]
            kp_beams_df['filterbank_path'] = kp_beams_df['filterbank_path'].apply(lambda x:x.replace(os.path.join(os.path.basename(os.path.dirname(x)), os.path.basename(x)), '%'))
            log.info("Number of known pulsar cands found: {}".format(kp_beams_df.shape[0]))
            kp_all = kp_all.append(kp_beams_df, ignore_index=True)


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
    second_round_path = opts.output_dir + '/{}_{}_round1'.format(os.path.basename(opts.main_dir),opts.tag)
    try:
        os.makedirs(second_round_path+'/plots')
        os.makedirs(second_round_path+'/metafiles')
        
    except FileExistsError:       
        shutil.rmtree(second_round_path)
        os.makedirs(second_round_path+'/plots') 
        os.makedirs(second_round_path+'/metafiles') 

    # Transfer all necessary  T1/T2 plots + metafiles
    classified_files =  [c for c in glob.glob("{}/**/*{}*.csv".format(opts.main_dir, opts.tag)) if 'autosave' not in c]
    all_candidate_csvs = glob.glob("{}/**/candidates.csv".format(opts.main_dir)) 

    if len(classified_files) != len(all_candidate_csvs):
        raise Exception("The total number of classification files does not match the total candidate.csv files! Please check for multiple redundant files")

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
            shutil.copyfile(os.path.dirname(classified_file)+'/'+t1_t2_meta_df['metafile_path'].values[0], 
                                            opts.output_dir+'/{}_{}_round1/metafiles/{}'.format(os.path.basename(opts.main_dir), opts.tag,                                              os.path.basename(t1_t2_meta_df['metafile_path'].values[0])))   
    
            for val in t1_t2_meta_df['png_path'].values:
                shutil.copyfile(os.path.dirname(classified_file)+'/'+val,                                
                                opts.output_dir+'/{}_{}_round1/plots/{}'.format(os.path.basename(opts.main_dir), opts.tag,                                                 os.path.basename(val)))
            
            t1_t2_all_classified_df = t1_t2_all_classified_df.append(t1_t2_classified_df, ignore_index=True)            
            t1_t2_all_meta_df = t1_t2_all_meta_df.append(t1_t2_meta_df, ignore_index=True)
          
    # Write out candidates.csv and first classification csv file
    log.info("Writing out candidates.csv file in {}".format(second_round_path))
    t1_t2_all_meta_df.to_csv(second_round_path+'/candidates.csv', index=False)
    log.info("Writing out first classification file in {}".format(second_round_path))
    t1_t2_all_classified_df.to_csv(second_round_path+'/{}_classified1.csv'.format(opts.tag), index=False)
    
    # Tar up the directory
    tar_name = "{}_{}_round1_T1_T2.tar.gz".format(os.path.basename(opts.main_dir), opts.tag)
    make_tarfile(opts.output_dir, second_round_path, tar_name)    

    
    
if __name__ == "__main__":
    # Select options
    parser = optparse.OptionParser()
    parser.add_option('--name_tag', type=str, help = 'Name tag for searching classified files', dest='tag', default='autosave')
    parser.add_option('--main_dir', type=str, help='Root directory where all pointing subdirectories are stored', dest='main_dir')
    parser.add_option('--output_dir',type=str, help='Output directory where csvs/tar files will be written out', dest='output_dir')
    parser.add_option('--separate_csvs',type=int, help='Flag for writing out separate csvs for T1/T2 and pulsars', dest='sep_csv',default=0)
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
        opts.output_dir = opts.main_dir

    t1_t2_df = write_t1_t2_beams(opts)
    if opts.sep_csv:
        # Write out T1,T2 cand beam path+names
        log.info("Writing out beams for {} T1/T2 cands in total found in {}".format(t1_t2_df.shape[0], opts.main_dir))
        t1_t2_df.to_csv('{}/{}_{}_keep_T1_T2_cands.csv'.format(opts.output_dir, os.path.basename(opts.main_dir), opts.tag), index = False)
    
    
    kp_df = write_known_pulsar_beams(opts) 
    if opts.sep_csv:
        # Write out known pulsar beam path+names
        log.info("Writing out beams for {} known pulsar cands in total found in {}".format(t1_t2_df.shape[0], opts.main_dir))
        kp_df.to_csv('{}/{}_{}_keep_pulsar_cands.csv'.format(opts.output_dir, os.path.basename(opts.main_dir), opts.tag), index = False)

    # Write out a csv file for beams to be retained
    all_beams_df = t1_t2_df.append(kp_df, ignore_index=True)
    log.info("Writing out beams for all known pulsar and T1/T2 cands (Total: {}) found in {}".format(all_beams_df.shape[0], opts.main_dir))
    all_beams_df.drop_duplicates().to_csv('{}/{}_{}_keep_beams.csv'.format(opts.output_dir, os.path.basename(opts.main_dir), opts.tag), index = False)
  
    # Write out a second revision file
    log.info("Writing out a tarred file containing T1/T2 candidates")
    write_out_second_revision_tar_file(opts) 
 
