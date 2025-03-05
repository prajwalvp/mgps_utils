import os
import sys
import glob
import shutil
import zipfile
import tarfile
import logging
import optparse
import pandas as pd
from pathlib import Path


log = logging.getLogger('candy_sweetener')
FORMAT = "[%(levelname)s - %(asctime)s - %(filename)s:%(lineno)s] %(message)s"
logging.basicConfig(format=FORMAT)
log.setLevel('INFO')


def unzip_file(opts):
    """ 
    Unzip the main downloaded zip folder named after viewer and containing Phase 2 viewing results 
    """
    with zipfile.ZipFile(opts.zip_path, 'r') as zip_ref:
        zip_ref.extractall(opts.output_path)


def untar_file(opts):
    """       
    Untar the tarred folder containing all T2 candidates of interest
    """
    tar_path = glob.glob(
        "{}/{}/*tar.gz".format(opts.output_path, opts.viewer))[0]
    my_tar = tarfile.open(tar_path)
    # specify which folder to extract to
    my_tar.extractall("{}/{}".format(opts.output_path, opts.viewer))
    my_tar.close()


if __name__ == "__main__":
    # Select options
    parser = optparse.OptionParser()
    parser.add_option('--input_zip', type=str,
                      help='Full path to zip file', dest='zip_path')
    parser.add_option('--epoch', type=str, help='Epoch of interest',
                      dest='epoch', default='20210701')
    parser.add_option('--output_dir', type=str,
                      help='Root directory where all pointing subdirectories are stored', dest='output_path')
    parser.add_option('--viewer_name', type=str,
                      help='Input name of viewer as specified in candidates repository', dest='viewer')
    parser.add_option('--viewer_threshold', type=int,
                      help='NUmber of approvals from viewers to select Tier2 candidates', dest='viewer_threshold', default=3)
    opts, args = parser.parse_args()

    if opts.zip_path is None:
        log.info("Zip path not specified")
        sys.exit(0)

    if opts.output_path is None:
        log.info("Output path not specified")
        log.info("Files will be written out in path of zip file")
        opts.output_path = os.path.dirname(opts.zip_path)

    # Unzip the folder
    unzip_file(opts)
    log.info("Folder unzipped")

    # Untar the T1/T2 candidates round1 tarred file
    untar_file(opts)
    log.info("Untarred T1/T2 candidates of interest")

    # Read in all the round2 csv files + round 1 csv file
    round2_csvs = glob.glob(
        '{}/{}/*round2*.csv'.format(opts.output_path, opts.viewer))
    round1_csv = glob.glob(
        "{}/{}/**/*classified1.csv".format(opts.output_path, opts.viewer))

    # Append all classifications to first classfication file
    classified_data = pd.read_csv(round1_csv[0])
    for i, csv in enumerate(round2_csvs):
        df = pd.read_csv(csv)
        classified_data['classification_round2_{}'.format(
            i+1)] = df['classification']

    # Apply conditions to select final viewing candidates
    final_cand_dir = "{}/{}_final_phase2_candidates".format(
        opts.output_path, opts.epoch)
    Path(final_cand_dir).mkdir(parents=True, exist_ok=True)
    log.info("Applying conditions: Looking for at least {} T2 ratings or 1 T1 rating per candidate".format(
        opts.viewer_threshold))
    for index, row in classified_data.iterrows():
        if (row.str.contains('T2_CAND|TIER2').sum()) >= opts.viewer_threshold or row.str.contains('T1_CAND|TIER1').sum() >= 1:
            png_file = classified_data['png'][index]
            log.info("{} meets threshold".format(os.path.basename(png_file)))
            png_path = glob.glob(
                '{}/{}/**/{}'.format(opts.output_path, opts.viewer, png_file))[0]
            shutil.copyfile(png_path, final_cand_dir +
                            '/'+os.path.basename(png_path))

    if len(os.listdir(final_cand_dir)) == 0:
        log.info("No candidates meet required threshold")
    else:
        log.info("Final candidates copied to {}".format(final_cand_dir))
