import re
import os
import sys
import glob
import tarfile
import optparse

"""
This script untars specific tarballs only based on range of pointing IDs specified by user.
Assumes that the candidate cutter outputs small tarballs with the following naming convention: Survey-epoch_pointingID_pointingname.tar.gz
"""


def sorter(item):
    """Get an item from the list (one-by-one) and return a pattern required"""
    return os.path.basename(item).split('_')[1] 


def untar_candidates(opts):
    """
    Untar multiple sets of candidates based on pointing id ranges 
    """   
    all_tar_files = sorted(glob.glob('{}/*.tar.gz'.format(opts.path)),key=sorter)    
    start,end = [index for index,value in enumerate(all_tar_files) if os.path.basename(value).split('_')[1] == opts.start_id or os.path.basename(value).split('_')[1] == opts.end_id]
    reduced_tar_set = all_tar_files[start:end+1]
    print("Found {} tarballs between specified pointing ID ranges".format(len(reduced_tar_set)))
    for tar_file in reduced_tar_set:
        my_tar = tarfile.open(tar_file) 
        pointing_id = os.path.basename(tar_file).split('_')[1]
        my_tar.extractall('{}/{}'.format(opts.path, re.search("{}_(.*?).tar.gz".format(pointing_id), tar_file).group(1)))
        print("Untarred {} to {}/{}".format(tar_file, opts.path, re.search("{}_(.*?).tar.gz".format(pointing_id), tar_file).group(1)))
        my_tar.close()
    


if __name__=="__main__":
    # Select options
    parser = optparse.OptionParser()
    
    parser.add_option('--epoch_directory',type=str, help='Path to epoch of interest',dest='path')
    parser.add_option('--start_id',type=str, help='Starting pointing ID of range',dest='start_id')
    parser.add_option('--end_id',type=str, help='Ending pointing ID of range',dest='end_id')
    opts, args = parser.parse_args()

    untar_candidates(opts)
