#Utility to make a list of beams to save from good candidates, known pulsars or rfi from save and candidates files.
#Uncomment/comment lines 6,7,8 and 21,22,23 to decide which of them. 
import sys
import glob
import optparse
import pandas as pd


def write_t1_t2_beams(opts):
    """
    Write out T1, T2 beams 
    """ 

    columns = ['filterbank_path','beam_name']

    classified_files =  glob.glob("{}/**/*{}*.csv".format(opts.main_dir, opts.tag))    
    all_candidate_csvs = glob.glob("{}/**/candidates.csv".format(opts.main_dir))

    #print (classified_files, all_candidate_csvs)
        
    for i,classified_file in enumerate(classified_files):
        #print (classified_file, all_candidate_csvs[i])
        df = pd.read_csv(classified_file)
        t1_t2_df =  (df[(df['classification']=='T1_CAND') | (df['classification']=='T2_CAND')])
        
        if t1_t2_df.empty:
            continue         
        else:
            # Put condition to retrieve beam path and beam name for given dataframe


        # Store csv as keep_t1_t2_cands_tag_file

  
         
    
    #    print(name,candidates_collection[j])
    #    save_file=open(name,"r")
    #    candidates=open(candidates_collection[j],"r")
    #    i=0
    #    j=j+1
#       write_file.write(name+"\n")
    #    for line in save_file:
    #            if i>0:
    #                    line_vector=line.split(",")
    #                    if line_vector[2]=="T1_CAND\n" or line_vector[2]=="T2_CAND\n":
#                       if line_vector[2]=="KNOWN_PSR\n":


#write_file=open("keep_miquel.csv","w")
#write_file=open("known_miquel.csv","w")
#write_file=open("rfi_miquel.csv","w")
#write_file.write("# .fil path, beam name\n")
    




if __name__ == "__main__":
    # Select options
    parser = optparse.OptionParser()
    parser.add_option('--name_tag', type=str, help = 'Name tag for searching classified files', dest='tag', default='save')
    parser.add_option('--main_dir', type=str, help='Root directory where all pointing subdirectories are stored', dest='main_dir')
    parser.add_option('--output_dir',type=str, help='Output directory where csvs/tar files will be written out', dest='output_dir')
    opts, args = parser.parse_args()

    if opts.main_dir is None:
        print ("Pointing directory not specified!")
        sys.exit(0) 
    if opts.tag is None:
        print ("Tag to find classified files not specified!")
        sys.exit(0) 

    if opts.output_dir is None:
        print ("Files will be written in main specified directory")
        opts.output_dir = opts.main_dir

    # Write out T1,T2 cand beam path+names
    write_t1_t2_beams(opts)

    # Write out known pulsar beam path+names
    #write_known_pulsar_beams(opts)

    # Write out a tar file with T1 and T2 cands for multiple viewings
    #write_t1_t2_tar(opts)

 

   
#j=0
#for name in save_files:
#	print(name,candidates_collection[j])
#	save_file=open(name,"r")
#	candidates=open(candidates_collection[j],"r")
#	i=0
#	j=j+1
##	write_file.write(name+"\n")
#	for line in save_file:
#		if i>0:
#			line_vector=line.split(",")
#			if line_vector[2]=="T1_CAND\n" or line_vector[2]=="T2_CAND\n":
##			if line_vector[2]=="KNOWN_PSR\n":
##			if line_vector[2]=="RFI\n":
#				for candidate in candidates:
#					candidate_vector=candidate.split(",")
#					if line_vector[1]==candidate_vector[29]:
#						trapum_score=candidate_vector[27]
#						palfa_score=candidate_vector[28]
#						sn_fold=candidate_vector[23]
#						file=candidate_vector[31]
#						period=str(1/float(candidate_vector[11]))
#						path=file.rsplit("/",2)[0]
#						print(path+"/%,"+candidate_vector[2])
##						write_file.write(path+"/%,"+candidate_vector[2]+" ,"+sn_fold+"\n")
#						write_file.write(path+"/%,"+candidate_vector[2]+"\n")
#				candidates.seek(0,0)
#		i=i+1
#save_file.close()
#write_file.close()
