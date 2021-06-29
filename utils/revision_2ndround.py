#Utility to create a small tarball with all the tier2 candidates
import glob

save_files=glob.glob("*/*miquel_save.csv")
candidates_collection=glob.glob("*/candidates.csv")
write_file=open("compress_2ndround.sh","w")
write_file.write("#!/bin/bash\n")
write_file.write("rm -r 2ndround\n")
write_file.write("mkdir 2ndround\n")

write_file.write("\n")
j=0
for name in save_files:
	print(name,candidates_collection[j])
	save_file=open(name,"r")
	candidates=open(candidates_collection[j],"r")
	i=0
	j=j+1
	for line in save_file:
		if i>0:
			line_vector=line.split(",")
			if line_vector[2]=="T1_CAND\n" or line_vector[2]=="T2_CAND\n":
				for candidate in candidates:
					candidate_vector=candidate.split(",")
					if line_vector[1]==candidate_vector[29]:
						pointing=candidate_vector[3]

						write_file.write("mkdir 2ndround/"+pointing+"\n")
						write_file.write("mkdir 2ndround/"+pointing+"/metafiles\n")
						write_file.write("mkdir 2ndround/"+pointing+"/plots\n")
						write_file.write("cp "+pointing+"/metafiles/*.meta 2ndround/"+pointing+"/metafiles/\n")
						write_file.write("cp "+pointing+"/"+line_vector[1]+" 2ndround/"+pointing+"/"+line_vector[1]+"\n")

						write_file.write("if test -f 2ndround/"+pointing+"/candidates.csv; then\n")
						write_file.write("	echo '"+candidate.split("\n")[0]+"' >> 2ndround/"+pointing+"/candidates.csv\n")
						write_file.write("else\n")
						write_file.write("	echo 'pointing_id,beam_id,beam_name,source_name,ra,dec,gl,gb,mjd_start,utc_start,f0_user,f0_opt,f0_opt_err,f1_user,f1_opt,f1_opt_err,acc_user,acc_opt,acc_opt_err,dm_user,dm_opt,dm_opt_err,sn_fft,sn_fold,pepoch,maxdm_ymw16,dist_ymw16,pics_trapum_ter5,pics_palfa,png_path,metafile_path,filterbank_path,candidate_tarball_path' > 2ndround/"+pointing+"/candidates.csv\n")
						write_file.write("	echo '"+candidate.split("\n")[0]+"' >> 2ndround/"+pointing+"/candidates.csv\n")
						write_file.write("fi\n")

						write_file.write("if test -f 2ndround/"+pointing+"/"+pointing+"_save.csv; then\n")
						write_file.write("	echo '"+line.split("\n")[0]+"' >> 2ndround/"+pointing+"/"+pointing+"_save.csv\n")
						write_file.write("else\n")
						write_file.write("	echo 'utc,png,classification' > 2ndround/"+pointing+"/"+pointing+"_save.csv\n")
						write_file.write("	echo '"+line.split("\n")[0]+"' >> 2ndround/"+pointing+"/"+pointing+"_save.csv\n")
						write_file.write("fi\n")
						write_file.write("\n")
				candidates.seek(0,0)
		i=i+1
write_file.write("tar -czvf 2ndround.tar.gz 2ndround\n")
save_file.close()
write_file.close()