#Script that reads save files and candidate files from CandyJar and makes a persto script to refold T1 and T2 candidates.
#It performs the peasoup-prepfold epoch correction on frequencies.
import glob
#Collect the save files.
save_files=sorted(glob.glob("*/*miquel_save.csv"))
#Open the file where the folding script will be written.
write_file=open("refold.sh","w")
#Write a header line onto it.
write_file.write("#!/bin/bash\n")
#Collect the candidates files.
candidates_collection=sorted(glob.glob("*/candidates.csv"))
j=0
for name in save_files:
	print(name,candidates_collection[j])
	#Open the save file.
	save_file=open(name,"r")
	#Open the candidate file.
	candidates=open(candidates_collection[j],"r")
	i=0
	j=j+1
	for line in save_file:
		#Loop through you classification save file.
		if i>0:
			#Spilt the classification line.
			line_string=line.split(",")
			#Select tier 1 and tier 2 candidates.
			if line_string[2]=="T1_CAND\n" or line_string[2]=="T2_CAND\n":
				#Loop over the candidates in the candidates file.
				for candidate in candidates:
					#Split the candidates line.
					candidate_string=candidate.split(",")
					#Select the candidates with coinciding png id.
					if line_string[1]==candidate_string[29]:
						file_name=candidate_string[31]
						path=file_name.rsplit("/",1)[0]
#						path_0th=path.rsplit("/",1)[0]+"/cfbf00000"
						name=file_name.rsplit("/",1)[1]
						name=name.rsplit("_",1)[0]
						name_mask=name.rsplit("_",1)[0]
						#Write a comment line indicating beam and pointing.
						write_file.write("#---------------------------------------------------------#\n")
						write_file.write("#					                                        #\n")
						write_file.write("#Candidate "+name+"----------------#\n")
						#Write the copy commmand to your local folder.
						files="$(find /data"+path+" -name '"+name_mask+"*')"
						write_file.write("if test -f "+name_mask+"_rfifind.mask; then\n")
						write_file.write("	echo '"+name_mask+"_rfifind.mask exists'\n")
						write_file.write("else\n")
						write_file.write("	ln -s "+files+" .\n")
						rfi_command="	rfifind -ncpus 8 -time 5 -freqsig 6 -o "+name_mask+" *.fil"
						#Write the rfifind command.
						write_file.write(rfi_command+"\n")
						write_file.write("	rm *.fil\n")
						write_file.write("fi\n")
						#Copy the files again.
						files="$(find /data"+path+" -name '"+name+"*')"
						write_file.write("ln -s "+files+" .\n")
						# pulsarix gives the f0 at the middle epoch, but presto assumes that you give it at the starting epoch.
						#Compute the f0 value at the observation start from the period derivative.
						frequency_derivative=str(candidate_string[14])
						epoch_time=str((float(candidate_string[24])-float(candidate_string[8]))*24*3600)
						print(epoch_time)
						folding_frequency=str(float(candidate_string[11])-float(frequency_derivative)*float(epoch_time))
						fold_command="prepfold -ncpus 8 -zerodm -nosearch -topo -mask "+name_mask+"_rfifind.mask -f "+folding_frequency+" -fd "+frequency_derivative+" -dm "+candidate_string[20]+" *.fil"
						#Write the folding commnad.
						write_file.write(fold_command+"\n")
						#Compute the frequencie (and its derivative) of the second harmonic for the candidate.
#						lower_f0=str(float(folding_frequency)/2)
#						lower_f1=str(float(frequency_derivative)/2)
#						fold_command="prepfold -ncpus 8 -zerodm -mask "+name_mask+" -f "+lower_f0+" -fd "+lower_f1+" -dm "+candidate_string[20]+"-nodmsearch -topo *.fil"
						#Write the folding command for the second harmonic.
#						write_file.write(fold_command+"\n")
						#Remove the files from your local directory.
						write_file.write("rm *.fil\n")
				candidates.seek(0,0)
		i=i+1
	save_file.close()
write_file.close()
