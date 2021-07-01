#Utility to make a list of beams to save from good candidates, known pulsars or rfi from save and candidates files.
#Uncomment/comment lines 6,7,8 and 21,22,23 to decide which of them. 
import glob
save_files=glob.glob("*/*miquel_save.csv")
candidates_collection=glob.glob("*/candidates.csv")
write_file=open("keep_miquel.csv","w")
#write_file=open("known_miquel.csv","w")
#write_file=open("rfi_miquel.csv","w")
write_file.write("# .fil path, beam name\n")
j=0
for name in save_files:
	print(name,candidates_collection[j])
	save_file=open(name,"r")
	candidates=open(candidates_collection[j],"r")
	i=0
	j=j+1
#	write_file.write(name+"\n")
	for line in save_file:
		if i>0:
			line_vector=line.split(",")
			if line_vector[2]=="T1_CAND\n" or line_vector[2]=="T2_CAND\n":
#			if line_vector[2]=="KNOWN_PSR\n":
#			if line_vector[2]=="RFI\n":
				for candidate in candidates:
					candidate_vector=candidate.split(",")
					if line_vector[1]==candidate_vector[29]:
						trapum_score=candidate_vector[27]
						palfa_score=candidate_vector[28]
						sn_fold=candidate_vector[23]
						file=candidate_vector[31]
						period=str(1/float(candidate_vector[11]))
						path=file.rsplit("/",2)[0]
						print(path+"/%,"+candidate_vector[2])
#						write_file.write(path+"/%,"+candidate_vector[2]+" ,"+sn_fold+"\n")
						write_file.write(path+"/%,"+candidate_vector[2]+"\n")
				candidates.seek(0,0)
		i=i+1
save_file.close()
write_file.close()
