#Script to write a presto refolding script of candidates from a peasoup search output file.
import sys
import glob
if len(sys.argv)<5 or sys.argv[1]=="h" or sys.argv[1]=="help":
	print("Not enough arguments, please indicate:")
	print("	1.- Path to observation filterbank files.")
	print("	2.- The candidate sumary files.")
	print("	3.- The file where you want to write the refolding commands.")
	print(" 4.- Range of periods to check for sperated by a coma. Give it a range, and only candidates with 1 harmonic up and 1 harmonic down will be selected and refolded.")
	print("You should also indicate:")
	print("	5.- The amount of desired cpus (Integer, default: 1)")
	print(' Example: >> python3 buildScriptAccelPrep.py "/beegfs/DATA/TRAPUM/SCI-20200703-MK-01/20210515-0015/20210519_223207/*fbf*" "/beegfs/PROCESSING/TRAPUM/peasoup_ddplan_20210401/2021/5/20/*_1306_6043_confirmation/*.xml" refold.sh 0.005665,0.005675 8')
	sys.exit()
#The filterbank file.
observation_files=sorted(glob.glob(sys.argv[1]))
#Open the file with candidates.
candidate_files=sorted(glob.glob(sys.argv[2]))
#Open the file where the folding script will be written.
write_file=open(sys.argv[3],"w")
#Write a header line onto it.
write_file.write("#!/bin/bash\n")
#Select the period ranges.
starting_period=float(sys.argv[4].split(",")[0])
ending_period=float(sys.argv[4].split(",")[1])
if len(sys.argv)>5:
	ncpus=str(int(sys.argv[5]))
else:
	ncpus="1"
#Loop throgh candidates files
i=0
for candidate_file in candidate_files:
	print(" ")
	print("Reading candidates from file "+candidate_file)
	print("From beam "+observation_files[i])
	candidates=open(candidate_file,"r")
	obs_file_path=observation_files[i]
	if i<=9:
		beam="cfbf0000"+str(i)
	elif i<=99:
		beam="cfbf000"+str(i)
	elif i==100:
		beam="ifbf00000"
	files="$(find /data"+obs_file_path+"/ -name '*.fil')"
	#Write the rfifind command
	if i==0:
		write_file.write("ln -s "+files+" .\n")
		rfi_command="rfifind -ncpus "+ncpus+" -time 5 -timesig 6 -freqsig 8 -o cfbf00000 *.fil"
		write_file.write(rfi_command+"\n")
		write_file.write("rm *.fil\n")
	#Loop through all the lines until you find a cadidate.
	for line in candidates:
		#Split the line by the "'" character.
		line=line.split(">")
		#Find the sampling time.
		if line[0]=="    <tsamp":
			tsamp=float((line[1]).split("<")[0])
		#Find the number of time samples
		if line[0]=="    <nsamples":
			nsamples=float((line[1]).split("<")[0])
			#Compute observation length.
			length=tsamp*nsamples
		#If the line corresponds to a candiate, record its porperties.
		line=(line[0]).split("'")
		if line[0]=="    <candidate id=":
			cand_id=line[1]
			period=(candidates.readline().split(">")[1]).split("<")[0]
			opt=candidates.readline()
			dm=(candidates.readline().split(">")[1]).split("<")[0]
			acc=(candidates.readline().split(">")[1]).split("<")[0]
			caca=candidates.readline()
			snr=(candidates.readline().split(">")[1]).split("<")[0]
			# accel_peasoup gives the p at the middle epoch, but presto assumes that you give it at the starting epoch
			#Compute the P0 value at the observation start from the acchel and length.
			period=str(float(period)*(1-(length*float(acc))/(2*299792458)))
			#Do the same with the derivative.
			p_dot=str(float(period)*float(acc)/299792458)
			#Write the folding commands for the correspondign harmonics.
			if float(period) > starting_period and float(period) < ending_period:
				write_file.write("#Candidate "+cand_id+" from beam "+beam+" with S/N "+snr+"\n")
				write_file.write("ln -s "+files+" .\n")
				fold_command="prepfold -ncpus "+ncpus+" -zerodm -p "+str(float(period))+" -pd "+str(float(p_dot))+" -dm "+dm+" -nosearch -topo -mask cfbf00000_rfifind.mask -o "+beam+"  *.fil"                                                 
				write_file.write(fold_command+"\n")
				write_file.write("rm *.fil\n")
				write_file.write("#---------------------------------------------#\n")
				break
			elif 2*float(period) > starting_period and 2*float(period) < ending_period:
				write_file.write("#Candidate "+cand_id+" from beam "+beam+" with S/N "+snr+"\n")
				write_file.write("ln -s "+files+" .\n")
				fold_command="prepfold -ncpus "+ncpus+" -zerodm -p "+str(2*float(period))+" -pd "+str(2*float(p_dot))+" -dm "+dm+" -nosearch -topo -mask cfbf00000_rfifind.mask -o "+beam+"  *.fil"                                      
				write_file.write(fold_command+"\n")
				write_file.write("rm *.fil\n")
				write_file.write("#---------------------------------------------#\n")
				break
			elif float(period)/2 > starting_period and float(period)/2 < ending_period:
				write_file.write("#Candidate "+cand_id+" from beam "+beam+" with S/N "+snr+"\n")
				write_file.write("ln -s "+files+" .\n")
				fold_command="prepfold -ncpus "+ncpus+" -zerodm -p "+str(float(period)/2)+" -pd "+str(float(p_dot)/2)+" -dm "+dm+" -nosearch -topo --mask cfbf00000_rfifind.mask -o "+beam+"  *.fil"                            
				write_file.write(fold_command+"\n")
				write_file.write("rm *.fil\n")
				write_file.write("#---------------------------------------------#\n")
				break
	candidates.close()
	i=i+1
	print(" ")
write_file.close()