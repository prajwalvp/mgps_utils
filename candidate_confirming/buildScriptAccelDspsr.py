#Script to write a dspsr/pdmp refolding script of candidates from a peasoup search output file.
import sys
import glob
from ast import literal_eval

if len(sys.argv)<5 or sys.argv[1]=="h" or sys.argv[1]=="help":
	print("Not enough arguments, please indicate:")
	print("	1.- Path to observation filterbank files. It has to include the metafile.")
	print("	2.- The candidate sumary files.")
	print("	3.- The file where you want to write the refolding commands.")
	print(" 4.- Range of periods to check for sperated by a coma. Give it a range, and only candidates with 1 harmonic up and 4 harmonics down will be selected and refolded.")
	print("You should also indicate:")
	print("	5.- List of extra channels that you zant to zap befor going to pdmp, given as a list separated by spaces from 0 to 63 (Default: none)")
	print(' Example: >> python3 buildScriptAccelDspsr.py "/beegfs/DATA/TRAPUM/SCI-20200703-MK-01/20210515-0015/20210519_223207/*" "/beegfs/PROCESSING/TRAPUM/peasoup_ddplan_20210401/2021/5/20/*_1306_6043_confirmation/*.xml" refold.sh 0.005665,0.005675 12')
	sys.exit()
#The filterbank file.
observation_files=sorted(glob.glob(sys.argv[1]))
metafile=observation_files[0]
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
	zap=True
	zapping_list=sys.argv[5]
else:
	zap=False
#Loop throgh candidates files
i=0
for candidate_file in candidate_files:
	print(" ")
	print("Reading candidates from file "+candidate_file)
	print("From beam "+observation_files[i])
	candidates=open(candidate_file,"r")
	obs_file_path=observation_files[i+1]
	if i<=9:
		beam="cfbf0000"+str(i)
	elif i<=99:
		beam="cfbf000"+str(i)
	elif i==100:
		beam="ifbf00000"
	files="$(find /data"+obs_file_path+"/ -name '*.fil')"
	#Transform the mask into a format readable by psrchive.
	if i==0:
		rfi_command="python getout_rfifind.py cfbf00000_rfifind.mask cfbf00000"
		write_file.write(rfi_command+"\n")
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
				with open(metafile,'r') as f:
					data = f.read()
				all_info = literal_eval(data)
				beam_info=all_info['beams']
				ra = beam_info[beam].split(',')[2]     
				dec = beam_info[beam].split(',')[3]
				name="J"+ra.split(":")[0]+ra.split(":")[1]+dec.split(":")[0]+dec.split(":")[1]
				name=name.replace(" ","")
				predictor=open(name+"_"+beam+"_predictor.txt","w")
				predictor.write("SOURCE: "+name+"_"+beam+"\n")
				predictor.write("PERIOD: "+period+"\n")
				predictor.write("DM: "+dm+"\n")
				predictor.write("PDOT: "+p_dot+"\n")
				predictor.write("RA: "+ra+"\n")
				predictor.write("DEC: "+dec)
				predictor.close()
				write_file.write("#Candidate "+cand_id+" from beam "+beam+" with S/N "+snr+"\n")
				write_file.write("ln -s "+files+" .\n")
				fold_command="dspsr -P "+name+"_"+beam+"_predictor.txt -L 30 -b 64 -O "+name+"_"+beam+" -e .ar *.fil"                                                 
				add_command="psradd -o "+name+"_"+beam+".ar "+name+"_"+beam+"_*.ar"
				edit_command="psredit -m -c site='MeerKAT' "+name+"_"+beam+".ar"
				paz_command="paz -e .paz -k cfbf00000.badchan "+name+"_"+beam+".ar"
				pam_command="pam -m --setnchn 64 "+name+"_"+beam+".paz"
				pdmp_command="pdmp -mc 64 -ao "+acc+" -ar 2 -g "+name+"_"+beam+"_scrunched.png/PNG "+name+"_"+beam+".paz"
				write_file.write(fold_command+"\n")
				write_file.write("rm *.fil\n")
				write_file.write(add_command+"\n")
				write_file.write("rm *_0*\n")
				write_file.write(edit_command+"\n")
				write_file.write(paz_command+"\n")
				write_file.write(pam_command+"\n")
				if zap==True:
					paz_command_2="paz -m -z '"+zapping_list+"' "+name+"_"+beam+".paz"
					write_file.write(paz_command_2+"\n")
				write_file.write(pdmp_command+"\n")
				write_file.write("#---------------------------------------------#\n")
				break
			elif 2*float(period) > starting_period and 2*float(period) < ending_period:
				with open(metafile,'r') as f:
					data = f.read()
				all_info = literal_eval(data)
				beam_info=all_info['beams']
				ra = beam_info[beam].split(',')[2]     
				dec = beam_info[beam].split(',')[3]
				name="J"+ra.split(":")[0]+ra.split(":")[1]+dec.split(":")[0]+dec.split(":")[1]
				name=name.replace(" ","")
				predictor=open(name+"_"+beam+"_predictor.txt","w")
				predictor.write("SOURCE: "+name+"_"+beam+"\n")
				predictor.write("PERIOD: "+str(2*float(period))+"\n")
				predictor.write("DM: "+dm+"\n")
				predictor.write("PDOT: "+str(2*float(p_dot))+"\n")
				predictor.write("RA: "+ra+"\n")
				predictor.write("DEC: "+dec)
				predictor.close()
				write_file.write("#Candidate "+cand_id+" from beam "+beam+" with S/N "+snr+"\n")
				write_file.write("ln -s "+files+" .\n")
				fold_command="dspsr -P "+name+"_"+beam+"_predictor.txt -L 30 -b 64 -O "+name+"_"+beam+" -e .ar *.fil"                                               
				add_command="psradd -o "+name+"_"+beam+".ar "+name+"_"+beam+"_*.ar"
				edit_command="psredit -m -c site='MeerKAT' "+name+"_"+beam+".ar"
				paz_command="paz -e .paz -k cfbf00000.badchan "+name+"_"+beam+".ar"
				pam_command="pam -m --setnchn 64 "+name+"_"+beam+".paz"
				pdmp_command="pdmp -mc 64 -ao "+acc+" -ar 2 -g "+name+"_"+beam+"_scrunched.png/PNG "+name+"_"+beam+".paz"
				write_file.write(fold_command+"\n")
				write_file.write("rm *.fil\n")
				write_file.write(add_command+"\n")
				write_file.write("rm *_0*\n")
				write_file.write(edit_command+"\n")
				write_file.write(paz_command+"\n")
				write_file.write(pam_command+"\n")
				if zap==True:
					paz_command_2="paz -m -z '"+zapping_list+"' "+name+"_"+beam+".paz"
					write_file.write(paz_command_2+"\n")
				write_file.write(pdmp_command+"\n")
				write_file.write("#---------------------------------------------#\n")
				break
			elif float(period)/2 > starting_period and float(period)/2 < ending_period:
				with open(metafile,'r') as f:
					data = f.read()
				all_info = literal_eval(data)
				beam_info=all_info['beams']
				ra = beam_info[beam].split(',')[2]     
				dec = beam_info[beam].split(',')[3]
				name="J"+ra.split(":")[0]+ra.split(":")[1]+dec.split(":")[0]+dec.split(":")[1]
				name=name.replace(" ","")
				predictor=open(name+"_"+beam+"_predictor.txt","w")
				predictor.write("SOURCE: "+name+"_"+beam+"\n")
				predictor.write("PERIOD: "+str(float(period)/2)+"\n")
				predictor.write("DM: "+dm+"\n")
				predictor.write("PDOT: "+str(float(p_dot)/2)+"\n")
				predictor.write("RA: "+ra+"\n")
				predictor.write("DEC: "+dec)
				predictor.close()
				write_file.write("#Candidate "+cand_id+" from beam "+beam+" with S/N "+snr+"\n")
				write_file.write("ln -s "+files+" .\n")
				fold_command="dspsr -P "+name+"_"+beam+"_predictor.txt -L 30 -b 64 -O "+name+"_"+beam+" -e .ar *.fil"                                                  
				add_command="psradd -o "+name+"_"+beam+".ar "+name+"_"+beam+"_*.ar"
				edit_command="psredit -m -c site='MeerKAT' "+name+"_"+beam+".ar"
				paz_command="paz -e .paz -k cfbf00000.badchan "+name+"_"+beam+".ar"
				pam_command="pam -m --setnchn 64 "+name+"_"+beam+".paz"
				pdmp_command="pdmp -mc 64 -ao "+acc+" -ar 2 -g "+name+"_"+beam+"_scrunched.png/PNG "+name+"_"+beam+".paz"
				write_file.write(fold_command+"\n")
				write_file.write("rm *.fil\n")
				write_file.write(add_command+"\n")
				write_file.write("rm *_0*\n")
				write_file.write(edit_command+"\n")
				write_file.write(paz_command+"\n")
				write_file.write(pam_command+"\n")
				if zap==True:
					paz_command_2="paz -m -z '"+zapping_list+"' "+name+"_"+beam+".paz"
					write_file.write(paz_command_2+"\n")
				write_file.write(pdmp_command+"\n")
				write_file.write("#---------------------------------------------#\n")
				break
	candidates.close()
	i=i+1
	print(" ")
write_file.close()