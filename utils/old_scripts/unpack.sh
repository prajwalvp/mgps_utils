#!bin/bash
#Utility to unpack a bucnh of tarballs into files with MGPS_L_{}.
for tarball in $(ls *.tar.gz); do
	IFS="_" read -a pointing <<< $tarball
	IFS="." read -a directory <<< ${pointing[2]}_${pointing[3]}_${pointing[4]}
	mkdir ${directory[0]}
	tar -xvzf ${tarball} -C ${directory[0]}
done
