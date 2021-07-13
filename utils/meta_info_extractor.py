import re



# Supply a meta file from a MeerKAT pointing and this script gives the following information:
# 1. PNG file showing the beam tiling with known pulsars within the incoherent beam
# 2. PNG file showing the a zoomed version of the known pulsar within a tiling and surrounding beams
# 3. Check for associated Fermi sources in DR2 release and report them in a csv file
# 4. Gives csv files of all pulsars in field, separation from boresight, best beam to find the pulsar
# 5. Query HESS, Gaia via webparser (?) 
