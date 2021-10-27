Running metaxtractor

```
Usage: metaxtractor.py [options]


Options:
  -h, --help            show this help message and exit
  --meta_path=META      Path to meta file
  --check_unpublished=UNPUBLISHED_FLAG
                        Flag for checking unpublished sources
  --sheet_id=SHEET_ID   Unique spreadheet ID
  --sheet_name=SHEET_NAME
                        Name of sheet of unpublished pulsars
  --check_fermi=FERMI_FLAG
                        Check for possible Fermi associations
  --fits_file=FITS_FILE
                        Full path for Fermi fits file
  --known_pulsar=KP_CATALOGUE
                        Cross match with known pulsar list as specified.
                        Options: ATNF, PSS
  --user_coordinate=USER_COORDS
                         Plot user specified coordinates  e.g. 12:08 -59:36
  --user_beam_radius=BEAM_RADIUS
                        Plot user specified beam radius in degrees
  --output_path=OUTPUT_PATH
                        Path to store output files
```

Prerequisite files:

```
4FGL_DR2_Ppsr.fits (Provided in repository)
get_psrs_in_field.py (Web parser for Pulsar survey scraper; Provided in repository)
```

Packages needed:
```
os, sys, math, logging, ast, optparse, subprocess, numpy, matplotlib, pandas, astropy, psrqpy, bs4, urllib, requests_html   


Example command:
```
python metaxtractor.py --meta_path /Users/mgps/candidates/20210411/MSGPS_L_2063/metafiles/2021-04-11T18:05:24.meta --known_pulsar PSS --output_path /test
```


Extra instructions/notes:

1. Unpublished source matching currently works for HTRU only. Make sure the '--check_unpublished' is 1 and provide sheet ID of Google spreadsheet URL
2. Default user_beam_radius is the survey beam radius at L-Band for MGPS.
3. There are constants defined biased to MGPS-L. These can be changed for other receiver runs.
