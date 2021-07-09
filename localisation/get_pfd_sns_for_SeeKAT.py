from presto.prepfold import pfd
from presto.presto import chi2_sigma
from scipy.stats import chi2, norm
import os
import sys
import numpy as np
import glob
from astropy.io import ascii

data = glob.glob('*.pfd')
data = np.sort(data)

RA=[]
DEC=[]
SNR1=[]
SNR2=[]
SNR3=[]
SNR4=[]
RA_all=[]
DEC_all=[]
SNR_all=[]
NAME_all=[]
SNRok=[]
NAME=[]

for i in range(0,len(data)):
    folded_data = pfd(data[i])
    redchi2 = folded_data.calc_redchi2()
    DoF = folded_data.DOFcor
    ra = folded_data.rastr
    dec = folded_data.decstr
    nsub = np.shape(folded_data.profs)[1]
    sigma = norm.isf(chi2.sf(redchi2 * DoF,DoF),scale=1.0)
    RA.append(ra)
    DEC.append(dec)
    SNRok.append(sigma)
    NAME.append(data[i])
    print(data[i],'Total S/N:',sigma)

    folded_data = pfd(data[i])
    folded_data.kill_subbands(range(nsub//4,nsub))
    redchi2 = folded_data.calc_redchi2()
    ra = folded_data.rastr
    dec = folded_data.decstr
    DoF = folded_data.DOFcor
    sigma = norm.isf(chi2.sf(redchi2 * DoF,DoF),scale=1.0)
    SNR1.append(sigma)
    RA_all.append(ra)
    DEC_all.append(dec)
    SNR_all.append(sigma)
    NAME_all.append(data[i] + '.0')
    print(data[i],'Subband 1 S/N:',sigma)

    folded_data = pfd(data[i])
    folded_data.kill_subbands(range(0,nsub//4))
    folded_data.kill_subbands(range(nsub//2,nsub))
    redchi2 = folded_data.calc_redchi2()
    DoF = folded_data.DOFcor
    sigma = norm.isf(chi2.sf(redchi2 * DoF,DoF),scale=1.0)
    SNR2.append(sigma)
    RA_all.append(ra)
    DEC_all.append(dec)
    SNR_all.append(sigma)
    NAME_all.append(data[i] + '.1')
    print(data[i],'Subband 2 S/N:',sigma)

    folded_data = pfd(data[i])
    folded_data.kill_subbands(range(0,nsub//2))
    folded_data.kill_subbands(range(3*nsub//4,nsub))
    redchi2 = folded_data.calc_redchi2()
    DoF = folded_data.DOFcor
    sigma = norm.isf(chi2.sf(redchi2 * DoF,DoF),scale=1.0)
    SNR2.append(sigma)
    RA_all.append(ra)
    DEC_all.append(dec)
    SNR_all.append(sigma)
    NAME_all.append(data[i] + '.2')
    print(data[i],'Subband 3 S/N:',sigma)

    folded_data = pfd(data[i])
    folded_data.kill_subbands(range(0,3*nsub//4))
    redchi2 = folded_data.calc_redchi2()
    DoF = folded_data.DOFcor
    sigma = norm.isf(chi2.sf(redchi2 * DoF,DoF),scale=1.0)
    SNR2.append(sigma)
    RA_all.append(ra)
    DEC_all.append(dec)
    SNR_all.append(sigma)
    NAME_all.append(data[i] + '.3')
    print(data[i],'Subband 4 S/N:',sigma)

ascii.write([RA,DEC,NAME,SNRok], 'beam_total_SNs.txt',names=['#ra','dec','name','snr'], fast_writer=False)                                                    
ascii.write([NAME_all,SNR_all], 'beam_subband_SNs.txt',names=['#name','snr'],fast_writer=False)
