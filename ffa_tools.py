import numpy as np
#import cython
import matplotlib.pyplot as plt
import matplotlib.colors
import FFA_cy as FFA
import math
import time
from scipy import *
import scipy.signal
import sys
import re


def get_timeseries(beam):
	""" beam : str file name of the .dat file
	 Return the time series and the name of the beam (no .dat)	
	"""
	name = beam
	ts = list(np.fromfile(name,dtype='float32'))
	if name.endswith('.dat'):
		name = re.sub('\.dat$',"",name)
	return ts, name 
	
	
def get_info_beam(beam):
	""" beam : str file name of the .inf file
	 Return the lenght of the observation (T in sec), the sampling interval (dt in sec) and the DM	
	"""
	inffile = open(beam,'r')
	for line in inffile:
        	if line.startswith(" Number of bins in the time series"):
           	 	N = int(line.split()[-1])
        	if line.startswith(" Width of each time series bin (sec)"):
            		dt = float(line.split()[-1])
		if line.startswith(" Dispersion measure (cm-3 pc)"):
			DM = float(line.split()[-1])
	T = N * dt
	return T,dt,DM
	
	
def write_inf(name,file_to_write):
	inf_file = open(name + '.inf','r')
	for line in inf_file:
    		file_to_write.write(line)
	inf_file.close()
	return file_to_write
		
		
def select_factor(ts, m, mm):
    """
	select_factor(ts, m, mm):
	ts :  array (time series usually)
	m  : minimum factor
	mm : maximum factor
	It will delete the last element of ts until it has a factor within [m,mm]
	There is a maximum number of element that you delete in order to get a factor.
	This function returns the array ts (possibly shorter) and the first factor that is within the range [minimum,maximum]
    """
    ts = np.array(ts)
    a=np.array(factors(len(ts)))
    x = a[(a >=m) & (a <=mm)]
    counts =0
    while len(x)==0 and counts <20 :
	ts = np.delete(ts,-1)
	a=np.array(factors(len(ts)))
	x = a[(a >=m) & (a <=mm)]
	counts+=1
	if counts >= 20 :
		print "Having a hard time finding a downsampling factor to match the desired time resolution. "
		print "Try changing : 1) the desired time resolution 	2) by how much you are welling to vary from that time resolution   3) How much sample you are willing to delete (default = max 50) "
		sys.exit()

    return ts,x[0] 


def factors(x):
    """
	factors(x):
	x : int
	This function takes a number and returns the factors in an array, increasing order
    """
    facts = []
    for i in range(1, x + 1):
        if x % i == 0:facts.append(i)
    return facts


def make_factor_10(x):
	"""
	   make_factor_10(x)"
	   x : array
	   Deletes the last element of x until x is a factor of 10
	   return array factor of 10. 
	
	"""
	while len(x)%10!=0:
		x=np.delete(x,-1)
	return x

def normalize(lst):
    """
	normalize(lst):
	lst: Input must be 1-D
	returns the normalized array 
    """
    s = sum(lst)
    return array(map(lambda x: float(x)/s, lst))


def downsample(vector, factor):
    """
    downsample(vector, factor):
        Downsample (i.e. co-add consecutive numbers) a short section
            of a vector by an integer factor.
    """
    if (len(vector) % factor):
        print "Lenght of 'vector' is not divisible by 'factor'=%d!" % factor
        sys.exit()
    newvector = np.reshape(vector, (len(vector)/factor, factor))
    return np.add.reduce(newvector, 1)


def set_dws(data,downsamp,minimum,maximum):
    """ it will delete last element of the time series so that it can be downsampled by the 		chosen amount, unless it requires deleting more than 50 elements.
    """
    i=0
    data=list(data)
    d = 1
    while d == 1 and i<=60:
    	 data.pop()
   	 n=len(data)
   	 downsamp = select_factor(factors(n),minimum,maximum)
	 i+=1
    if i>60:
	 print'downsamp=1'
	 sys.exit() 												
    return array(data),downsamp	

def forced_dws_2phase(data):
		"""
		Returns TWO arrays : TWO different phases (i.e, 1+2,3+4,.. and 2+3,4+5,..
		"""	
		if len(data)%2 == 1:
			data = np.delete(data,-1)
		data1=downsample(data,2)
		f=data[0]
		l=data[-1]
		f_l=(f+l)/2
		data2=np.delete(data,0)
		data2=np.delete(data2,-1)
		data2=downsample(data2,2)
		data2=np.append(data2,f_l)
		return data1,data2


def forced_dws_3phase(data):
		"""
		Returns 3 arrays : 3 different phases (i.e, 1+2+3,4+5+6,.. and 2+3+4,5+6+7,.. 			and 3+4+5,6+7+8,...
		"""	
		if len(data)%3 == 1:
			data = np.delete(data,-1)
		if len(data)%3 == 2:
			data = np.delete(data,-1)
			data = np.delete(data,-1)
		data1=downsample(data,3)
		f=data[0]
		ff=data[1]
		l=data[-1]
		ll=data[-2]
		f_l_ll=(f+l+ll)/3
		f_ff_l=(f+ff+l)/3

		data2=np.delete(data,0)
		data2=np.delete(data2,-1)
		data2=np.delete(data2,-1)
		data2=downsample(data2,3)
		data2=np.append(data2,f_l_ll)

		data3=np.delete(data,0)
		data3=np.delete(data3,0)
		data3=np.delete(data3,-1)
		data3=downsample(data3,3)
		data3=np.append(data3,f_ff_l)


		return np.array(normalize(data1)),np.array(normalize(data2)),np.array(normalize(data2))		


def forced_dws(data,factor):
     while(len(data)%factor!=0):
         data=np.delete(data,-1)
     data=downsample(data,factor)
     return np.array(data)


def look_for_nan(SN):
	if SN.any()=='nan':
		print 'Signal to noise = nan'
		sys.exit()


#-----------------------	Signal-to-noise function	-----------------------
def simple_SNR(folds, sigma_total, w):
    """Return a very simple signal-to-noise for a profile.
       Works for narrow duty-cycle since  the S/N is max_value/std
	For each M profiles, returns a value of SNR (i.e, output is a list of lenght M)
    """
    M, P0 = folds.shape
    prof_std = 1.0/(np.sqrt(M)*sigma_total*np.sqrt(w))
    snr = (folds.max(axis=1)-np.median(folds, axis=1))*prof_std
    look_for_nan(snr)
    return snr



