import ffa_tools as f
#import ffa_code as c
import ffa_code_withplots as c
import FFA_cy as FFA
import numpy as np
import math
import time
import scipy.signal
import inspect
import sys
import os
import ConfigParser
import ast
import re
total_time=time.time()


class ffa_cands(object):
    """
    FFA candidates has 3 parameters: the period (sec), the SNR and the sampling interval, dt (sec). 
    """
    def __init__(self):
	self.periods = np.array([])
	self.SNRs =  np.array([])
	self.dts = np.array([])

    def add_cand(self,p,SNR,dt):
	self.periods = np.insert(self.periods,0, p)
	self.SNRs = np.insert(self.SNRs,0, SNR)
	self.dts = np.insert(self.dts,0,dt)

    def print_content(self):
	print "  Period : ", self.periods,  ", SNR: ", self.SNRs, ",dt: ", self.dts
	print "	len(P): ",len(self.periods), ",   len(SNR): ",len(self.SNRs)," ,   len(dts): ",len(self.dts)


# ------------- 	Initial		------------------
cands=ffa_cands()
cfg = ConfigParser.ConfigParser()
cfg.read('config_ffa.cfg')
#cfg.read(str(sys.argv[2]))


# ------------- 	Read configuration file		------------------

p_ranges = eval(cfg.get('FFA_settings','p_ranges'))
dt_list = eval(cfg.get('FFA_settings','dt_list'))
SN_tresh = float(cfg.get('FFA_settings','SN_tresh'))


plot_period = ast.literal_eval(cfg.get('extra_settings','plot_period'))
plot_Ms = ast.literal_eval(cfg.get('extra_settings','plot_Ms'))
write_cands = ast.literal_eval(cfg.get('extra_settings','write_cands'))
do_large_dc = ast.literal_eval(cfg.get('extra_settings','do_large_dc'))
min_dc = ast.literal_eval(cfg.get('extra_settings','min_dc'))



# ------------- 	Select beam	------------------
beam = sys.argv[1]
ts,name = f.get_timeseries(beam)
T,dt,DM = f.get_info_beam(name+'.inf')

# ------------- 	Downsampling	------------------

if min_dc > 0.5:
	dt_list = 2*dt_list  
dwn_ideal = int(dt_list[0]/dt)
minimum_dwn = dwn_ideal-int(dwn_ideal*0.05)
maximum_dwn = dwn_ideal+int(dwn_ideal*0.15)

ts,dwn = f.select_factor(ts,minimum_dwn,maximum_dwn)
print 'Downsampling is being performed'
ts = f.downsample(ts, dwn)

# ------------- 	Detrending	------------------
print "Doing the detrending"
window_size = 30*int(len(ts)/T)		#window size over which statistics are computed 
break_points = np.arange(0,len(ts),window_size) 
ts = scipy.signal.detrend(ts,bp=break_points)

ts2=ts
sigma_total = np.std(ts)
# ------------- 		FFA		------------------
dt= T/len(ts)
count_lim = 2		# used in stage 2 and 3; how many consecutive downsamplings
print "Initial time resolution: ",'%.4f'%(dt),'/n'
print "For period range of :", p_ranges[0], ", 	  the sampling interval is : ", '%.5f'%(dt)
c.ffa_code_stage1(ts, dt,T, sigma_total,p_ranges[0][0],p_ranges[0][1],SN_tresh, count_lim,name, cands)
c.ffa_code_stage2(ts, dt, T, sigma_total, p_ranges[0][0],p_ranges[0][1],SN_tresh, count_lim,name, cands)
c.ffa_code_stage3(ts, dt, T, sigma_total, p_ranges[0][0],p_ranges[0][1],SN_tresh, count_lim,name, cands)

#Going through subsets of periods 
for num_p_ranges in range(len(p_ranges)-1):
	rang = num_p_ranges+1 
	dwn_ideal = int(dt_list[rang]/dt)
	if dwn_ideal ==1 :
		dwn_ideal =2 
	minimum_dwn = dwn_ideal-int(dwn_ideal*0.05)
	maximum_dwn = dwn_ideal+int(dwn_ideal*0.15)

	ts,dwn = f.select_factor(ts,minimum_dwn,maximum_dwn)
	ts = f.downsample(ts, dwn)
	dt = T/len(ts)
	print "For period range of :", p_ranges[rang], ",    the sampling interval is : ", '%.5f'%(dt)
	c.ffa_code_stage1(ts, dt, T, sigma_total, p_ranges[rang][0],p_ranges[rang][1], SN_tresh, count_lim,name, cands)
	c.ffa_code_stage2(ts, dt, T, sigma_total, p_ranges[rang][0],p_ranges[rang][1], SN_tresh, count_lim,name, cands)
	c.ffa_code_stage3(ts, dt, T, sigma_total, p_ranges[rang][0],p_ranges[rang][1], SN_tresh, count_lim,name, cands)

# ------------- 		end of FFA		------------------



if do_large_dc:		

	# __________________________ Start of DO_LARGE_DC ____________________________

	print  '\n','\n'
	print  '	------------------     					-----------------------'
	print  '	------------------     Doing larger duty-cycles 	-----------------------'
	print  '	------------------     					-----------------------'


	# ------------- 	Downsampling - extra	------------------
	dt_list2=10*dt_list
	dwn_ideal = 10
	print "Doing the downsampling"
	ts2 = f.forced_dws(ts2, dwn_ideal)
	time_d = time.time() - time_downsamp
	 
	
	# ------------- 	Detrending - extra------------------
	print "Doing the detrending"
	window_size = 30*int(len(ts2)/T)
	break_points = np.arange(0,len(ts2),window_size) 
	ts2 = scipy.signal.detrend(ts2,bp=break_points)
	
	
	# ------------- 		FFA - extra	------------------
	count_lim = 1		# stops at 4*dt in stage 2 and at 9*dt in stage 3, otherwise duty-cycles~100%
	print "For period range of :", p_ranges[0], ",    the sampling interval is : ", '%.5f'%(dt)
	c.ffa_code_stage1(ts2, dt, T, p_ranges[0][0],p_ranges[0][1],SN_tresh, count_lim,name, cands)
	c.ffa_code_stage2(ts2, dt, T, p_ranges[0][0],p_ranges[0][1],SN_tresh, count_lim,name, cands)
	c.ffa_code_stage3(ts2, dt, T, p_ranges[0][0],p_ranges[0][1],SN_tresh, count_lim,name, cands)

	#Going through subsets of periods 
	for num_p_ranges in range(len(p_ranges)-1):
		rang = num_p_ranges+1 
		dwn_ideal = int(np.ceil(dt_list[rang]/dt))
		minimum_dwn = dwn_ideal-int(dwn_ideal*0.15)
		maximum_dwn = dwn_ideal+int(dwn_ideal*0.35)
		ts2,dwn = f.select_factor(ts2,minimum_dwn,maximum_dwn)
		if dwn !=1 or dwn !=0:
			ts2 = f.downsample(ts2, dwn)
		dt = T/len(ts2)
		print "For period range of :", p_ranges[rang], ",    the sampling interval is : ", '%.5f'%(dt)
		c.ffa_code_stage1(ts2, dt, T, p_ranges[rang][0],p_ranges[rang][1],SN_tresh, count_lim,name,cands)
		c.ffa_code_stage2(ts2, dt, T, p_ranges[rang][0],p_ranges[rang][1],SN_tresh, count_lim,name, cands)
		c.ffa_code_stage3(ts2, dt, T, p_ranges[rang][0],p_ranges[rang][1],SN_tresh, count_lim,name, cands)


	# __________________________ End of DO_LARGE_DC ____________________________




# ------------- 		Write candidates to file		------------------
if write_cands:
	k=int(1)
	fo = open(name+'_cands.ffa','w')
	for i in range(len(cands.periods)):
		fo.write(str(k)+'\t'+'\t'+str(cands.periods[i])+'\t'+'\t'+str(cands.SNRs[i])	+'\t'+'\t'+str(cands.dts[i])+'\n')
		k+=1
	print "Wrote ", str(k), "candidates in : "+name+'_cands.ffa'
	f.write_inf(name,fo)
	fo.close()


# ------------- 		End		------------------
print "FINISHED with ",name
time_tot =  (time.time() - total_time)
print ( " --- %.7s seconds is the FFA total time ---" % time_tot),'\n'
	
	











