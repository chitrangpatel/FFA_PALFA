import ffa_tools as f
import FFA_cy as FFA
import numpy as np
#import cython
import matplotlib.pyplot as plt
import math
import time
from scipy import *
import scipy.signal
import inspect
import sys 
import os


script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, 'Plots_FFA/')
colors = ['green','mediumturquoise','dodgerblue','mediumblue','blueviolet','magenta','crimson']
c=0
fill_value=0

def ffa_code_stage1(data ,dt , T , p_min, p_max, SN_tresh , count_lim, name, cands):
	"""
	ffa_code_stage1 (data , dt , T, period , N , p_min , p_max , medfilt_SN , write_sns, SN_tresh , count_lim , name , cands):
		- data		:  Time series 
		- dt		:  Sampling interval (s)
		- T		:  Total observative time (s)
		- p_min		:  Minimum period in the subset of trial periods (s)
		- p_max		:  Maximum period in the subset of trial periods (s)
		- SN_tresh 	:  S/N treshold when selecting candidates 
		- count_lim	:  int 1 or 2, used in stage2 and stage 3
				   if count_lim =1, goes to 4*dt and 9*dt 
				   if count_lim =2, goes to 8*dt and 27*dt
		- name		:  Name of the beam (without the extension)
		- cands 	:  ffa candidates to be written, belong to class_ffa (see ffa_run.py)

	Returns nothing, adds candidates period with S/N > SN_tresh to cands 

	"""
	# --------------------	   FFA Stage 1	----------------------
	N = T/dt
	w=int(1)
	P0_start, P0_end = np.floor(p_min/dt), np.ceil(p_max/dt)
	P0s = np.arange(P0_start,P0_end,1)
	sigma_total = np.std(data)
	FFA_time1 = time.time()
	print '\n',"Folding for periods from ", p_min, ' to ',p_max, 'sec   with sampling interval',dt,'\n'
	for p0 in P0s:
		p0 = int(p0)
		M_real = float(float(N)/p0)
		added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
		if p0==0 or p0 ==1:
			print 'It tried to fold with period = 0 bin or 1 bin'
			continue
		p0_sec = p0*dt
		xwrap = FFA.XWrap2(data,p0, fill_value=fill_value, pow2=True)
		folds = FFA.FFA(xwrap)
		M = folds.shape[0]

		SN = f.SNR_func(folds, sigma_total, w)
		i = SN >= SN_tresh   
		P = p0 + (np.arange(M, dtype=np.float) / (M-1))
		Psec=P*dt
		if len(Psec[i]) != 0:
			width = [w]*len(Psec[i])
			cands.add_cand( np.around(Psec[i],4), np.around(SN[i],2), dt)
		SNs.extend(SN)
		all_Ps.extend(Psec)

	print '		First FFA is done'

	

#______________________________________________________________________


# --------------------	   FFA Stage 2	------------------------------
# -------------------- Extra downsampling : 2 	----------------------

def ffa_code_stage2(data ,dt ,T, p_min, p_max, SN_tresh, count_lim, name, cands):
	"""
	ffa_code_stage2 (data , dt , T, period , N , p_min , p_max , medfilt_SN , write_sns, SN_tresh , count_lim , name , cands):
		- data		:  Time series 
		- dt		:  Sampling interval (s)
		- T		:  Total observative time (s)
		- p_min		:  Minimum period in the subset of trial periods (s)
		- p_max		:  Maximum period in the subset of trial periods (s)
		- SN_tresh 	:  S/N treshold when selecting candidates 
		- count_lim	:  int 1 or 2, used in stage2 and stage 3
				   if count_lim =1, goes to 4*dt and 9*dt 
				   if count_lim =2, goes to 8*dt and 27*dt
		- name		:  Name of the beam (without the extension)
		- cands 	:  ffa candidates to be written, belong to class_ffa (see ffa_run.py)

	Returns nothing, adds candidates period with S/N > SN_tresh to cands 

	"""	
	
	P0_start, P0_end = np.floor(p_min/dt), np.ceil(p_max/dt)
	P0s = np.arange(P0_start,P0_end,1)
	count=0
	w=int(2)
	while count<=count_lim:
		print '\n','\n','-                Extra downsampling of 2 is being performed'
		if count==0:	
			data_1,data_2 = f.forced_dws_2phase(data)
			new_dt2 = dt*2
			print '\n',"-    	* * *	 	Sampling interval : ",new_dt2*1000," ms	* * * "
			N = T/(new_dt2)
			P0_start, P0_end = np.floor(p_min/new_dt2), np.ceil(p_max/new_dt2)
			P0s2=np.arange(P0_start,P0_end,1)
			sigma_total = np.mean([np.std(data_1),np.std(data_2)])
			print "		Folding ..."
			for p0 in P0s2:
				p0=int(p0)
				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
				if p0==0 or p0 ==1:
					print 'It tried to fold with period = 0 bin or 1 bin','\n'
					continue
				p0_sec = p0*new_dt2
				xwrap_1 = FFA.XWrap2(data_1,p0, fill_value=fill_value, pow2=True)
				xwrap_2 = FFA.XWrap2(data_2,p0, fill_value=fill_value, pow2=True)
		
				folds_1 = FFA.FFA(xwrap_1)
				folds_2 = FFA.FFA(xwrap_2)
					
				SN_1 = f.SNR_func(folds_1, sigma_total,w, added_profs)
				SN_2 = f.SNR_func(folds_2, sigma_total,w, added_profs)

				j = SN_1 >= SN_tresh
				k = SN_2 >= SN_tresh

				M = folds_1.shape[0]
				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt2 
				if len(Psec[j]) != 0:
					width = [w]*len(Psec[j])
 					cands.add_cand( np.around(Psec[j],4), np.around(SN_1[j],2), new_dt2)
				if len(Psec[k]) != 0:
					width = [w]*len(Psec[k])
 					cands.add_cand( np.around(Psec[k],4), np.around(SN_2[k],2), new_dt2)

			

		if count == 1:

			data_11, data_12 = f.forced_dws_2phase(data_1)
			data_21, data_22 = f.forced_dws_2phase(data_2)	
			new_dt2 = new_dt2*2
			w=int(w*2)
			N = T/(new_dt2)
			P0_start, P0_end = np.floor(p_min/new_dt2), np.ceil(p_max/new_dt2)
			P0s2=np.arange(P0_start,P0_end,1)
			sigma_total = np.mean([np.std(data_11),np.std(data_12),np.std(data_21),np.std(data_22)])
			print '\n',"-           * * * 	Sampling interval : ",new_dt2*1000," ms		* * *"
			print "-            Folding ..."
			for p0 in P0s2:
				p0=int(p0)
				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
				if p0==0 or p0 ==1:
					print 'It tried to fold with period = 0 bin or 1 bin'
					continue
				p0_sec = p0*new_dt2
				xwrap_11 = FFA.XWrap2(data_11,p0, fill_value=fill_value, pow2=True)
				xwrap_12 = FFA.XWrap2(data_12,p0, fill_value=fill_value, pow2=True)
				xwrap_21 = FFA.XWrap2(data_21,p0, fill_value=fill_value, pow2=True)
				xwrap_22 = FFA.XWrap2(data_22,p0, fill_value=fill_value, pow2=True)
		
				folds_11 = FFA.FFA(xwrap_11)
				folds_12 = FFA.FFA(xwrap_12)
				folds_21 = FFA.FFA(xwrap_21)
				folds_22 = FFA.FFA(xwrap_22)

				SN_11 = f.SNR_func(folds_11, sigma_total,w, added_profs)
				SN_12 = f.SNR_func(folds_12, sigma_total,w, added_profs)
				SN_21 = f.SNR_func(folds_21, sigma_total,w, added_profs)
				SN_22 = f.SNR_func(folds_22, sigma_total,w, added_profs)

				j = SN_11 >= SN_tresh
				k = SN_21 >= SN_tresh
				l = SN_12 >= SN_tresh
				m = SN_22 >= SN_tresh

				M = folds_11.shape[0]
				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt2 

				if len(Psec[j]) != 0:
					width = [w]*len(Psec[j])
 					cands.add_cand( np.around(Psec[j],4), np.around(SN_11[j],2), new_dt2)
				if len(Psec[k]) != 0:
					width = [w]*len(Psec[k])
 					cands.add_cand( np.around(Psec[k],4), np.around(SN_21[k],2), new_dt2)
				if len(Psec[l]) != 0:
					width = [w]*len(Psec[l])
 					cands.add_cand( np.around(Psec[l],4), np.around(SN_12[l],2), new_dt2)
				if len(Psec[m]) != 0:
					width = [w]*len(Psec[m])
 					cands.add_cand( np.around(Psec[m],4), np.around(SN_22[m],2), new_dt2)



		if count == 2:

			data_111, data_112 = f.forced_dws_2phase(data_11)
			data_121, data_122 = f.forced_dws_2phase(data_12)
			data_211, data_212 = f.forced_dws_2phase(data_21)
			data_221, data_222 = f.forced_dws_2phase(data_22)	
			new_dt2 = new_dt2*2
			w=int(w*2)
			N = T/(new_dt2)

			P0_start, P0_end = np.floor(p_min/new_dt2), np.ceil(p_max/new_dt2)
			P0s2=np.arange(P0_start,P0_end,1)
			sigma_total = np.mean([np.std(data_111),np.std(data_112),np.std(data_211),np.std(data_221),np.std(data_222),np.std(data_212),np.std(data_122),np.std(data_121)])
			print '\n',"*  * * 	Sampling interval : ",new_dt2*1000," ms		* * *"
			print "-            Folding ..."
			for p0 in P0s2:
				p0=int(p0)
				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
				if p0==0 or p0 ==1:
					print 'It tried to fold with period = 0 bin or 1 bin'
					continue
				p0_sec = p0*new_dt2

				xwrap_111 = FFA.XWrap2(data_111,p0, fill_value=fill_value, pow2=True)
				xwrap_121 = FFA.XWrap2(data_121,p0, fill_value=fill_value, pow2=True)
				xwrap_211 = FFA.XWrap2(data_211,p0, fill_value=fill_value, pow2=True)
				xwrap_221 = FFA.XWrap2(data_221,p0, fill_value=fill_value, pow2=True)
				xwrap_112 = FFA.XWrap2(data_112,p0, fill_value=fill_value, pow2=True)
				xwrap_122 = FFA.XWrap2(data_122,p0, fill_value=fill_value, pow2=True)
				xwrap_212 = FFA.XWrap2(data_212,p0, fill_value=fill_value, pow2=True)
				xwrap_222 = FFA.XWrap2(data_222,p0, fill_value=fill_value, pow2=True)
		
				folds_111 = FFA.FFA(xwrap_111)
				folds_121 = FFA.FFA(xwrap_121)
				folds_211 = FFA.FFA(xwrap_211)
				folds_221 = FFA.FFA(xwrap_221)
				folds_112 = FFA.FFA(xwrap_112)
				folds_122 = FFA.FFA(xwrap_122)
				folds_212 = FFA.FFA(xwrap_212)
				folds_222 = FFA.FFA(xwrap_222)


				SN_111 = f.SNR_func(folds_111, sigma_total,w, added_profs)
				SN_121 = f.SNR_func(folds_121, sigma_total,w, added_profs)
				SN_211 = f.SNR_func(folds_211, sigma_total,w, added_profs)
				SN_221 = f.SNR_func(folds_221, sigma_total,w, added_profs)
				SN_112 = f.SNR_func(folds_112, sigma_total,w, added_profs)
				SN_122 = f.SNR_func(folds_122, sigma_total,w, added_profs)
				SN_212 = f.SNR_func(folds_212, sigma_total,w, added_profs)
				SN_222 = f.SNR_func(folds_222, sigma_total,w, added_profs)


				j = SN_111 >= SN_tresh
				k = SN_211 >= SN_tresh
				l = SN_121 >= SN_tresh
				m = SN_221 >= SN_tresh
				n = SN_112 >= SN_tresh
				o = SN_212 >= SN_tresh
				q = SN_122 >= SN_tresh
				r = SN_222 >= SN_tresh

				M = folds_111.shape[0]
				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt2

				if len(Psec[j]) != 0:
					width = [w]*len(Psec[j])
 					cands.add_cand( np.around(Psec[j],4), np.around(SN_111[j],2), new_dt2)
				if len(Psec[k]) != 0:
					width = [w]*len(Psec[k])
 					cands.add_cand( np.around(Psec[k],4), np.around(SN_211[k],2), new_dt2)
				if len(Psec[l]) != 0:
					width = [w]*len(Psec[l])
 					cands.add_cand( np.around(Psec[l],4), np.around(SN_121[l],2), new_dt2)
				if len(Psec[m]) != 0:
					width = [w]*len(Psec[m])
 					cands.add_cand( np.around(Psec[m],4), np.around(SN_221[m],2), new_dt2)
				if len(Psec[n]) != 0:
					width = [w]*len(Psec[n])
 					cands.add_cand( np.around(Psec[n],4), np.around(SN_112[n],2), new_dt2)
				if len(Psec[o]) != 0:
					width = [w]*len(Psec[o])
 					cands.add_cand( np.around(Psec[o],4), np.around(SN_212[o],2), new_dt2)
				if len(Psec[q]) != 0:
					width = [w]*len(Psec[q])
 					cands.add_cand( np.around(Psec[q],4), np.around(SN_122[q],2), new_dt2)
				if len(Psec[r]) != 0:
					width = [w]*len(Psec[r])
 					cands.add_cand( np.around(Psec[r],4), np.around(SN_222[r],2), new_dt2)

		count+=1
	

	
# --------------------	   FFA Stage 3	----------------------
# -------------------- -                Extra downsampling : 3 	----------------------
def ffa_code_stage3(data ,dt ,T, p_min,p_max, SN_tresh,count_lim,name, cands):	
	"""
	ffa_code_stage3 (data , dt , T, period , N , p_min , p_max , medfilt_SN , write_sns, SN_tresh , count_lim , name , cands):
		- data		:  Time series 
		- dt		:  Sampling interval (s)
		- T		:  Total observative time (s)
		- p_min		:  Minimum period in the subset of trial periods (s)
		- p_max		:  Maximum period in the subset of trial periods (s)
		- SN_tresh 	:  S/N treshold when selecting candidates 
		- count_lim	:  int 1 or 2, used in stage2 and stage 3
				   if count_lim =1, goes to 4*dt and 9*dt 
				   if count_lim =2, goes to 8*dt and 27*dt
		- name		:  Name of the beam (without the extension)
		- cands 	:  ffa candidates to be written, belong to class_ffa (see ffa_run.py)

	Returns nothing, adds candidates period with S/N > SN_tresh to cands 

	"""
	P0_start, P0_end = np.floor(p_min/dt), np.ceil(p_max/dt)
	P0s = np.arange(P0_start,P0_end,1)
	count=0
	plot_num=311
	larges_dt3=[]
	w=int(3)
	while count<=count_lim:
		if count==0:	
			print '\n','\n','-                Extra downsampling of 3 is being performed'
			data_1,data_2,data_3=f.forced_dws_3phase(data)
			new_dt3 = dt*3
			print '\n','-              * * * 	Sampling interval : ',new_dt3*1000, ' ms	* * * '
			N = T/(new_dt3)
			P0_start, P0_end = np.floor(p_min/new_dt3), np.ceil(p_max/new_dt3)
			P0s3=np.arange(P0_start,P0_end,1)
			sigma_total = np.mean([np.std(data_1),np.std(data_2),np.std(data_3)])
			print "-            Folding ..."
			for p0 in P0s3:
				p0=int(p0)
				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
				if p0==0 or p0 ==1:
					print 'It tried to fold with period = 0 bin or 1 bin'
					continue
				p0_sec = p0*new_dt3
				xwrap_1 = FFA.XWrap2(data_1,p0, fill_value=fill_value, pow2=True)
				xwrap_2 = FFA.XWrap2(data_2,p0, fill_value=fill_value, pow2=True)
				xwrap_3 = FFA.XWrap2(data_3,p0, fill_value=fill_value, pow2=True)
			
				folds_1 = FFA.FFA(xwrap_1)
				folds_2 = FFA.FFA(xwrap_2)
				folds_3 = FFA.FFA(xwrap_3)

				SN_1 = f.SNR_func(folds_1, sigma_total,w, added_profs)
				SN_2 = f.SNR_func(folds_2, sigma_total,w, added_profs)
				SN_3 = f.SNR_func(folds_3, sigma_total,w, added_profs)

				j = SN_1 >= SN_tresh
				k = SN_2 >= SN_tresh
				l = SN_3 >= SN_tresh
	
				M = folds_1.shape[0]

				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt3 
				if len(Psec[j]) != 0:
					width = [w]*len(Psec[j])
 					cands.add_cand( np.around(Psec[j],4), np.around(SN_1[j],2), new_dt3)
				if len(Psec[k]) != 0:
					width = [w]*len(Psec[k])
 					cands.add_cand( np.around(Psec[k],4), np.around(SN_2[k],2), new_dt3)
				if len(Psec[l]) != 0:
					width = [w]*len(Psec[l])
 					cands.add_cand( np.around(Psec[l],4), np.around(SN_3[l],2), new_dt3)

			count+=1
		


		if count==1:
			print '\n','\n','-                Extra downsampling of 3 is being performed'
			w=int(w*3)	
			data_11, data_12, data_13 = f.forced_dws_3phase(data_1)
			data_21, data_22, data_23 = f.forced_dws_3phase(data_2)
			data_31, data_32, data_33 = f.forced_dws_3phase(data_3)
			new_dt3 = new_dt3*3
			N = T/(new_dt3)
			P0_start, P0_end = np.floor(p_min/new_dt3), np.ceil(p_max/new_dt3)
			P0s3=np.arange(P0_start,P0_end,1)
			sigma_total = np.mean([np.std(data_11),np.std(data_12),np.std(data_13),np.std(data_21),np.std(data_22),np.std(data_23),np.std(data_31),np.std(data_32),np.std(data_33)])
			print '\n','-              * * * 	Sampling interval : ',new_dt3*1000, ' ms	* * *  '
			print "-            Folding ..."
			for p0 in P0s3:
				p0=int(p0)
				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
				if p0==0 or p0 ==1:
					print 'It tried to fold with period = 0 bin or 1 bin'
					continue
				p0_sec = p0*new_dt3
				xwrap_11 = FFA.XWrap2(data_11,p0, fill_value=fill_value, pow2=True)
				xwrap_21 = FFA.XWrap2(data_21,p0, fill_value=fill_value, pow2=True)
				xwrap_31 = FFA.XWrap2(data_31,p0, fill_value=fill_value, pow2=True)
				xwrap_12 = FFA.XWrap2(data_12,p0, fill_value=fill_value, pow2=True)
				xwrap_22 = FFA.XWrap2(data_22,p0, fill_value=fill_value, pow2=True)
				xwrap_32 = FFA.XWrap2(data_32,p0, fill_value=fill_value, pow2=True)
				xwrap_13 = FFA.XWrap2(data_13,p0, fill_value=fill_value, pow2=True)
				xwrap_23 = FFA.XWrap2(data_23,p0, fill_value=fill_value, pow2=True)
				xwrap_33 = FFA.XWrap2(data_33,p0, fill_value=fill_value, pow2=True)
			
				folds_11 = FFA.FFA(xwrap_11)
				folds_21 = FFA.FFA(xwrap_21)
				folds_31 = FFA.FFA(xwrap_31)
				folds_12 = FFA.FFA(xwrap_12)
				folds_22 = FFA.FFA(xwrap_22)
				folds_32 = FFA.FFA(xwrap_32)
				folds_13 = FFA.FFA(xwrap_13)
				folds_23 = FFA.FFA(xwrap_23)
				folds_33 = FFA.FFA(xwrap_33)

				SN_11 = f.SNR_func(folds_11, sigma_total,w, added_profs)
				SN_21 = f.SNR_func(folds_21, sigma_total,w, added_profs)
				SN_31 = f.SNR_func(folds_31, sigma_total,w, added_profs)
				SN_12 = f.SNR_func(folds_12, sigma_total,w, added_profs)
				SN_22 = f.SNR_func(folds_22, sigma_total,w, added_profs)
				SN_32 = f.SNR_func(folds_32, sigma_total,w, added_profs)
				SN_13 = f.SNR_func(folds_13, sigma_total,w, added_profs)
				SN_23 = f.SNR_func(folds_23, sigma_total,w, added_profs)
				SN_33 = f.SNR_func(folds_33, sigma_total,w, added_profs)

				j = SN_11 >= SN_tresh
				k = SN_21 >= SN_tresh
				l = SN_31 >= SN_tresh
				m = SN_12 >= SN_tresh
				n = SN_22 >= SN_tresh
				o = SN_32 >= SN_tresh
				q = SN_13 >= SN_tresh
				r = SN_23 >= SN_tresh
				s = SN_33 >= SN_tresh

				M = folds_11.shape[0]
				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt3 

				if len(Psec[j]) != 0:
					width = [w]*len(Psec[j])
					cands.add_cand( np.around(Psec[j],4), np.around(SN_11[j],2), new_dt3)
				if len(Psec[k]) != 0:
					width = [w]*len(Psec[k])
 					cands.add_cand( np.around(Psec[k],4), np.around(SN_21[k],2), new_dt3)
				if len(Psec[l]) != 0:
					width = [w]*len(Psec[l])
 					cands.add_cand( np.around(Psec[l],4), np.around(SN_31[l],2), new_dt3)
				if len(Psec[m]) != 0:
					width = [w]*len(Psec[m])
 					cands.add_cand( np.around(Psec[m],4), np.around(SN_12[m],2),new_dt3)
				if len(Psec[n]) != 0:
					width = [w]*len(Psec[n])
 					cands.add_cand( np.around(Psec[n],4), np.around(SN_22[n],2),new_dt3)
				if len(Psec[o]) != 0:
					width = [w]*len(Psec[o])
 					cands.add_cand( np.around(Psec[o],4), np.around(SN_32[o],2),new_dt3)
				if len(Psec[q]) != 0:
					width = [w]*len(Psec[q])
 					cands.add_cand( np.around(Psec[q],4), np.around(SN_13[q],2),new_dt3)
				if len(Psec[r]) != 0:
					width = [w]*len(Psec[r])
 					cands.add_cand( np.around(Psec[r],4), np.around(SN_23[r],2),new_dt3)
				if len(Psec[s]) != 0:
					width = [w]*len(Psec[s])
 					cands.add_cand( np.around(Psec[s],4), np.around(SN_33[s],2),new_dt3)

			count+=1
			# ----------------     Periodogram plot  --------------------

		
		
		if count==2 :
			print '\n','\n','-                Extra downsampling of 3 is being performed'
			w=int(w*3)	
			data_111, data_112, data_113 = f.forced_dws_3phase(data_11)
			data_211, data_212, data_213 = f.forced_dws_3phase(data_21)
			data_311, data_312, data_313 = f.forced_dws_3phase(data_31)
			data_121, data_122, data_123 = f.forced_dws_3phase(data_12)
			data_221, data_222, data_223 = f.forced_dws_3phase(data_22)
			data_321, data_322, data_323 = f.forced_dws_3phase(data_32)
			data_131, data_132, data_133 = f.forced_dws_3phase(data_13)
			data_231, data_232, data_233 = f.forced_dws_3phase(data_23)
			data_331, data_332, data_333 = f.forced_dws_3phase(data_33)
			
			new_dt3 = new_dt3*3
			print '\n','-              * * * 	Sampling interval : ',new_dt3*1000, ' ms	* * *  '
			N = T/(new_dt3)
			P0_start, P0_end = np.floor(p_min/new_dt3), np.ceil(p_max/new_dt3)
			P0s3=np.arange(P0_start,P0_end,1)
			sigma_total = np.mean([np.std(data_111),np.std(data_112),np.std(data_113),np.std(data_121),np.std(data_122),np.std(data_123),np.std(data_131),np.std(data_132),np.std(data_133),np.std(data_211),np.std(data_212),np.std(data_213),np.std(data_221),np.std(data_222),np.std(data_223),np.std(data_231),np.std(data_232),np.std(data_233),np.std(data_311),np.std(data_312),np.std(data_313),np.std(data_321),np.std(data_322),np.std(data_323),np.std(data_331),np.std(data_332),np.std(data_333)])

			print '\n',"-            Folding ..."
			for p0 in P0s3:
				p0=int(p0)
				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
				if p0==0 or p0 ==1:
					print 'It tried to fold with period = 0 bin or 1 bin'
					continue
				p0_sec = p0*new_dt3
				xwrap_111 = FFA.XWrap2(data_111,p0, fill_value=fill_value, pow2=True)
				xwrap_112 = FFA.XWrap2(data_112,p0, fill_value=fill_value, pow2=True)
				xwrap_113 = FFA.XWrap2(data_113,p0, fill_value=fill_value, pow2=True)	
				xwrap_121 = FFA.XWrap2(data_121,p0, fill_value=fill_value, pow2=True)
				xwrap_122 = FFA.XWrap2(data_122,p0, fill_value=fill_value, pow2=True)
				xwrap_123 = FFA.XWrap2(data_123,p0, fill_value=fill_value, pow2=True)
				xwrap_131 = FFA.XWrap2(data_131,p0, fill_value=fill_value, pow2=True)
				xwrap_132 = FFA.XWrap2(data_132,p0, fill_value=fill_value, pow2=True)
				xwrap_133 = FFA.XWrap2(data_133,p0, fill_value=fill_value, pow2=True)
				xwrap_211 = FFA.XWrap2(data_211,p0, fill_value=fill_value, pow2=True)
				xwrap_212 = FFA.XWrap2(data_212,p0, fill_value=fill_value, pow2=True)
				xwrap_213 = FFA.XWrap2(data_213,p0, fill_value=fill_value, pow2=True)
				xwrap_221 = FFA.XWrap2(data_221,p0, fill_value=fill_value, pow2=True)
				xwrap_222 = FFA.XWrap2(data_222,p0, fill_value=fill_value, pow2=True)
				xwrap_223 = FFA.XWrap2(data_223,p0, fill_value=fill_value, pow2=True)
				xwrap_231 = FFA.XWrap2(data_231,p0, fill_value=fill_value, pow2=True)
				xwrap_232 = FFA.XWrap2(data_232,p0, fill_value=fill_value, pow2=True)
				xwrap_233 = FFA.XWrap2(data_233,p0, fill_value=fill_value, pow2=True)
				xwrap_311 = FFA.XWrap2(data_311,p0, fill_value=fill_value, pow2=True)
				xwrap_312 = FFA.XWrap2(data_312,p0, fill_value=fill_value, pow2=True)
				xwrap_313 = FFA.XWrap2(data_313,p0, fill_value=fill_value, pow2=True)
				xwrap_321 = FFA.XWrap2(data_321,p0, fill_value=fill_value, pow2=True)
				xwrap_322 = FFA.XWrap2(data_322,p0, fill_value=fill_value, pow2=True)
				xwrap_323 = FFA.XWrap2(data_323,p0, fill_value=fill_value, pow2=True)
				xwrap_331 = FFA.XWrap2(data_331,p0, fill_value=fill_value, pow2=True)
				xwrap_332 = FFA.XWrap2(data_332,p0, fill_value=fill_value, pow2=True)
				xwrap_333 = FFA.XWrap2(data_333,p0, fill_value=fill_value, pow2=True)
				
				folds_111 = FFA.FFA(xwrap_111)
				folds_112 = FFA.FFA(xwrap_112)
				folds_113 = FFA.FFA(xwrap_113)
				folds_211 = FFA.FFA(xwrap_211)
				folds_212 = FFA.FFA(xwrap_212)
				folds_213 = FFA.FFA(xwrap_213)
				folds_311 = FFA.FFA(xwrap_311)
				folds_312 = FFA.FFA(xwrap_312)
				folds_313 = FFA.FFA(xwrap_313)
				folds_121 = FFA.FFA(xwrap_121)
				folds_122 = FFA.FFA(xwrap_122)
				folds_123 = FFA.FFA(xwrap_123)
				folds_221 = FFA.FFA(xwrap_221)
				folds_222 = FFA.FFA(xwrap_222)
				folds_223 = FFA.FFA(xwrap_223)
				folds_321 = FFA.FFA(xwrap_321)
				folds_322 = FFA.FFA(xwrap_322)
				folds_323 = FFA.FFA(xwrap_323)
				folds_131 = FFA.FFA(xwrap_131)
				folds_132 = FFA.FFA(xwrap_132)
				folds_133 = FFA.FFA(xwrap_133)
				folds_231 = FFA.FFA(xwrap_231)
				folds_232 = FFA.FFA(xwrap_232)
				folds_233 = FFA.FFA(xwrap_233)
				folds_331 = FFA.FFA(xwrap_331)
				folds_332 = FFA.FFA(xwrap_332)
				folds_333 = FFA.FFA(xwrap_333)
			

				SN_111 = f.SNR_func(folds_111, sigma_total,w, added_profs)
				SN_211 = f.SNR_func(folds_211, sigma_total,w, added_profs)
				SN_311 = f.SNR_func(folds_311, sigma_total,w, added_profs)
				SN_121 = f.SNR_func(folds_121, sigma_total,w, added_profs)
				SN_221 = f.SNR_func(folds_221, sigma_total,w, added_profs)
				SN_321 = f.SNR_func(folds_321, sigma_total,w, added_profs)
				SN_131 = f.SNR_func(folds_131, sigma_total,w, added_profs)
				SN_231 = f.SNR_func(folds_231, sigma_total,w, added_profs)
				SN_331 = f.SNR_func(folds_331, sigma_total,w, added_profs)
				SN_112 = f.SNR_func(folds_112, sigma_total,w, added_profs)
				SN_212 = f.SNR_func(folds_212, sigma_total,w, added_profs)
				SN_312 = f.SNR_func(folds_312, sigma_total,w, added_profs)
				SN_122 = f.SNR_func(folds_122, sigma_total,w, added_profs)
				SN_222 = f.SNR_func(folds_222, sigma_total,w, added_profs)
				SN_322 = f.SNR_func(folds_322, sigma_total,w, added_profs)
				SN_132 = f.SNR_func(folds_132, sigma_total,w, added_profs)
				SN_232 = f.SNR_func(folds_232, sigma_total,w, added_profs)
				SN_332 = f.SNR_func(folds_332, sigma_total,w, added_profs)
				SN_113 = f.SNR_func(folds_113, sigma_total,w, added_profs)
				SN_213 = f.SNR_func(folds_213, sigma_total,w, added_profs)
				SN_313 = f.SNR_func(folds_313, sigma_total,w, added_profs)
				SN_123 = f.SNR_func(folds_123, sigma_total,w, added_profs)
				SN_223 = f.SNR_func(folds_223, sigma_total,w, added_profs)
				SN_323 = f.SNR_func(folds_323, sigma_total,w, added_profs)
				SN_133 = f.SNR_func(folds_133, sigma_total,w, added_profs)
				SN_233 = f.SNR_func(folds_233, sigma_total,w, added_profs)
				SN_333 = f.SNR_func(folds_333, sigma_total,w, added_profs)


				j = SN_111 >= SN_tresh
				k = SN_211 >= SN_tresh
				l = SN_311 >= SN_tresh
				m = SN_121 >= SN_tresh
				n = SN_221 >= SN_tresh
				o = SN_321 >= SN_tresh
				q = SN_131 >= SN_tresh
				r = SN_231 >= SN_tresh
				s = SN_331 >= SN_tresh
				jj = SN_112 >= SN_tresh
				kk = SN_212 >= SN_tresh
				ll = SN_312 >= SN_tresh
				mm = SN_122 >= SN_tresh
				nn = SN_222 >= SN_tresh
				oo = SN_322 >= SN_tresh
				qq = SN_132 >= SN_tresh
				rr = SN_232 >= SN_tresh
				ss = SN_332 >= SN_tresh
				jjj = SN_113 >= SN_tresh
				kkk = SN_213 >= SN_tresh
				lll = SN_313 >= SN_tresh
				mmm = SN_123 >= SN_tresh
				nnn = SN_223 >= SN_tresh
				ooo = SN_323 >= SN_tresh
				qqq = SN_133 >= SN_tresh
				rrr = SN_233 >= SN_tresh
				sss = SN_333 >= SN_tresh		

				M = folds_111.shape[0]
				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt3 

				if len(Psec[j]) != 0:
					width = [w]*len(Psec[j])
 					cands.add_cand( np.around(Psec[j],4), np.around(SN_111[j],2),new_dt3)
				if len(Psec[k]) != 0:
					width = [w]*len(Psec[k])
 					cands.add_cand( np.around(Psec[k],4), np.around(SN_211[k],2),new_dt3)
				if len(Psec[l]) != 0:
					width = [w]*len(Psec[l])
 					cands.add_cand( np.around(Psec[l],4), np.around(SN_311[l],2),new_dt3)
				if len(Psec[m]) != 0:
					width = [w]*len(Psec[m])
 					cands.add_cand( np.around(Psec[m],4), np.around(SN_121[m],2),new_dt3)
				if len(Psec[n]) != 0:
					width = [w]*len(Psec[n])
 					cands.add_cand( np.around(Psec[n],4), np.around(SN_221[n],2),new_dt3)
				if len(Psec[o]) != 0:
					width = [w]*len(Psec[o])
 					cands.add_cand( np.around(Psec[o],4), np.around(SN_321[o],2),new_dt3)
				if len(Psec[q]) != 0:
					width = [w]*len(Psec[q])
 					cands.add_cand( np.around(Psec[q],4), np.around(SN_131[q],2),new_dt3)
				if len(Psec[r]) != 0:
					width = [w]*len(Psec[r])
 					cands.add_cand( np.around(Psec[r],4), np.around(SN_231[r],2),new_dt3)
				if len(Psec[s]) != 0:
					width = [w]*len(Psec[s])
 					cands.add_cand( np.around(Psec[s],4), np.around(SN_331[s],2),new_dt3)
				if len(Psec[jj]) != 0:
					width = [w]*len(Psec[jj])
 					cands.add_cand( np.around(Psec[jj],4), np.around(SN_112[jj],2),new_dt3)
				if len(Psec[kk]) != 0:
					width = [w]*len(Psec[kk])
 					cands.add_cand( np.around(Psec[kk],4), np.around(SN_212[kk],2),new_dt3)
				if len(Psec[ll]) != 0:
					width = [w]*len(Psec[ll])
 					cands.add_cand( np.around(Psec[ll],4), np.around(SN_312[ll],2),new_dt3)
				if len(Psec[mm]) != 0:
					width = [w]*len(Psec[mm])
 					cands.add_cand( np.around(Psec[mm],4), np.around(SN_122[mm],2),new_dt3)
				if len(Psec[nn]) != 0:
					width = [w]*len(Psec[nn])
 					cands.add_cand( np.around(Psec[nn],4), np.around(SN_222[nn],2),new_dt3)
				if len(Psec[oo]) != 0:
					width = [w]*len(Psec[oo])
 					cands.add_cand( np.around(Psec[oo],4), np.around(SN_322[oo],2),new_dt3)
				if len(Psec[qq]) != 0:
					width = [w]*len(Psec[qq])
 					cands.add_cand( np.around(Psec[qq],4), np.around(SN_132[qq],2),new_dt3)
				if len(Psec[rr]) != 0:
					width = [w]*len(Psec[rr])
 					cands.add_cand( np.around(Psec[rr],4), np.around(SN_232[rr],2),new_dt3)
				if len(Psec[ss]) != 0:
					width = [w]*len(Psec[ss])
 					cands.add_cand( np.around(Psec[ss],4), np.around(SN_332[ss],2),new_dt3)
				if len(Psec[jjj]) != 0:
					width = [w]*len(Psec[jjj])
 					cands.add_cand( np.around(Psec[jjj],4), np.around(SN_113[jjj],2),new_dt3)
				if len(Psec[kkk]) != 0:
					width = [w]*len(Psec[kkk])
 					cands.add_cand( np.around(Psec[kkk],4), np.around(SN_213[kkk],2),new_dt3)
				if len(Psec[lll]) != 0:
					width = [w]*len(Psec[lll])
 					cands.add_cand( np.around(Psec[lll],4), np.around(SN_313[lll],2),new_dt3)
				if len(Psec[mmm]) != 0:
					width = [w]*len(Psec[mmm])
 					cands.add_cand( np.around(Psec[mmm],4), np.around(SN_123[mmm],2),new_dt3)
				if len(Psec[nnn]) != 0:
					width = [w]*len(Psec[nnn])
 					cands.add_cand( np.around(Psec[nnn],4), np.around(SN_223[nnn],2),new_dt3)
				if len(Psec[ooo]) != 0:
					width = [w]*len(Psec[ooo])
 					cands.add_cand( np.around(Psec[ooo],4), np.around(SN_323[ooo],2),new_dt3)
				if len(Psec[qqq]) != 0:
					width = [w]*len(Psec[qqq])
 					cands.add_cand( np.around(Psec[qqq],4), np.around(SN_133[qqq],2),new_dt3)
				if len(Psec[rrr]) != 0:
					width = [w]*len(Psec[rrr])
 					cands.add_cand( np.around(Psec[rrr],4), np.around(SN_233[rrr],2),new_dt3)
				if len(Psec[sss]) != 0:
					width = [w]*len(Psec[sss])
 					cands.add_cand( np.around(Psec[sss],4), np.around(SN_333[sss],2),new_dt3)

			count+=1


