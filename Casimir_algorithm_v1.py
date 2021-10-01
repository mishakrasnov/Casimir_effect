import numpy as np
import matplotlib
import mpmath as mp
import sympy as sp
import plot_graph
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import casimir_m_s_p as msp
from consts import *

def find_zeros(ar_y,ar_x,epsilon,zero_separation,zeros_count = 10,label = ''):
	if ar_y.shape != ar_x.shape:
		print("error in find_zeros")
		return
	result = []
	vals = np.abs(ar_y)
	j = 0
	for i in range(ar_x.size):
		if len(result) == zeros_count:
			return result	
		if vals[i] < epsilon and j==0:
			result.append(ar_x[i])
			j = zero_separation
		elif j:
			j-=1	
	if(result == []):	 
		return [ar_x[-1]]
	return result
def find_max(Vals,zeros_count):
	Result = np.zeros((Vals.shape[0],zeros_count-1))
	for j in range(Vals.shape[0]):
		row = Temp[j]
		result = []
		for i in range(row.size-2):
			if row[i+1]-row[i] > 0 and row[i+2]-row[i+1] < 0:
				result.append(np.maximum(np.abs(row[i]),0))
			if row[i+1]-row[i] < 0 and row[i+2]-row[i+1] > 0:
				result.append(np.maximum(np.abs(row[i]),0))
		Result[j] = np.array(result)
	Result = np.c_[Temp[:,0],Result]
	return Result

def eps_l(w,f,w0,const):
	w = w*h
	return const + (f[0]*w0[0]**2/(w0[0]**2 + w**2 +w)+ 
                       f[1]*w0[1]**2/(w0[1]**2 + w**2 +w)  )
def eps_d(w):
	return np.exp(4)

def eps_sub(w):
    w0_l_sub = [0.1,sp.symbols('w0')]
    f_l_sub = [0.7,0.3]
    return eps_l(w,f_l_sub,w0_l_sub,1.6)
def eps_med(w):
    w0_l_med= [0.05,2]
    f_l_med = [0.8,0.29]
    return eps_l(w,f_l_med,w0_l_med,1.7)

et_freqs =  np.array([1.62/1000,4.32/1000,1.12/10,6.87,15.2,15.6,43.8])
et_strength =  np.array([9.57/10,1.62,1.4/10,1.26/10,4.16/10,2.44/10,7.1/100])

def eps_et(w):
	w = w*h
	eps = 1.5
	for i in range(7):
		eps += (et_strength[i])/(1 + (w/et_freqs[i])**2)
	return eps
	
solver = msp.Casimir_m_s_p([eps_l_sub,eps_l_med,eps_d])

Deviration = solver.get_deviration('numpy',sp.symbols('w0'))

def Dev(w,k,d,w0):
	return Deviration([w/h,k/t_nm_m,d*t_nm_m,w0])

W = np.array([ 10**(-30) + i/10 for i in range(2000)],dtype=np.cdouble)
D = np.array(list(range(1,100)),dtype = np.cdouble)
W0 = np.array([5 + i for i in range(10)])
All_DATA = np.zeros((W0.size,D.size,W.size))
ALL_ZEROS = []
ALL_MAX = []
ALL_DELTA_ZEROS = []
W_for_eps = np.exp(np.array([-8 + i/50 for i in range(600)]))
Integrlals = np.zeros((W0.size,D.size))
W_,D_ = np.meshgrid(W,D,indexing = 'xy')

w0_l_sub = [0.1,8]
plt.loglog(W_for_eps,eps_l_sub(W_for_eps/h))
plt.loglog(W_for_eps,eps_l_medium(W_for_eps/h))
plt.show()


def eps_diff(w):
	return eps_l_sub(w)-eps_l_medium(w) 

def elmentary_transform(A):
	Temp = np.zeros(A.shape)
	Temp[:,0] = A[:,0]
	for i in range(1,A.shape[1]):
		Temp[:,i] = (-1)**(i)*(A[:,i]-A[:,i-1])	
	return Temp

for j in range(W0.size):
	w0_l_sub[0] = W0[j]
	try:
		vals = eps_diff(W_for_eps/h)
		zeros = find_zeros(vals,W_for_eps,0.05,20) + [0]
		print(zeros)
		Zeros = np.array(zeros*D.size).reshape(D.size,-1)
	except Exception as e:
		raise e
	else:
		Temp = All_DATA[j] = np.real(Dev(W_,0.001,D_,W0[j])) # Temp[i][j] = Dev(W[j],D[i])
		Max = find_max(Temp,Zeros.shape[1])
		ALL_MAX.append(Max)
		for i in range(D.size):
			H = Max[i][-1]/4
			data = Temp[i,W > Zeros[i,-2]+10]
			Zeros[i][-1]=np.real(find_zeros(data,W[W > Zeros[i,-2]+10],H,20,1,str(i))[0])
		ALL_ZEROS.append(Zeros)
		Zeros = elmentary_transform(Zeros)
		ALL_DELTA_ZEROS.append(Zeros)
		Integrlal = np.sum(Zeros*Max,axis = 1).reshape(-1,D.size)
		Integrlals[j] = Integrlal


BRIGHT = 10**(20)
Mask = np.zeros(Integrlals.shape)
enter = input('I\'m ready\n')

def masking(M,DATA,cond):
	try:
		for j in range(M.shape[0]):
			row = DATA[j]
			row_ = M[j]
			temp = []
			for i in range(1,row.size-1):
				if eval(cond):
					temp.append(i)
			row_[temp] = np.array([BRIGHT]*len(temp))
	except Exception as e:
		print(e)	

while enter != 'exit':
	if enter == 'i_t':
		enter = input('...')
		while enter != 'up':
			if enter == 'mask':
				enter = input('...')
				Mask[:,:] = 0
				masking(Mask,Integrlals,enter)
			else:	
				try:
					(i,j),(k,e)= eval(enter)
					plt.matshow(np.maximum(Mask[i:j,k:e],np.abs(Integrlals[i:j,k:e])))
					plt.title('W0 = [' + str(np.round_(W0[0],2)) +','+str(np.round_(W0[-1],2)) + ']')
				except Exception as e:
					print(e)
				else:
					plt.colorbar()
					plt.show()
			enter = input('...')		
	elif enter == 'f_g':
		enter = input('...')
		while enter != 'up':
			try:	
				w0,d = eval(enter)
				plt.plot(np.real(W),All_DATA[w0][d])
				print('...zeros: ',ALL_ZEROS[w0][d])
				print('...maxs: ',ALL_MAX[w0][d])
				print('...delta_zeros: ',ALL_DELTA_ZEROS[w0][d])
			except:
				print('...not 2*int')
			else:
				plt.show()
			enter = input('...')
	enter = input()