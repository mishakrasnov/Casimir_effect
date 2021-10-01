import numpy as np
import matplotlib
import mpmath as mp
import sympy as sp
import plot_graph
from consts import *
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import stp as tm

class Casimir_m_s_p:
	def __init__(self,args):
		self.eps_1,self.eps_2,self.eps_3 = args #epsilons: s,m,p

	def get_f(self,lib,*params):
		e = sp.symbols('e',positive = 1)
		k = sp.symbols('k',positive = 1)
		d = sp.symbols('d',positive = 1)
		eps = [self.eps_1,self.eps_2,self.eps_3,self.eps_2]
		n1 = 0	
		n2 = 1
		name1 = 'sub'
		name2 = 'plate'
		solver1 = tm.Transfer_matrix(n1,name1)
		solver2 = tm.Transfer_matrix(n2,name2)
		Coeffs1 = solver1.get_coeffs(1j*e,k)
		Coeffs2 = solver2.get_coeffs(1j*e,k)
		r_s1 = Coeffs1[1]
		r_p1 = Coeffs1[3]
		r_s2 = Coeffs2[1]
		r_p2 = Coeffs2[3]
		f_s = (1 - r_s1*r_s2*sp.exp( -2*d*sp.sqrt((e/c)**2*self.eps_2(e) + k**2)))
		f_p = (1 - r_p2*r_p1*sp.exp( -2*d*sp.sqrt((e/c)**2*self.eps_2(e) + k**2)))

		f = ((sp.log(f_s)+sp.log(f_p))*k)

		args = []
		for i in range(n1+2):
			args.append((solver1.S_permit[i],eps[i](e)))
			args.append((solver1.S_permet[i],1))
		for i in range(n2+2):
			args.append((solver2.S_permit[i],eps[i+1](e)))
			args.append((solver2.S_permet[i],1))
		args.append((solver2.S_length[0],5*10**(-9)))
		f = sp.lambdify([(e,k,d,*params)],f.subs(args),lib)
		return f

	def get_deviration(self,lib,*params):
		e = sp.symbols('e',positive = 1)
		k = sp.symbols('k',positive = 1)
		d = sp.symbols('d',positive = 1)
		eps = [self.eps_1,self.eps_2,self.eps_3,self.eps_2]
		n1 = 0	
		n2 = 1
		name1 = 'sub'
		name2 = 'plate'
		solver1 = tm.Transfer_matrix(n1,name1)
		solver2 = tm.Transfer_matrix(n2,name2)
		Coeffs1 = solver1.get_coeffs(1j*e,k)
		Coeffs2 = solver2.get_coeffs(1j*e,k)
		r_s1 = Coeffs1[1]
		r_p1 = Coeffs1[3]
		r_s2 = Coeffs2[1]
		r_p2 = Coeffs2[3]
		kp = sp.sqrt(k**2+(e/c)**2*self.eps_2(e))
		f_s = k*kp*r_s1*r_s2*sp.exp(-2*d*kp)/(1-r_s1*r_s2*sp.exp(-2*d*kp))
		f_p = k*kp*r_p1*r_p2*sp.exp(-2*d*kp)/(1-r_p1*r_p2*sp.exp(-2*d*kp))

		f = f_s + f_p

		args = []
		for i in range(n1+2):
			args.append((solver1.S_permit[i],eps[i](e)))
			args.append((solver1.S_permet[i],1))
		for i in range(n2+2):
			args.append((solver2.S_permit[i],eps[i+1](e)))
			args.append((solver2.S_permet[i],1))
		args.append((solver2.S_length[0],5*10**(-9)))
		f = sp.lambdify([(e,k,d,*params)],f.subs(args),lib)
		return f
	
