
import sympy as sp
class Transfer_matrix:

	def get_list_of_symbols(self,N,name):
		name = self.system_name + name
		return [sp.symbols(name+str(i),complex =1) for i in range(N)] 

	def get_3d_list_of_symbols(self,N):

		return [sp.zeros(2,2,complex =1) for i in range(N)]
	
	def __init__(self,N,system_name):

		self.system_name = system_name
		
		c = 3*10**(8)

		self.N = N

		self.S_t_matrixs_s = self.get_3d_list_of_symbols(2*N+1) 

		self.S_t_matrixs_p = self.get_3d_list_of_symbols(2*N+1)
 
		self.S_wave_vectors = sp.zeros(1,N+2,complex =1)
 
		self.S_k_z = sp.zeros(1,N+2,complex =1 )

		self.S_permit = self.get_list_of_symbols(N+2,'permit')
		
		self.S_permet = self.get_list_of_symbols(N+2,'permet')

		self.S_length = self.get_list_of_symbols(N,'length')

	def get_r_s_bound(self,k_in_z, k_out_z, p_in, p_out):

		return sp.simplify((p_out*k_in_z - p_in*k_out_z)/(p_out*k_in_z + p_in*k_out_z))

	def get_t_s_bound(self,k_in_z, k_out_z, p_in, p_out):

		return sp.simplify((2*p_out*k_in_z)/(p_out*k_in_z + p_in*k_out_z))

	def get_r_p_bound(self,k_in_z, k_out_z, e_in, e_out):

		return sp.simplify((e_out*k_in_z - e_in*k_out_z)/(e_out*k_in_z + e_in*k_out_z))

	def get_t_p_bound(self,k_in_z, k_out_z, e_in, e_out,p_in,p_out):

		return sp.simplify(2*(e_out*k_in_z)/(e_out*k_in_z + e_in*k_out_z)*sp.sqrt(p_out*e_in/(p_in*e_out)))

	def get_matrix_bound(self,t_in,t_out,r_in,r_out):

		elements = [[t_in*t_out - r_in*r_out, r_out],[(-r_in),1]]
		T = sp.Matrix(elements,complex = 1) 
		return T/t_out

	def get_matrix_slab(self,length,k_z):

		elements = [[sp.exp(length*k_z*1j),0],[0,sp.exp(-length*k_z*1j)]]
		T = sp.Matrix(elements,complex =1) 
		T = sp.simplify(T)
		return T

	def single_det(self,T):
		return(T.row(0)[0]*T.row(1)[1] - T.row(1)[0]*T.row(0)[1])
	def get_multy_det(self,T):   #det(A*B) = det(A)*det(B)
		
		N = self.N
		result = 1
		for i in [j for j in range(2*N+1) if j%2 == 0]: # must be j%2 == 0
			result = result*self.single_det(T[i])

		return result

	def get_t_from_matrix(self,T,T_full):
		return self.get_multy_det(T_full)/T.row(1)[1]

	def get_r_from_matrix(self,T):
		return -T.row(1)[0]/T.row(1)[1]

	def get_coeffs(self,S_w,S_k_x):
		
		N = self.N
		
		S_k = S_w/(3*10**(8))
		
		result = sp.zeros(1,4,complex =1) 

		for i in range(N+2):

			self.S_wave_vectors[i] = sp.simplify(sp.sqrt(self.S_permit[i]*self.S_permet[i])*S_k)

		for i in range(N+2):

			self.S_k_z[i] = sp.simplify(sp.sqrt(self.S_wave_vectors[i]**2 - S_k_x**2))

		for i in range(N+1):	

			t_in_s = self.get_t_s_bound(self.S_k_z[i],self.S_k_z[i+1],self.S_permet[i],self.S_permet[i+1])
			r_in_s = self.get_r_s_bound(self.S_k_z[i],self.S_k_z[i+1],self.S_permet[i],self.S_permet[i+1])
			t_out_s = self.get_t_s_bound(self.S_k_z[i+1],self.S_k_z[i],self.S_permet[i+1],self.S_permet[i])
			r_out_s = self.get_r_s_bound(self.S_k_z[i+1],self.S_k_z[i],self.S_permet[i+1],self.S_permet[i])
			
			t_in_p = self.get_t_p_bound(self.S_k_z[i],self.S_k_z[i+1], self.S_permit[i],self.S_permit[i+1],
																		self.S_permet[i],self.S_permet[i+1])

			r_in_p = self.get_r_p_bound(self.S_k_z[i],self.S_k_z[i+1],self.S_permit[i],self.S_permit[i+1])
			
			t_out_p = self.get_t_p_bound(self.S_k_z[i+1],self.S_k_z[i],self.S_permit[i+1],self.S_permit[i],
																		self.S_permet[i+1],self.S_permet[i])

			r_out_p = self.get_r_p_bound(self.S_k_z[i+1],self.S_k_z[i],self.S_permit[i+1],self.S_permit[i])

			self.S_t_matrixs_s[2*i] = self.get_matrix_bound(t_in_s,t_out_s,r_in_s,r_out_s)
			self.S_t_matrixs_p[2*i] = self.get_matrix_bound(t_in_p,t_out_p,r_in_p,r_out_p)

			if(i < N):

				self.S_t_matrixs_s[2*i+1] = self.get_matrix_slab(self.S_length[i],self.S_k_z[i+1])
				self.S_t_matrixs_p[2*i+1] = self.get_matrix_slab(self.S_length[i],self.S_k_z[i+1])
		
		S_t_result_s = self.S_t_matrixs_s[2*N] # the result S-matrix for all structure
		S_t_result_p = self.S_t_matrixs_p[2*N]
		
		for i in range(1,2*N+1):	
			S_t_result_s = S_t_result_s*self.S_t_matrixs_s[2*N-i]
			S_t_result_p = S_t_result_p*self.S_t_matrixs_p[2*N-i]

		result[0] = self.get_t_from_matrix(S_t_result_s,self.S_t_matrixs_s)		
		result[1] = self.get_r_from_matrix(S_t_result_s)
		result[2] = self.get_t_from_matrix(S_t_result_p,self.S_t_matrixs_p)		
		result[3] = self.get_r_from_matrix(S_t_result_p)

		return result
