
import toml
import os
import sys
#import subprocess
dirname = os.path.dirname(__file__)
import seir
#sys.path.append(dirname+"/../../KaXim-Tutorial/script")
#import plotting_kappa as kappa
if(dirname+"/../../lib" not in sys.path):
	sys.path.insert(1,dirname+"/../../lib")
#import KaXimDebug as KaXim
import KaXim as KaXim
import numpy as np
import json
import toml
import tempfile

#applies a factor to n dependening if it is a percent(float) or quantity(int)
def p_or_q(n,base_p,base_q=1):
	if isinstance(n,float) :	#is quantity
		return n*base_p
	elif isinstance(n,int):		#is percent
		return n*base_q
	else:
		raise Exception("p_or_q: not a number.")

#def call_wparams(**kwargs):

#returns a fraction init value using global or local fraction
#fractions can be quantities or percents
#local fraction used only if global is 0
def init_frac(fglob,floc,total_glob,total_loc):
	if n > 1:		#glob is quantity -> 
		return fglob*total_loc / total_glob
	elif n:
		return fglob*total_loc
	else:
		return p_or_q(floc,total_loc)
		
is_frac = lambda val : val >= 0.0 and val <= 1.0
is_rate = lambda val : val >= 0.0
is_count = lambda val : val >= 0

class Param:
	def __init__(self,name,dflt = None,deps = None,dyn = False,func = None,is_valid = None):
		self.name = name			#Name of the parameter
		self.dflt = dflt			#default value (None if mandatory)
		self.deps = deps			#dependencies from toml params (sorted by prevalence, space-separated) 
		self.dyn = dyn				#True if dynamic parameter
		self.func = func			#value = func(v = toml[param_name],p = param_name, d = values_dict)
		self.is_valid = is_valid	#if not checks(value) then bad-parameter-definition. 
		

class RBM_SEIR(seir.SEIR):

	STATES = ["S","E","I","R"]
	STATES_NAMES = ["Susceptible","Exposed","Infected","Removed"]
	
	PARAMETERS = [
#Initial Conditions
		Param("population",is_valid = is_count),
		Param("I_0",1,"I",func = lambda v,p,d : p_or_q(v,d["population"]),is_valid = is_count),
		Param("E_0",0,"mu E",func = lambda v,p,d : v*d['I_0'] if p == "mu" else p_or_q(v,d['population']),is_valid = is_count),
		Param("R_0",0,"R",func = lambda v,p,d : p_or_q(v,d["population"]),is_valid = is_count),
		Param("H_0",0,"H",func = lambda v,p,d : p_or_q(v,d["population"]),is_valid = is_count),
		Param("D_0",0,"D",func = lambda v,p,d : p_or_q(v,d["population"]),is_valid = is_count),
		Param("S_0",0,"population",func = lambda v,p,d : v-d["E_0"]-d["I_0"]-d["R_0"]-d["H_0"]-d["D_0"],is_valid = is_count), #useful tu calc other params
	#vaccs
		Param("Sv_0",0,"V_0 Sv",is_valid = is_count,
			func = lambda v,p,d : (lambda pop,s0: p_or_q(v,s0,s0/pop) if p == 'V_0' else p_or_q(v,s0,1))(d['population']-d['H_0']-d['D_0'],d['S_0'],),
		),
		Param("Ev_0",0,"V_0 Ev",is_valid = is_count,
			func = lambda v,p,d : (lambda pop,e0: p_or_q(v,e0,e0/pop) if p == 'V_0' else p_or_q(v,e0,1))(d['population']-d['H_0']-d['D_0'],d['E_0'])
		),
		Param("Iv_0",0,"V_0 Iv",is_valid = is_count,
			func = lambda v,p,d : (lambda pop,i0: p_or_q(v,i0,i0/pop) if p == 'V_0' else p_or_q(v,i0,1))(d['population']-d['H_0']-d['D_0'],d['I_0'])
		),
		Param("Rv_0",0,"V_0 Rv",is_valid = is_count,
			func = lambda v,p,d : (lambda pop,r0: p_or_q(v,r0,r0/pop) if p == 'V_0' else p_or_q(v,r0,1))(d['population']-d['H_0']-d['D_0'],d['R_0'])
		),
#		Param("Hv_0",0,"V_0 Hv",func = lambda v,p,d : (lambda pop,h0: p_or_q(v,h0,h0/pop) if p == 'V_0' else p_or_q(v,h0,1))(d['population'],d['H_0'])),
#Flux
		Param("S_f",0,"Flux S_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,s0: p_or_q(v,s0,s0/pop) if p == 'Flux' else p_or_q(v,s0,1))(d['population']-d['H_0']-d['D_0'],d['S_0'])
		),
		Param("E_f",0,"Flux E_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,e0: p_or_q(v,e0,e0/pop) if p == 'Flux' else p_or_q(v,e0,1))(d['population']-d['H_0']-d['D_0'],d['E_0'])
		),
		Param("I_f",0,"Flux I_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,i0: p_or_q(v,i0,i0/pop) if p == 'Flux' else p_or_q(v,i0,1))(d['population']-d['H_0']-d['D_0'],d['I_0'])
		),
		Param("R_f",0,"Flux R_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,r0: p_or_q(v,r0,r0/pop) if p == 'Flux' else p_or_q(v,r0,1))(d['population']-d['H_0']-d['D_0'],d['R_0'])
		),
	#vaccs
		Param("Sv_f",0,"V_f S_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,s0: p_or_q(v,s0,s0/pop) if p == 'Flux' else p_or_q(v,s0,1))(d['population']-d['H_0']-d['D_0'],d['Sv_0'])
		),
		Param("Ev_f",0,"V_f S_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,e0: p_or_q(v,e0,e0/pop) if p == 'Flux' else p_or_q(v,e0,1))(d['population']-d['H_0']-d['D_0'],d['Ev_0'])
		),
		Param("Iv_f",0,"V_f S_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,i0: p_or_q(v,i0,i0/pop) if p == 'Flux' else p_or_q(v,i0,1))(d['population']-d['H_0']-d['D_0'],d['Iv_0'])
		),
		Param("Rv_f",0,"V_f S_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,r0: p_or_q(v,r0,r0/pop) if p == 'Flux' else p_or_q(v,r0,1))(d['population']-d['H_0']-d['D_0'],d['Rv_0'])
		),
#Hosp
		Param("H_cap",0,func = lambda v,p,d : p_or_q(v,d["population"]),is_valid = is_count),
		
#Rates
		Param("alpha",1.0,dyn=True,is_valid = is_rate),
		Param("beta",dyn=True,is_valid = is_rate),
		Param("beta_v",0.0,dyn=True,is_valid = is_rate),
		Param("vac_eff",1.0,dyn=True,is_valid = is_frac),
		Param("vac_d",0,dyn=True,is_valid = is_count),
		Param("pE_Im",is_valid = is_frac),
		Param("tE_Im",is_valid = is_count),
		Param("pE_Icr",is_valid = is_frac),
		Param("tE_Icr",is_valid = is_count),
		Param("tEv_Iv",is_valid = is_count),# Infectious (asymptomatic + mild + severe)
		Param("tIm_R",is_valid = is_count),
		Param("tIcr_H",is_valid = is_count),	# Infectious (critical)
		Param("pIv_R",is_valid = is_frac),	# Infectious (vaccinated)
		Param("tIv_R",is_valid = is_count),
		Param("pIv_H",is_valid = is_frac),
		Param("tIv_H",is_valid = is_count),
		Param("pH_R",is_valid = is_frac),		# Hospitalized (IMV)
		Param("tH_R",is_valid = is_count),
		Param("pH_D",is_valid = is_frac),
		Param("tH_D",is_valid = is_count),
		Param("pR_S",0.0,is_valid = is_frac),
		Param("tR_S",90,is_valid = is_count)
	]

	def __init__(self,config,**kwargs):
		#super().__init__(config,**kwargs)
		self.cfg = toml.load(config)
		self.model = self.cfg["model"]
		self.initialized = False
		
	
	def init_kargs(self):
		self.kappa = ""
		self.ka_params = {}
		try:
			params = dict({"t_init":0.0,
				"I_det":False,"pI_det":1.0,
				"E_d":False},
				**self.cfg['initialconditions']
			)
			params.update(self.cfg['parameters']['static'])
			
			params.update(self.cfg['parameters'].get('dynamic',dict()))
		except KeyError as e:
			raise("Parameter section missing:\n\t"+str(e))
		
		for opt in self.PARAMETERS:
			#try to call value function of the param
			name = None
			if opt.deps == None:			#the param has the same name for toml and kappa
				name = opt.name
			elif isinstance(opt.deps,str) :	#the param has another name(s) in toml
				ignore = []
				for dep in opt.deps.split(" "):	#iter param names
					if dep in params:
						if name == None:		#first param-name defined is used
							name = dep
						else:	#if a param-name was choosen, warn about others
							ignore.append(dep)
					if ignore != []:
						print("Warning: parameters ",ignore," will be ignored cause ",name," was defined.")
				if name == None:
					if opt.dflt == None:
						raise Exception("At least one of ["+opt.deps+"] parameters is missing.")
					name = dep	
			else:
				raise Exception("Not valid dependency list (space-separated string)")
			
			if callable(opt.func):		#define value-function
				val_func = opt.func
			else:
				val_func = lambda val,name,params : val	#identity function
			
			try:
				value = params[name]
				del params[name]
				func = json.loads(value)
				if not opt.dyn:
					raise Exception("Parameter "+opt.name+" do not allow dynamic values.")
				try:
					#A dynamic value was given for the parameter
					#print("try func")
					if func['function'].upper() == 'EVENTS':
						#print("inside if")
						for (t0,tf),v in zip(func['days'],func['values']):
							val = val_func(float(v),name,self.ka_params)
							#print(t0,' ',tf,' ',v,' ',val)
							if float(t0) == 0.0:
								value = float(val)
							self.kappa = self.kappa + '%mod: [T] = '+str(t0)+" do $UPDATE '"+opt.name+"' "+str(val)+"\n"
					else:
						raise Exception("not a function: "+func['function'])
				except Exception as e:
					raise Exception("Malformed dynamic function for parameter "+opt.name+":\n\t"+e)
					
			except KeyError as e:	#the param is not defined in toml
				if opt.dflt == None:
					raise Exception("Parameter " + name + " is mandatory.")
				value = opt.dflt	#but its not mandatory (has default value)
			except (json.JSONDecodeError,TypeError) as e:		#not a dynamic value
				value = val_func(value,name,self.ka_params)
			
			try:
				if not opt.is_valid(value):
					raise Exception("Parameter "+name+" does not have a valid value for the model ("+opt.name+" = "+str(value)+")")
			except:
				raise Exception("Parameter "+name+" does not have a valid value for the model ("+opt.name+" = "+str(value)+")")
			self.ka_params[opt.name] = value
		
		self.initialized  = True
		self.eval_params = params
	
	
	def init_args(self):
		#model = toml.load(toml_file)
		self.kappa = ""
		#default-values/
		params = dict({"t_init":0.0,
			"I_det":False,"I":0,"pI_det":1.0,
			"E_d":False,"E":False,"mu":0.0,
			"R":0},
			**self.cfg['initialconditions']
		)
		try:
			params.update(self.cfg['parameters']['static'])
			
			for name,value in self.cfg['parameters'].get('dynamic',dict()).items() :
				try:
					func = json.loads(value)
				except:
					params[name] = value
					continue
				try:
					params[name] = "0.00"
					if func['function'].upper() == 'EVENTS':
						for (t0,tf),v in zip(func['days'],func['values']):
							self.kappa = self.kappa + '%mod: [T] > '+str(t0)+" do $UPDATE '"+name+"' "+str(v)+"\n"
				except e:
					print("Malformed function:\n\t",e)

			self.eval_params = dict(params)
			
			I = params['I_det']
			if(not I):    #0 or False
				I = params['I']*params['pI_det']
			E = params['E_d']
			if(not E):
				E = params['E']
				if(E != False):
					E = E * params['pI_det']#todo -???
				else:
					E = params['mu']*I
			self.eval_params['E_0'] = E  #E
			self.eval_params['I_0'] = I  #I
			self.eval_params['R_0'] = self.eval_params['R']
			del self.eval_params['I_det']
			del self.eval_params['E_d']
			del self.eval_params['pI_det']
			del self.eval_params['E']		
			del self.eval_params['I']
			del self.eval_params['R']

		
		except KeyError as e:
			print("Some SEIR parameter is missing:\n\t",e)
			
		self.initialized = True

	def kappa_sim(self,verbose = None,kappa="",runs = 0,time = 0.):
		if(not self.initialized):
			self.init_kargs()
		self.kappa = self.kappa + kappa
			
		if(self.kappa != ""):
			hd,path = tempfile.mkstemp()
			f = open(hd,"w")
			f.write(self.kappa)
			path = " " + path
			f.close()
		else:
			path = ""
		
		if time == 0.:
			time = self.eval_params['t_end'] - self.eval_params['t_init']
		del self.eval_params['t_init']
		del self.eval_params['t_end']

		out_folder = dirname+"/runs"
		model = os.path.relpath(dirname+"/kappa/SEIR-base.xka")
		if(self.model['name'] == "SEIR-byage"):
			model = os.path.relpath(dirname+"/kappa/SEIR-by_age.xka")
		elif(self.model['name'] == "SEIR-byage-vac1"):
			model = os.path.relpath(dirname+"/kappa/SEIR-by_age-vac1.xka")
		elif(self.model['name'] == "SEIR-byage-vac2"):
			model = os.path.relpath(dirname+"/kappa/SEIR-by_age-vac2.xka")
		elif(self.model['name'] == "SEIRHVD"):
			model = os.path.relpath(dirname+"/kappa/SEIRHVD.xka")

		if runs == 0:
			runs  = self.cfg['parameters'].get('RBM',dict()).get('runs',1)
		run_params = {
			"-i"	:	model+path,
			"-r"	:	str(runs),
			"-t"	:	str(time),
			"-p"	:	str(int(time)),
			"--verbose"	:	str(self.cfg['parameters'].get('RBM',dict()).get('verbose_lvl',1))#,
			#"--params"	:	" ".join([str(x) for x in eval_params])
		}
		if(type(verbose) == int):
			run_params["--verbose"] = str(verbose)
			
		if(verbose):
			print(run_params)
			print(self.kappa)
			KaXim.noisy_func()
		try:			
			self.result = KaXim.run(run_params,self.ka_params)
		except Exception as e:
			print(e)
			raise e
		try:
			os.unlink(path)
		except:
			print("Temp-file was already closed.")
		
		self.avg = self.result.getAvgTrajectory()[0].asDataFrame()
		
		self.S = self.avg['Susceptible'].values
		self.E = self.avg['Exposed'].values
		self.E_d = self.avg['Daily Exposed'].values
		self.I = self.avg['Infected'].values
		self.I_d = self.avg['Daily Infected'].values
		self.R = self.avg['Removed'].values
		self.R_d = self.avg['Daily Removed'].values
		#self.Flux = self.avg[].values

		self.E_ac = np.cumsum(self.E_d)
		#self.I_ac = np.cumsum(self.I_d) + self.I_ac
		self.R_ac = np.cumsum(self.R_d)

		#self.I_det = self.I*self.pI_det
		#self.I_d_det = self.I_d*self.pI_det
		#self.I_ac_det = self.I_ac*self.pI_det
	
		return self.result
		
	
	def solve(self,t0=0,T=None,h=0.01):
		return self.kappa_sim(False)


