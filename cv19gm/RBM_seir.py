
import toml
import os
import sys
dirname = os.path.dirname(__file__)
if(dirname+"/../lib" not in sys.path):
	sys.path.insert(1,dirname+"/../lib")
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
	def __init__(self,name,dflt = None,deps = None,dyn = False,func = None,is_valid = None,desc = ""):
		self.name = name			#Name of the parameter
		self.dflt = dflt			#default value (None if mandatory)
		self.deps = deps			#dependencies from toml params (sorted by prevalence, space-separated) 
		self.dyn = dyn				#True if dynamic parameter
		self.func = func			#value = func(v = toml[param_name],p = param_name, d = values_dict)
		self.is_valid = is_valid	#if not checks(value) then bad-parameter-definition. 
		self.desc = desc			#A brief description of the parameter 
		

class RBM_SEIR():

	STATES = ["S","E","I","R"]
	STATES_NAMES = ["Susceptible","Exposed","Infected","Removed"]
	
	PARAMETERS = [
#Initial Conditions
		Param("population",is_valid = is_count,
			desc = "Initial population"),
		Param("I_0",1,"I",func = lambda v,p,d : p_or_q(v,d["population"]),is_valid = is_count,
			desc = "Initial infected, [I] as % of population if a [0,1] value is given"),
		Param("E_0",0,"mu E",func = lambda v,p,d : v*d['I_0'] if p == "mu" else p_or_q(v,d['population']),is_valid = is_count,
			desc = "Initial exposed, [mu] as % of I_0, [E]  as % of population if a [0,1] value is given"),
		Param("R_0",0,"R",func = lambda v,p,d : p_or_q(v,d["population"]),is_valid = is_count,
			desc = "Initial recovered, [R] as % of population if a [0,1] value is given"),
		Param("H_0",0,"H",func = lambda v,p,d : p_or_q(v,d["population"]),is_valid = is_count,
			desc = "Initial hostpitalized, [H] as % of population if a [0,1] value is given" ),
		Param("D_0",0,"D",func = lambda v,p,d : p_or_q(v,d["population"]),is_valid = is_count,
			desc = "Initial dead, [D] as % of population if a [0,1] value is given" ),
		Param("S_0",0,"population",func = lambda v,p,d : v-d["E_0"]-d["I_0"]-d["R_0"]-d["H_0"]-d["D_0"],is_valid = is_count), #useful to calc other params
	#vaccs
		Param("Sv_0",0,"V_0 Sv",is_valid = is_count,
			func = lambda v,p,d : (lambda pop,s0: p_or_q(v,s0,s0/pop) if p == 'V_0' else p_or_q(v,s0))(d['population']-d['H_0']-d['D_0'],d['S_0']),
			desc = "Initial vaccinated suscp., [V_0] as a % of S_0 if [0,1] value else as S_0/population % of V_0, [Sv] as % of S_0 if [0,1] value"
		),
		Param("Ev_0",0,"V_0 Ev",is_valid = is_count,
			func = lambda v,p,d : (lambda pop,e0: p_or_q(v,e0,e0/pop) if p == 'V_0' else p_or_q(v,e0))(d['population']-d['H_0']-d['D_0'],d['E_0']),
			desc = "Initial vaccinated Expsd., [V_0] as a % of E_0 if [0,1] value else as E_0/population % of V_0, [Ev] as % of E_0 if [0,1] value"
		),
		Param("Iv_0",0,"V_0 Iv",is_valid = is_count,
			func = lambda v,p,d : (lambda pop,i0: p_or_q(v,i0,i0/pop) if p == 'V_0' else p_or_q(v,i0))(d['population']-d['H_0']-d['D_0'],d['I_0']),
			desc = "Initial vaccinated infct., [V_0] as a % of I_0 if [0,1] value else as I_0/population % of V_0, [Iv] as % of I_0 if [0,1] value"
		),
		Param("Rv_0",0,"V_0 Rv",is_valid = is_count,
			func = lambda v,p,d : (lambda pop,r0: p_or_q(v,r0,r0/pop) if p == 'V_0' else p_or_q(v,r0))(d['population']-d['H_0']-d['D_0'],d['R_0']),
			desc = "Initial vaccinated rcvrd., [V_0] as a % of R_0 if [0,1] value else as R_0/population % of V_0, [Rv] as % of R_0 if [0,1] value"
		),
#		Param("Hv_0",0,"V_0 Hv",func = lambda v,p,d : (lambda pop,h0: p_or_q(v,h0,h0/pop) if p == 'V_0' else p_or_q(v,h0,1))(d['population'],d['H_0'])),
#Flux
		Param("S_f",0,"Flux S_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,s0: p_or_q(v,s0,s0/pop) if p == 'Flux' else p_or_q(v,s0))(d['population']-d['H_0']-d['D_0'],d['S_0']),
			desc = "Suscp. flux, [Flux] as a % of S_0 if [0,1] value else as S_0/population % of Flux, [S_f] as % of S_0 if [0,1] value"
		),
		Param("E_f",0,"Flux E_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,e0: p_or_q(v,e0,e0/pop) if p == 'Flux' else p_or_q(v,e0))(d['population']-d['H_0']-d['D_0'],d['E_0']),
			desc = "Expsd. flux, [Flux] as a % of E_0 if [0,1] value else as E_0/population % of Flux, [E_f] as % of E_0 if [0,1] value"
		),
		Param("I_f",0,"Flux I_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,i0: p_or_q(v,i0,i0/pop) if p == 'Flux' else p_or_q(v,i0))(d['population']-d['H_0']-d['D_0'],d['I_0']),
			desc = "Infct. flux, [Flux] as a % of I_0 if [0,1] value else as I_0/population % of Flux, [I_f] as % of I_f if [0,1] value"
		),
		Param("R_f",0,"Flux R_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,r0: p_or_q(v,r0,r0/pop) if p == 'Flux' else p_or_q(v,r0))(d['population']-d['H_0']-d['D_0'],d['R_0']),
			desc = "Rcvrd. flux, [Flux] as a % of R_0 if [0,1] value else as R_0/population % of Flux, [R_f] as % of R_f if [0,1] value"
		),
	#vaccs
		Param("Sv_f",0,"V_f Sv_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,s0: p_or_q(v,s0,s0/pop) if p == 'Flux' else p_or_q(v,s0))(d['population']-d['H_0']-d['D_0'],d['Sv_0']),
			desc = "Vacc. Suscp. flux, [V_f] as a % of Sv_0 if [0,1] value else as Sv_0/population % of Flux, [Sv_f] as % of Sv_0 if [0,1] value"
		),
		Param("Ev_f",0,"V_f Ev_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,e0: p_or_q(v,e0,e0/pop) if p == 'Flux' else p_or_q(v,e0))(d['population']-d['H_0']-d['D_0'],d['Ev_0']),
			desc = "Vacc. Suscp. flux, [V_f] as a % of Ev_0 if [0,1] value else as Ev_0/population % of Flux, [Ev_f] as % of Ev_0 if [0,1] value"
		),
		Param("Iv_f",0,"V_f Iv_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,i0: p_or_q(v,i0,i0/pop) if p == 'Flux' else p_or_q(v,i0))(d['population']-d['H_0']-d['D_0'],d['Iv_0']),
			desc = "Vacc. Suscp. flux, [V_f] as a % of Iv_0 if [0,1] value else as Iv_0/population % of Flux, [Iv_f] as % of Iv_0 if [0,1] value"
		),
		Param("Rv_f",0,"V_f Rv_f",dyn=True,is_valid = is_count,
			func = lambda v,p,d : (lambda pop,r0: p_or_q(v,r0,r0/pop) if p == 'Flux' else p_or_q(v,r0))(d['population']-d['H_0']-d['D_0'],d['Rv_0']),
			desc = "Vacc. Suscp. flux, [V_f] as a % of Rv_0 if [0,1] value else as Rv_0/population % of Flux, [Rv_f] as % of Rv_0 if [0,1] value"
		),
#Hosp
		Param("H_cap",0,func = lambda v,p,d : p_or_q(v,d["population"]),is_valid = is_count,
			desc = "Hospital bed capability, as a % of population if [0,1] value"
		),
#Rates
		Param("alpha",1.0,dyn=True,is_valid = is_rate,desc = "Movility rate"),
		Param("beta",dyn=True,is_valid = is_rate,desc = "Contagious rate"),
		Param("beta_v",0.0,dyn=True,is_valid = is_rate,desc = "Vaccinated contagious rate"),
		Param("vac_eff",1.0,dyn=True,is_valid = is_frac,desc = "Vaccine effectivity rate"),
		Param("vac_d",0,dyn=True,is_valid = is_count,desc = "Daily vaccination"),
		Param("pE_Im",is_valid = is_frac,desc = "% to develop mild symptoms"),
		Param("tE_Im",is_valid = is_count,desc = "Avg. days for mild symptoms incubation"),
		Param("pE_Icr",is_valid = is_frac,desc = "% to develop critical symptoms"),
		Param("tE_Icr",is_valid = is_count,desc = "Avg. days for critical symptoms incubation"),
		Param("tEv_Iv",is_valid = is_count,desc = "Avg. days for vaccinated incubation"),# Infectious (asymptomatic + mild + severe)
		Param("tIm_R",is_valid = is_count,desc = "Avg. days for mild symp. recovery"),
		Param("tIcr_H",is_valid = is_count,desc = "Avg. days of crit. symp. hospitalization"),	# Infectious (critical)
		Param("pIv_R",is_valid = is_frac,desc = "% of vaccinated recovery"),	# Infectious (vaccinated)
		Param("tIv_R",is_valid = is_count,desc = "Avg. days for vaccinated recovery"),
		Param("pIv_H",is_valid = is_frac,desc = "% of vaccinated hospitalization"),
#		Param("tIv_H",is_valid = is_count,desc = "Avg. days for vaccinated hospitalization"),
		Param("pH_R",is_valid = is_frac,desc = "% of hosp. recovery"),		# Hospitalized (IMV)
		Param("tH_R",is_valid = is_count,desc = "Avg. days for hosp. recovery"),
		Param("pH_D",is_valid = is_frac,desc = "% of hosp. deaths"),
		Param("tH_D",is_valid = is_count,desc = "Avg. days for hosp. death"),
		Param("pR_S",0.0,is_valid = is_frac,desc = "% of loss immunity"),
		Param("tR_S",90,is_valid = is_count,desc = "Avg. days for immunity loss")
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

	@classmethod
	def paramInfo(cls):
		import pandas as pd
		tab = []
		for param in cls.PARAMETERS:
			if param.desc == "":
				continue
			l = [param.name,param.deps,param.dflt,param.dyn,param.desc]
			tab.append(l)
		df = pd.DataFrame(tab,columns = ["Name","Alias","Default","Dynamic","Description"])
		import tabulate
		print(tabulate.tabulate(df,headers = 'keys'))
		return df





