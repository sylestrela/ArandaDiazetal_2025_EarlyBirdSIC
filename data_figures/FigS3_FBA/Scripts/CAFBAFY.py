__author__ = 'Jean Vila(jeanccvila@gmail.com)'
__version__ = '2.0.1'

''' Extension to Cobras model class to make models ready for CAFBA
	by default we deal with homogenous case as in mori et al. 
	This code is only compatabile with models using Biggs Notation. NB may work with KBASE models if flag KBASE =True set. 
	This feature was tested on April 9th 2019 and is not maintained
    NB version 2 was updated on 27/09/23 to coincide with submission of metabolic homology paper. updated transporter list to make sure it covers both notations.
'''
import cobra
import sys
import numpy as np
import os
import random

class CAFBA_Model(cobra.Model):
		
	def set_minimal_media(self,ions=[],auto_ions=True,ignore =[],KBASE=False):
		''' set media to minimal (M9) FBA by default it provides a predeterimined set of ions (auto_ions =True), you can specificy which metaoblites not to change 
			for example if you have calibrated the vmaxes'''
		if auto_ions == True:
			ions = ions + ['ca2_e', 'cl_e','cobalt2_e','cu2_e','fe2_e','fe3_e','h_e','h2o_e','k_e','mg2_e','mn2_e','mobd_e','na1_e','nh4_e','ni2_e',
							'o2_e','pi_e','so4_e','zn2_e',
							'tungs_e','sel_e','slnt_e','cbl1_e']
			if KBASE:
				ions = ions +  ['cpd29674_e0','cpd00099_e0','cpd00149_e0','cpd00058_e0','cpd10515_e0','cpd10516_e0','cpd00067_e0','cpd00001_e0',
								'cpd00205_e0','cpd00254_e0','cpd20863_e0','cpd11574_e0','cpd00971_e0','cpd19013_e0',
								'cpd00244_e0','cpd00007_e0','cpd00009_e0','cpd00048_e0','cpd00034_e0','cpd15574_e0',
								'cpd03396_e0','cpd03387_e0','cpd00635_e0']
			
		exch = self.exchanges   					
		#exch = [x for x in self.reactions  if len(x.metabolites) == 1 and next(iter(x.metabolites.values())) == -1 and next(iter(x.metabolites.keys())).compartment == 'e']
		#If compartments are not asigned to the model 
		#if len(exch) ==0:
		#	exch = [x for x in self.reactions  if len(x.metabolites) == 1 and next(iter(x.metabolites.values())) == -1 and next(iter(x.metabolites.keys())).id[-2:] == '_e' ]

		# .. media
		for l in exch:
			if next(iter(l.metabolites.keys())).id in ions:
				self.reactions.get_by_id(l.id).upper_bound = 1000
				self.reactions.get_by_id(l.id).lower_bound = -1000
			elif next(iter(l.metabolites.keys())).id not in ignore:
				self.reactions.get_by_id(l.id).upper_bound = 1000
				self.reactions.get_by_id(l.id).lower_bound = 0
			else:
				pass
				
	def set_full_media(self,default_vmax=1000,ions =[], auto_ions=True,ignore=[],KBASE=False):
		''' set model to allow uptake of all metabolits at vmax FBA. By default it provides a predetermined set of ions (auto_ions =True). 
		These models should be ready for COMETS , you can specificy which metaoblites not to change for example if you have calibrated the vmaxes'''
		if auto_ions == True:
			ions = ions + ['ca2_e', 'cl_e','cobalt2_e','cu2_e','fe2_e','fe3_e','h_e','h2o_e','k_e','mg2_e','mn2_e','mobd_e','na1_e','nh4_e','ni2_e',
							'o2_e','pi_e','so4_e','zn2_e',
							'tungs_e','sel_e','slnt_e','cbl1_e']
			if KBASE:
				ions = ions +  ['cpd29674_e0','cpd00099_e0','cpd00149_e0','cpd00058_e0','cpd10515_e0','cpd10516_e0','cpd00067_e0','cpd00001_e0',
								'cpd00205_e0','cpd00254_e0','cpd20863_e0','cpd11574_e0','cpd00971_e0','cpd19013_e0',
								'cpd00244_e0','cpd00007_e0','cpd00009_e0','cpd00048_e0','cpd00034_e0','cpd15574_e0',
								'cpd03396_e0','cpd03387_e0','cpd00635_e0','cpd00027_e0']
		exch = [x for x in self.reactions  if len(x.metabolites) == 1 and next(iter(x.metabolites.values())) == -1 and next(iter(x.metabolites.keys())).compartment == 'e']
		# .. media
		for l in exch:
			if next(iter(l.metabolites.keys())).id in ions:
				self.reactions.get_by_id(l.id).lower_bound = -1000
			elif next(iter(l.metabolites.keys())).id not in ignore:
				self.reactions.get_by_id(l.id).lower_bound = -default_vmax
			else:
				pass
				
	def calculate_cost(self,metabolite,grate):
		#calculates the cafba cost that would force the vmax to match the growth rate.
		exch = [x.id for x in self.metabolites.get_by_id(metabolite).reactions  if len(x.metabolites) == 1 
				and next(iter(x.metabolites.values())) == -1][0]
		obj = [x.id for x in self.reactions if x.objective_coefficient ==1][0]
		obj2 = 'cafba_rxn'
		no_change =False
		with self as m:
			m.set_minimal_media(ignore=[metabolite])
			m.reactions.get_by_id(obj).upper_bound =grate
			# The first optimization is necessary incase the emperically observed growth rate is infeasible given the CAFBA constraint.
			sol1 = m.optimize()
			if sol1.objective_value < grate:
				print('warning maximum growth rate for ' + metabolite + ' is less the observed')
				no_change = True
			m.reactions.get_by_id(obj).lower_bound =sol1.objective_value
			
			# The second optimization is to maximize the amount of budget left available whilst achieving maixmum growth
			m.objective = obj2
			sol2 = m.optimize(objective_sense='minimize')
			free_cafba = m.reactions.cafba_rxn.upper_bound- sol2.objectiviae_value
			if free_cafba ==0:
				print('no available budget at this growth rate for ' + metabolite)
				no_change = True	
			m.reactions.cafba_rxn.upper_bound = sol2.objective_value
				
			# The third optimization is to find the minimum uptake that delivers total usage of cafba and achieves maximum growth
			m.reactions.cafba_rxn.upper_bound = sol2.objective_value

			m.objective = exch
			sol3 = m.optimize()
			vmax = -sol3.objective_value	
			cost1 = free_cafba/vmax				
		if no_change:
			# if either criterion holds don't add a cost
			free_cafba = 0
			vmax = 1
		return(free_cafba/vmax)

	def transporters_list(self,ignore_list = ['DMSOR1e', 'DMSOR2e', 'GLCDe', 'TMAOR1e', 'TMAOR2e'],auto_ignore =True):
		'''take models and returns list of transporters it. To do this 
        it looks for reactions that move metabolites across compartments
		it ignores protons when deciding if a reaction is a transporters as well as sodium (also the cafba_met if model is cafbafied)
		ignore list is list of reacts that are not transporters, but are involved in oxidative phosphorylation / intracellular metabolism. '''	
		def is_transport(react):
			mets = {key:value for key, value in react.metabolites.items() if str(key)[:-1] !='h_' 
                and str(key)[:-2] !='cpd00067_'
                and str(key)[:-1] !='na1_'  
                and str(key)[:-2] != 'cpd00971_'
                and str(key) != 'CAFBA_met'}
			con_mets = [key.compartment for key, value in mets.items() if value < 0 ]
			prod_mets = [key.compartment for key, value in mets.items() if value > 0 ]
			#does metabolite go from cytoplasm or periplasm to extracellular
			test = ('e' in prod_mets and  ('c' in con_mets or 'p' in con_mets)) or ('C_e' in prod_mets and  ('C_c' in con_mets or 'C_p' in con_mets))
			#does metabolite go cytoplasm to periplasm
			test2 = ('p' in prod_mets  and  ('c' in con_mets)) or ('c_p' in prod_mets  and  ('c_c' in con_mets))
			#does metabolite go from extracellular compt to cytoplasm or periplasm
			test3 = ('e' in con_mets and  ('c' in prod_mets or 'p' in prod_mets)) or ('C_e' in con_mets and  ('C_c' in prod_mets or 'C_p' in prod_mets))
			#does metabolite go periplasm to cytoplasm
			test4 = 'p' in con_mets  and  ('c' in prod_mets) or ('C_p' in con_mets  and  ('C_c' in prod_mets))
			test5 = len(mets) < 2 #Sinks and Exchanges
			return test or test2 or test3 or test4 or test5
		rxn = [y.id for y in self.reactions if is_transport(y) or y.objective_coefficient ==1]
		if auto_ignore ==True:
			rxn = [x for x in rxn if x not in ignore_list]
		return rxn
		
	def irreversabilize(self,unbound_objective=True):
		'''	Function makes all metabolic reactions foward and irreversible
		Reactions with a negative lower bound and an upper bound of 0.0 are flipped
		Reaction with a negative lower bound and positive upper bound are split in 2 and a new copy is added with a new ID ('_rev')'''

		if unbound_objective:
			''' some models have the objective bound for some reason (i.e ijN764) by default we don't want this to be the case'''
			self.reactions.get_by_id([x.id for x in self.reactions if x.objective_coefficient ==1][0]).lower_bound = 0.0
			self.reactions.get_by_id([x.id for x in self.reactions if x.objective_coefficient ==1][0]).upper_bound = 1000.0
		for i in range(0, len(self.reactions)):
			if self.reactions[i].lower_bound < -999:
				self.reactions[i].lower_bound = -1000.0
			if self.reactions[i].upper_bound > 999:
				self.reactions[i].upper_bound = 1000.0
			if (self.reactions[i].lower_bound <0 and  len(self.reactions[i].metabolites)!=1 and self.reactions[i].objective_coefficient==0):  #Exchange and sink reactions
				pass
			else:
				continue
				
			# copy reversible reaction and split
			temp_rxn = self.reactions[i].copy()
			self.reactions[i].lower_bound = 0
			# change bounds of new reaction
			temp_rxn.upper_bound = -1*temp_rxn.lower_bound
			temp_rxn.lower_bound = 0
			
			#change metabolite signs of new reaction
			new_mets = temp_rxn.metabolites
			new_mets = {k: -v for (k, v) in new_mets.items()}
			temp_rxn.subtract_metabolites(temp_rxn.metabolites)
			temp_rxn.add_metabolites(new_mets)
			# add "rev" to new rxn name and add it to model
			temp_rxn.id = temp_rxn.id + '_rev'
			self.add_reactions([temp_rxn])
			
	def CAFBAfy(self,wc =0.0, wr = 0.169,wi = 0.00083, phiQ_0 =0.45,phiR_0 =0.066,transporters =[],auto_transporters =True):
		''' Implements a the CAFBA constraints every reaction in a given sector has the same cost Custom Costs can be introduced down the line. 
		By default transporters are identified automatically but you can parse them in as a list of reaction id.
		wc = c-sector cost
		wr = r-sector cost
		wi = e-sector cost
		phiR_0  = growth indepdnent cost of ribosome production
		phiQ_0 = growth indepdent cost of housekeeping proteins.
		default paramaters values are taken from the CAFBA paper.'''
		CAFBA_met = cobra.Metabolite('CAFBA_met',formula='',name='CAFBA cost fake metabolite', compartment='c')
		if auto_transporters ==True or len(transporters) >0:
			#Automatically identifies transporters by looking at metabolite compartments
			transporters = self.transporters_list() 
		for rxn in self.reactions:
			if len(rxn.metabolites) == 1 and next(iter(rxn.metabolites.values())) ==-1:
				continue
			elif rxn.id =='ATPM':
				continue #ATP maintenance is not included in the global constraint
			elif rxn.id in transporters:
				rxn.add_metabolites({CAFBA_met : -wc},combine=False) # not E sector (i.e. transport)
			elif rxn.objective_coefficient != 0.0:
				rxn.add_metabolites({CAFBA_met : -wr},combine=False)
				rxn.lower_bound = 0.0
				rxn.upper_bound = 1000.0
			else:
				rxn.add_metabolites({CAFBA_met : -wi},combine=False) # E sector

		
		#finishing by adding cafba budget
		cafba_rxn = cobra.Reaction('cafba_rxn')
		cafba_rxn.name = 'CAFBA budget reaction'
		cafba_rxn.lower_bound = 0
		cafba_rxn.upper_bound = (1-(phiQ_0 + phiR_0))
		cafba_rxn.add_metabolites({CAFBA_met : 1})		
		self.add_reactions([cafba_rxn])

	def customCost(self,custom_costs = {}):
		''' set custom costs from dictionary objects (keys are reaction ID item is the costs'''
		um = self.metabolites.CAFBA_met
		for keys in custom_costs:
			t =  self.reactions.get_by_id(keys)
			t.add_metabolites({um : custom_costs[keys]},combine=False)
				
		             
	def prune_model(self,rxn_list=[],AUTO = True, AUTO_GLCDe = True , AUTO_Univ=False):
		''' Remove unused set of reactions and unusued metabolitess. 
		Auto provides the option to remove remove periplasmic glucose dehydrogenase as in the original CAFBA paper  and 
		
		For CARVE me Universal remove a cdefined set automatically (this forces the model to behave in an Ecoli like manner(
		
		AUTO alone just removes fully bounded reactions
		
		At some point i will add capacity to remove duplicate reactions''' 
	    
		if AUTO_GLCDe:
			rxn_list = rxn_list + ['GLCDe','GLCDpp'] # as in the original CAFBA paper.
			# Glycerol dehydrogenase is also removed to force glycerol to go through glycerol kinase (Murarka 2008)
			rxn_list = rxn_list + ['rxn08608_p0'] # for KBASE
		if AUTO_Univ:
			rxn_list = rxn_list + ['OOR2', 'OOR2r', # To force TCA cycle early on 
				'ALC19','ALCD19_rev','GLYCDx', # Glyceorl dehyodrgenases
				'NADHNQR','POR', 'POR_2','CBFCpp','CBFC2pp','AGTi_rev','SPT_syn_rev'] 	# Force glyxolyate shunt off ea;ry on
		if AUTO:
			rxn_list = rxn_list + [y.id for y in self.reactions if y.lower_bound ==0.0 and y.upper_bound ==0.0]
		for i in rxn_list:
			try:
				self.reactions.get_by_id(i).delete()
			except:
				continue
				
		for j in self.metabolites:
			if len(j.reactions) ==0:
				j.remove_from_model()
				

	def calibrate_single_uptake_cost(self,csource,grate,thresh=50):
		''' Calibrate model to growth rates in a single carbon sources by tuning uptake costs.'''
		#For each carbon source calculate cost that would deliver corresponding growth rate
		# Calculate cost
		cost = self.calculate_cost(csource,grate)
		# First lets find the transporters for this metabolite
		t = self.metabolites.get_by_id(csource).reactions
		t = [y for y in t]
		t = [y.id for y in t if len(y.metabolites) != 1 and y.get_coefficient(csource) <0 and y.objective_coefficient ==0]       
		# Lets set find the transporters for this metabolite
		self.customCost(dict(zip(t,[-cost for x in t])))
		newcost = cost
		with self as m:
			counter =0
			m.set_minimal_media(ignore=[csource])
			f = m.slim_optimize()
			while f - grate > 1e-6:
				counter = counter+1
				m.set_minimal_media(ignore=[csource])
				f = m.slim_optimize()
				newcost = newcost*(f/grate)
				m.customCost(dict(zip(t,[-newcost for x in t])))
				if counter>thresh:
					break
					print('Failed to converge for ' + csource)
		self.customCost(dict(zip(t,[-newcost for x in t])))

				
		print(csource +' ' + str(cost))
		
	def calibrate_uptake_cost(self,csources,grates):
		''' Takes a list of carbon sources and the growth rate on each carbon source, and sets the costs of uptake to reproduce the observed growth rates.'''
		for i in range(0,len(csources)):			
			self.calibrate_single_uptake_cost(csources[i],grates[i])
   
	def set_we(self,cost = 0.00083):
		''' Takes a list of carbon sources and the growth rate on each carbon source, and sets the costs of uptake to reproduce the observed growth rates.'''
		transporters = self.transporters_list() 
		e_sector_reactions_g = [x for x in self.reactions if x.id not in transporters and x.id != 'ATPM' and x.objective_coefficient ==0.0 and not (len(x.metabolites) == 1 and next(iter(x.metabolites.values())) ==-1) and x.id !='cafba_rxn']
		um = self.metabolites.CAFBA_met
		for t in e_sector_reactions_g:
			t.add_metabolites({um : -cost},combine=False)			
															
	def calibrate_e_sector(self,csource = 'glc__D_e',grate =1.0,thresh=50,starting_cost = 0.00083,KBASE=False):
		'''So the cost of the E sector will depend on the size of the model.
		 As in Mori et al we set the E-sector per reaction cost, such thatThe homogeneous case: 
		 here, wi’s are uniformly set to the same value, denoted by wE, for each reaction i. 
		 wE is chosen so that the fastest growth rate achievable λmax, 
		 corresponding to wC = 0, matches the corresponding empirical value.  I use a simple convergence algorthm to find a cost that delivers the corresponding grate.'''
		#Non_Transporter_Reactions
		if KBASE:
			csource = 'cpd00027_e0'
		transporters = self.transporters_list() 
		cost = starting_cost
		self.set_minimal_media(ignore=[csource],KBASE=KBASE)
		self.slim_optimize()
		with self as m:
			e_sector_reactions = [x for x in m.reactions if x.id not in transporters 
				and x.id != 'ATPM' 
				and x.objective_coefficient ==0.0 
				and not (len(x.metabolites) == 1
				and next(iter(x.metabolites.values())) ==-1) 
				and x.id !='cafba_rxn']
			um = m.metabolites.CAFBA_met
			counter =0
			m.set_minimal_media(ignore=[csource],KBASE=KBASE)
			f = m.slim_optimize()
			while np.abs(f - grate) > 1e-6:
				counter = counter+1
				m.set_minimal_media(ignore=[csource],KBASE=KBASE)
				cost = cost*(f/grate)
				for t in e_sector_reactions:
					t.add_metabolites({um : -cost},combine=False)
					if counter>thresh:
						print('Failed to converge for ' + csource)
						break
					else:
						f = m.slim_optimize()
			print('w_e = ' + str(cost))
		e_sector_reactions_g = [x for x in self.reactions if x.id not in transporters and x.id != 'ATPM' and x.objective_coefficient ==0.0 and not (len(x.metabolites) == 1 and next(iter(x.metabolites.values())) ==-1) and x.id !='cafba_rxn']
		um = self.metabolites.CAFBA_met
		for t in e_sector_reactions_g:
			t.add_metabolites({um : -cost},combine=False)			
																																										
	def save(self,filename):
		cobra.io.write_sbml_model(self,filename)
		
	def reassign_compartments(self):
		''' if the metabolites in the model are not allocated to compartment, do this manually.'''
		if set(self.compartments.keys()) == set(['e','p','c']):
			next
		elif set(self.compartments.keys()) == set(['c0', 'e0', 'p0']):
			print('compartments are kbase messed up, manually reassigning compartments')
			for i in self.metabolites:
				i.compartment = i.compartment[0] 
		else:
			print('compartments are not default, manually reassigning compartments')
			for i in self.metabolites:
				i.compartment = i.compartment[-1] 
				
	def set_ATPM(self,default_ATPM = 3.15):
		''' set the ATPM lower bound.'''
		self.reactions.ATPM.lower_bound=default_ATPM
	
	def quick_cafbafy(self,wr=0.169,phi_q=0.45,grate =1,i=False,KBASE=False):
		self.reassign_compartments()
		self.irreversabilize()
		self.prune_model()
		self.CAFBAfy(wr = wr,phiQ_0 =phi_q)
		self.set_full_media()
		self.calibrate_e_sector(grate=grate,KBASE=KBASE)
		
def main(xml_filename,output_file=None,AUTO_Univ=False,AUTO_GLCDe=True,KBASE=False):		
	x = CAFBA_Model(cobra.io.read_sbml_model(xml_filename))
	x.reassign_compartments()
	x.irreversabilize()
	x.prune_model(AUTO_Univ=AUTO_Univ, AUTO_GLCDe=AUTO_GLCDe,KBASE=KBASE)
	x.CAFBAfy()
	x.set_full_media()
	x.calibrate_e_sector(KBASE=KBASE)
	if output_file==None:
		filename, file_extension = os.path.splitext(xml_filename)
		x.save(filename + '_cafba' + file_extension)
	else: 
		filename, file_extension = os.path.splitext(output_file)
		print(filename + '_cafba' + file_extension)
		x.save(filename + '_cafba' + file_extension)

if __name__ == "__main__":
	if len(sys.argv) == 3:
		main(sys.argv[1],sys.argv[2])
	else:
		main(sys.argv[1])

