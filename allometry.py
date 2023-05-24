# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 10:44:18 2022

@author: alauren
"""
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

class Allometry():
    def __init__(self):
        pass
    
    def get_motti(self, ifile, return_spe = False):
    #---read the Motti-simulation to be used as a basis for the Susi-simulation

        cnames=['yr', 'age', 'N', 'BA', 'Hg', 'Dg', 'hdom', 'vol', 'logs', 'pulp', 'loss', 'yield','mortality', 
                'stem', 'stemloss', 'branch_living', 'branch_dead', 'leaves', 'stump', 'roots_coarse', 'roots_fine']
        df = pd.read_excel(ifile, sheet_name=0, usecols=range(22), skiprows=1, header=None )
        df = df.drop([0], axis=1)
        df.columns=cnames
        cname=['idSpe']
        df2 = pd.read_excel(ifile, sheet_name=1, usecols=[4], skiprows=1, header=None )
        df2.columns = cname
        
        #---- find thinnings and add a small time to lines with the age to enable interpolation---------
        df = df.loc[df['age']!=0]    
        
        steps = np.array(np.diff(df['age']), dtype=float)
        idx = np.ravel(np.argwhere(steps<1.))+1
        df['age'][idx]=df['age'][idx]+5./365.
        
        if return_spe: 
            return df, df2['idSpe'][0]
        else:
            return df
    
    def motti_development(self, ifile, sfc):
        """
        Input:
            Motti-input file name including the folder path
        Out:
            ALL UNITS converted to /tree, except for number of stems, which is /ha 
            interpolation functions: 
                age in annual [yrs]  
                age [yrs] to variables: 
                    ageToHdom, [m] 
                    ageToBa, [m2/tree]
                    ageToVol, [m3/tree]
                    ageToYield, [m3/tree]
                    ageToBm [kg dry mass / tree]
                biomass [kg dry mass / tree] to variables:
                    bmToLeafMass, [kg/tree] 
                    bmToLAI, [m2/m2/tree]
                    bmToHdom, [m]
                    bmToYi, [m3/tree]
                    bmToBa, [m2/tree]
                    bmToLitter, [kg/ha/tree]
                    bmToStems [number/tree]
                volume or yield to variables:
                    yiToVol [m3]
                    yiToBm, [kg dry mass/tree]
                    volToLogs, [m3/tree]
                    volToPulp, [m3/tree]
                    sp    species
            Biomass models in Motti (Repola 2008, 2009) have been develped for mineral soils. In peatlands the 
            leaf mass is 35.5% lower. This bias is corrected in construction of the interpolation function
        Modifications needed:
            
            create new litter scheme: see Dec 21 esom model development
            biomass to nutrient interpolation functions: N, P, K
            nutrient to biomass interpolation functions
            litter: woody, nonwoody, locate to interpolation function 
            locate interpolation functions to dictionaty 
        """
        cnames=['yr', 'age', 'N', 'BA', 'Hg', 'Dg', 'hdom', 'vol', 'logs', 'pulp', 'loss', 'yield','mortality', 
                'stem', 'stemloss', 'branch_living', 'branch_dead', 'leaves', 'stump', 'roots_coarse', 'roots_fine']
        species_codes ={1:'Pine', 2:'Spruce', 3:'Birch'}
        df, sp = self.get_motti(ifile, return_spe=True)
        sp = sp if sp < 4 else 3
        spe = species_codes[sp]
        #leaf_scale ={1: 1.1, 2: 1.2, 3: 1.355, 4:1.4, 5: 1.5, 6: 1.6 }    # scales the leaf mass down from mineral soil, key is the site fertility class
        leaf_scale ={1: 1.0, 2: 1.0, 3: 1.355, 4:1.4, 5: 1.45, 6: 1.5 }    # scales the leaf mass down from mineral soil, key is the site fertility class

    #----modify the data frame, include age = 0-------------------------
        row = np.zeros((np.shape(df)[1]), dtype=float)
        dfrow = pd.DataFrame([row])
        dfrow.columns = cnames
        df = pd.concat([dfrow, df], axis=0)
        nrows = df.shape[0]
        df['new_index']= range(nrows)
        df = df.set_index('new_index')

        df.at[0, 'N'] = df.at[1, 'N']                                             # modify stem number at 0 yrs 
        df[['stem', 'branch_living', 'branch_dead', 'leaves', 
            'stump', 'roots_coarse', 'roots_fine']] = df[['stem','branch_living', 
                                                          'branch_dead', 'leaves', 
                                                          'stump', 'roots_coarse', 
                                                          'roots_fine']] *1000.     # unit conversion for all biomass components tn/ha -> kg/ha
        
        # Convert all to variables to values/stem
                                                    
        df[['BA',  'vol', 'logs', 'pulp', 'loss', 'yield','mortality',
            'stem', 'branch_living', 'branch_dead', 'leaves',
            'stump', 'roots_coarse', 'roots_fine']] = df[['BA',  'vol', 'logs', 'pulp', 'loss', 'yield','mortality',
                                                                   'stem', 'branch_living', 'branch_dead', 'leaves',
                                                                   'stump', 'roots_coarse', 'roots_fine']].values / df[['N',
                                                                  'N','N','N','N','N','N','N','N','N','N','N','N','N']].values
                                                                                                        
        #**********************************************************************
        # Volume and biomass variables are expressed as per mean stem: this is valid if the canopy layer diameter range is small ecough
        # attn, attn, attn, attn, attn, attn, attn, attn, attn, attn, attn,
        #**********************************************************************                                                  
        a_arr = np.arange(0, max(df['age'].values), 1.)                             # stand age array from 0 to max in Motti simulation, time step year
    
        
        #---Nutrient concentrations in tree biomass components: Palviainen & Finer 2012 Eur J For Res 131: 945-964
        #*********** Parameters *********************************************
        #---Concentrations in mg/g
        nuts = {'Pine':{
                    'Foliage':{'N': 12.5, 'P': 1.3, 'K': 4.25},
                    'Stem':{'N':1.17, 'P':0.08, 'K':0.45}},
                'Spruce':{
                    'Foliage':{'N': 12.5, 'P': 1.3, 'K': 4.25},
                    'Stem':{'N':1.12, 'P':0.09, 'K':0.64}},
                'Birch':{
                    'Foliage':{'N': 12.5, 'P': 1.3, 'K': 4.25},
                    'Stem':{'N':1.51, 'P':0.15, 'K':0.58}},
                }
                
        retrans = {'N': 0.69, 'P': 0.73, 'K':0.8}                                   # Nieminen Helmisaari 1996 Tree Phys
        rho = {'Pine': 400., 'Spruce': 380., 'Birch':480.}                          # wood density kg/m3
        sla= {'Pine': 6.8, 'Spruce': 7.25, 'Birch':14.0}                            # one-sided specific leaf area Härkönen et al. 2015 BER 20, 181-195      
    
        longevityLeaves = {'Pine':2.0, 'Spruce':4., 'Birch':1.}                     # yrs, life span of leaves and fine roots    
        longevityFineRoots ={'Pine':0.7, 'Spruce':1., 'Birch':1.}                    # Yuan & Chen 2010, turnover 1.07 times per year    
        # longevityBranch ={'Pine':22., 'Spruce':22., 'Birch':22.}                    # Pine Mäkinen 1999
        # longevityCoarseRoots ={'Pine':22., 'Spruce':22., 'Birch':22.}               # assumption as branches
        longevityBranch ={'Pine':15., 'Spruce':20., 'Birch':20.}                    # Pine Mäkinen 1999
        longevityCoarseRoots ={'Pine':15., 'Spruce':20., 'Birch':20}               # assumption as branches
    
        #********** Interpolation data ****************************************
        df['leaves'] = df['leaves'] / leaf_scale[sfc]                               # adjusting to peatland sites (Data: Hannu Hökkä 2022)
        df['leafarea'] = df['leaves'].values/10000. * sla[spe]                      # leaf area index m2 m-2
        df['stem_mass'] = df['stem'] #df['yield'] * rho[spe]                        # stem biomass
        df['stem_and_stump'] = df[['stem_mass', 'stump']].sum(axis=1)          
        df['bm'] = df[['stem_mass', 'branch_living', 'branch_dead', 'leaves',
            'stump', 'roots_coarse', 'roots_fine']].sum(axis=1)
        df['bm_noleaves'] =  df[['stem_mass','branch_living', 'branch_dead',  
            'stump', 'roots_coarse', 'roots_fine']].sum(axis=1)
        df['N_leaves'] =  df[ 'leaves'] *nuts[spe]['Foliage']['N']/1000.   
        df['Nbm_noleaves'] =  df[['stem_mass','branch_living', 'branch_dead', 'stump', 'roots_coarse']].sum(axis=1) * nuts[spe]['Stem']['N'] /1000.\
                                     + df[['roots_fine']].sum(axis=1) *nuts[spe]['Foliage']['N']/1000.       
        df['P_leaves'] =  df[ 'leaves'] *nuts[spe]['Foliage']['P']/1000.   
        df['Pbm_noleaves'] =  df[['stem_mass','branch_living', 'branch_dead', 'stump', 'roots_coarse']].sum(axis=1) * nuts[spe]['Stem']['P'] /1000.\
                                     + df[['roots_fine']].sum(axis=1) *nuts[spe]['Foliage']['P']/1000.       
        df['K_leaves'] =  df[ 'leaves'] *nuts[spe]['Foliage']['K']/1000.   
        df['Kbm_noleaves'] =  df[['stem_mass','branch_living', 'branch_dead', 'stump', 'roots_coarse']].sum(axis=1) * nuts[spe]['Stem']['K'] /1000.\
                                     + df[['roots_fine']].sum(axis=1) *nuts[spe]['Foliage']['K']/1000.       

        df['woody_logging_residues'] = df[['branch_living', 'branch_dead', 'roots_coarse', 'stump']].sum(axis=1)
        df['N_woody_logging_residues'] = df['woody_logging_residues']*nuts[spe]['Stem']['N']/1000.
        df['P_woody_logging_residues'] = df['woody_logging_residues']*nuts[spe]['Stem']['P']/1000.
        df['K_woody_logging_residues'] = df['woody_logging_residues']*nuts[spe]['Stem']['K']/1000.

        df['N_fine_roots'] = df['roots_fine']*nuts[spe]['Foliage']['N'] / 1000.
        df['P_fine_roots'] = df['roots_fine']*nuts[spe]['Foliage']['P'] / 1000.
        df['K_fine_roots'] = df['roots_fine']*nuts[spe]['Foliage']['K'] / 1000.

        df['N_leaf_demand'] = df['N_leaves']/longevityLeaves[spe]*(1.-retrans['N']) 
        df['P_leaf_demand'] = df['P_leaves']/longevityLeaves[spe]*(1.-retrans['P']) 
        df['K_leaf_demand'] = df['K_leaves']/longevityLeaves[spe]*(1.-retrans['K']) 
                
        #********** Interpolation functions ******************************************
        ageToHdom = interp1d(df['age'].values, df['hdom'].values, fill_value='extrapolate')
        ageToLAI= interp1d(df['age'].values, df['leafarea'].values,fill_value='extrapolate')        
        ageToYield = interp1d(df['age'].values, df['yield'].values,fill_value='extrapolate')
        ageToVol = interp1d(df['age'].values,df['vol'].values, fill_value='extrapolate')
        ageToBa = interp1d(df['age'].values,df['BA'].values,fill_value='extrapolate')
        ageToBm = interp1d(df['age'].values,df['bm'].values, fill_value='extrapolate')
        ageToBmNoLeaves = interp1d(df['age'].values,df['bm_noleaves'].values, fill_value='extrapolate')
        ageToStems = interp1d(df['age'].values,df['N'].values, fill_value=(df['N'].values[0], df['N'].values[-1]), bounds_error=False)
        ageToLeaves = interp1d(df['age'].values, df['leaves'].values, fill_value='extrapolate')
        ageToFineRoots = interp1d(df['age'].values, df['roots_fine'].values, fill_value='extrapolate')
        ageToBranchLiving = interp1d(df['age'].values, df['branch_living'].values, fill_value='extrapolate')
        ageToBranchDead = interp1d(df['age'].values, df['branch_dead'].values, fill_value='extrapolate')
        ageToCoarseRoots = interp1d(df['age'].values, df['roots_coarse'].values, fill_value='extrapolate')
        ageToStemStump = interp1d(df['age'].values, df['stem_and_stump'].values, fill_value='extrapolate')
        ageToNNoLeaves = interp1d(df['age'].values, df['Nbm_noleaves'].values, fill_value='extrapolate')
        ageToPNoLeaves = interp1d(df['age'].values, df['Pbm_noleaves'].values, fill_value='extrapolate')
        ageToKNoLeaves = interp1d(df['age'].values, df['Kbm_noleaves'].values, fill_value='extrapolate')

            
        volToLogs = interp1d(df['vol'].values,df['logs'].values, fill_value='extrapolate')  
        volToPulp = interp1d(df['vol'].values,df['pulp'].values,fill_value='extrapolate')
        
        yiToVol = interp1d(df['yield'].values,df['vol'].values, fill_value='extrapolate')
        yiToBm = interp1d(df['yield'].values, df['bm'].values, fill_value='extrapolate')
    
        bmToYi= interp1d(df['bm_noleaves'].values,df['yield'].values, fill_value='extrapolate')
        bmToVol= interp1d(df['bm_noleaves'].values,df['vol'].values, fill_value='extrapolate')
        
        bmToBa= interp1d(df['bm_noleaves'].values,df['BA'].values, fill_value='extrapolate')    
        bmToLeafMass = interp1d(df['bm_noleaves'].values, df['leaves'].values, fill_value='extrapolate')
        bmToLAI = interp1d(df['bm_noleaves'].values, df['leaves'].values* sla[spe]/10000.,  fill_value='extrapolate')
        bmToHdom = interp1d(df['bm_noleaves'].values,df['hdom'].values, fill_value='extrapolate')
        bmToStems = interp1d(df['bm_noleaves'].values, df['N'].values, fill_value=(df['N'].values[0], df['N'].values[-1]), bounds_error=False)
        
 
        bmToFineRoots = interp1d(df['bm_noleaves'].values, df['roots_fine'].values, fill_value='extrapolate')
        bmToNFineRoots = interp1d(df['bm_noleaves'].values, df['N_fine_roots'].values, fill_value='extrapolate')
        bmToPFineRoots = interp1d(df['bm_noleaves'].values, df['P_fine_roots'].values, fill_value='extrapolate')
        bmToKFineRoots = interp1d(df['bm_noleaves'].values, df['K_fine_roots'].values, fill_value='extrapolate')


        bmToWoodyLoggingResidues = interp1d(df['bm_noleaves'].values, df['woody_logging_residues'].values, fill_value='extrapolate')        
        bmToNWoodyLoggingResidues =  interp1d(df['bm_noleaves'].values, df['N_woody_logging_residues'].values, fill_value='extrapolate')        
        bmToPWoodyLoggingResidues =  interp1d(df['bm_noleaves'].values, df['P_woody_logging_residues'].values, fill_value='extrapolate')        
        bmToKWoodyLoggingResidues =  interp1d(df['bm_noleaves'].values, df['K_woody_logging_residues'].values, fill_value='extrapolate')        
        

        bmToNLeafDemand = interp1d(df['bm_noleaves'].values, df['N_leaf_demand'].values, fill_value='extrapolate')
        bmToPLeafDemand = interp1d(df['bm_noleaves'].values, df['P_leaf_demand'].values, fill_value='extrapolate')
        bmToKLeafDemand = interp1d(df['bm_noleaves'].values, df['K_leaf_demand'].values, fill_value='extrapolate')
        
    
        #**********************************************************************
        # We need demand functions for mass, N,P,K: 
        #              net change + fineroot_litter + woody_litter + add foliage net demand litter from other function
        # Nutrient contents in Litterfall: 
        #              fineroot litter + woody_litter + mortality_fine_root + mortality_woody + foliage litter and net demand from another function
        # Arrange: 
        #              demand functions; mass, N, P, K 
        #              litter functions: mass, N, P, K
        #***********************************************************************
        
        """ Litter functions """
        #---- Litter arrays------------------
        fineroot_litter = ageToFineRoots(a_arr)/longevityFineRoots[spe]*np.gradient(a_arr)                  # unit kg / tree / yr 
        mortality_fineroot = -np.gradient(ageToStems(a_arr))/ageToStems(a_arr)*(ageToFineRoots(a_arr))      # unitkg / tree / yr 
        mortality_leaf = -np.gradient(ageToStems(a_arr))/ageToStems(a_arr)*(ageToLeaves(a_arr))
    
        woody_litter = ageToBranchLiving(a_arr)/longevityBranch[spe]*np.gradient(a_arr) \
                    + ageToBranchDead(a_arr)/longevityBranch[spe]*np.gradient(a_arr)  \
                    + ageToCoarseRoots(a_arr)/longevityCoarseRoots[spe]*np.gradient(a_arr)                  # litterfall kg/tree in timestep
    
        mortality_woody = -np.gradient(ageToStems(a_arr))/ageToStems(a_arr) *(ageToBranchDead(a_arr) + ageToBranchLiving(a_arr) 
                                                                         + ageToStemStump(a_arr) + ageToCoarseRoots(a_arr))    
        
        dbm = np.gradient(ageToBmNoLeaves(a_arr)) + fineroot_litter + woody_litter                          # biomass change without leaves kg/ha/yr   
        
        #---- Interpolation functions -----------------
        bmToDbm = interp1d(ageToBmNoLeaves(a_arr), dbm, fill_value='extrapolate')   # from biomass to biomass change 
        bmToFinerootLitter = interp1d(ageToBmNoLeaves(a_arr), fineroot_litter, fill_value='extrapolate' )    
        bmToWoodyLitter = interp1d(ageToBmNoLeaves(a_arr), woody_litter, fill_value='extrapolate' )
        bmToMortalityFineRoot = interp1d(ageToBmNoLeaves(a_arr), mortality_fineroot, fill_value='extrapolate' )
        bmToMortalityLeaves = interp1d(ageToBmNoLeaves(a_arr), mortality_leaf, fill_value='extrapolate' )
        bmToMortalityWoody =  interp1d(ageToBmNoLeaves(a_arr), mortality_woody, fill_value='extrapolate' )
    
    
        """ Demand functions """
        #-------- arrays--------------------------
        N_fineroot_litter = (1.0-retrans['N'])*nuts[spe]['Foliage']['N'] / 1000. * fineroot_litter + mortality_fineroot*nuts[spe]['Foliage']['N'] / 1000.
        P_fineroot_litter = (1.0-retrans['P'])*nuts[spe]['Foliage']['P'] / 1000. * fineroot_litter + mortality_fineroot*nuts[spe]['Foliage']['P'] / 1000.
        K_fineroot_litter = (1.0-retrans['K'])*nuts[spe]['Foliage']['K'] / 1000. * fineroot_litter + mortality_fineroot*nuts[spe]['Foliage']['K'] / 1000.
    
        N_woody_litter = (1.0-retrans['N'])*nuts[spe]['Stem']['N'] / 1000. * woody_litter + mortality_woody*nuts[spe]['Stem']['N'] / 1000.
        P_woody_litter = (1.0-retrans['P']) *nuts[spe]['Stem']['N'] / 1000.* woody_litter + mortality_woody*nuts[spe]['Stem']['P'] / 1000.
        K_woody_litter = (1.0-retrans['K']) *nuts[spe]['Stem']['N'] / 1000.* woody_litter + mortality_woody*nuts[spe]['Stem']['K'] / 1000.
    
        N_mortality_leaves = mortality_leaf *  nuts[spe]['Foliage']['N'] / 1000.
        P_mortality_leaves = mortality_leaf *  nuts[spe]['Foliage']['P'] / 1000.
        K_mortality_leaves = mortality_leaf *  nuts[spe]['Foliage']['K'] / 1000.
    
        N_mortality_fineroot = mortality_fineroot *  nuts[spe]['Foliage']['N'] / 1000.
        P_mortality_fineroot = mortality_fineroot *  nuts[spe]['Foliage']['P'] / 1000.
        K_mortality_fineroot = mortality_fineroot *  nuts[spe]['Foliage']['K'] / 1000.
    
        N_demand = np.gradient(ageToNNoLeaves(a_arr)) + (1.0-retrans['N'])*nuts[spe]['Foliage']['N'] / 1000. * fineroot_litter\
                                                     +  (1.0-retrans['N'])*nuts[spe]['Stem']['N'] / 1000. * woody_litter
    
        P_demand = np.gradient(ageToPNoLeaves(a_arr)) + (1.0-retrans['P'])*nuts[spe]['Foliage']['P'] / 1000. * fineroot_litter\
                                                     +  (1.0-retrans['P'])*nuts[spe]['Stem']['P'] / 1000. * woody_litter
    
        K_demand = np.gradient(ageToKNoLeaves(a_arr)) + (1.0-retrans['K'])*nuts[spe]['Foliage']['K'] / 1000. * fineroot_litter\
                                                     +  (1.0-retrans['K'])*nuts[spe]['Stem']['K'] / 1000. * woody_litter
    
        #---- Interpolation functions -----------------
        bmToNdemand = interp1d(ageToBmNoLeaves(a_arr), N_demand, fill_value='extrapolate')   # from biomass to N demand No leaves here  
        bmToPdemand = interp1d(ageToBmNoLeaves(a_arr), P_demand, fill_value='extrapolate')   # from biomass to P demand 
        bmToKdemand = interp1d(ageToBmNoLeaves(a_arr), K_demand, fill_value='extrapolate')   # from biomass to K demand 
        
        bmToNFineRootLitter = interp1d(ageToBmNoLeaves(a_arr), N_fineroot_litter, fill_value='extrapolate')
        bmToPFineRootLitter = interp1d(ageToBmNoLeaves(a_arr), P_fineroot_litter, fill_value='extrapolate')
        bmToKFineRootLitter = interp1d(ageToBmNoLeaves(a_arr), K_fineroot_litter, fill_value='extrapolate')
    
        bmToNWoodyLitter = interp1d(ageToBmNoLeaves(a_arr), N_woody_litter, fill_value='extrapolate')
        bmToPWoodyLitter = interp1d(ageToBmNoLeaves(a_arr), P_woody_litter, fill_value='extrapolate')
        bmToKWoodyLitter = interp1d(ageToBmNoLeaves(a_arr), K_woody_litter, fill_value='extrapolate')
    
        bmToNMortalityLeaves =  interp1d(ageToBmNoLeaves(a_arr), N_mortality_leaves, fill_value='extrapolate')
        bmToPMortalityLeaves =  interp1d(ageToBmNoLeaves(a_arr), P_mortality_leaves, fill_value='extrapolate')
        bmToKMortalityLeaves =  interp1d(ageToBmNoLeaves(a_arr), K_mortality_leaves, fill_value='extrapolate')
    
        bmToNMortalityFineRoot =  interp1d(ageToBmNoLeaves(a_arr), N_mortality_fineroot, fill_value='extrapolate')
        bmToPMortalityFineRoot =  interp1d(ageToBmNoLeaves(a_arr), P_mortality_fineroot, fill_value='extrapolate')
        bmToKMortalityFineRoot =  interp1d(ageToBmNoLeaves(a_arr), K_mortality_fineroot, fill_value='extrapolate')
    
   
    
        allometry_f={}
        allometry_f['ageToHdom'] = ageToHdom 
        allometry_f['ageToBa'] = ageToBa
        allometry_f['ageToVol'] = ageToVol
        allometry_f['ageToYield'] = ageToYield
        allometry_f['ageToBm'] = ageToBm
        allometry_f['ageToBmNoLeaves'] = ageToBmNoLeaves
        allometry_f['ageToLeaves']= ageToLeaves
        
        allometry_f['bmToLeafMass'] = bmToLeafMass
        allometry_f['bmToLAI'] = bmToLAI
        allometry_f['bmToHdom'] = bmToHdom
        allometry_f['bmToYi'] = bmToYi
        allometry_f['bmToVol'] = bmToVol
        allometry_f['bmToBa'] = bmToBa
        allometry_f['bmToDbm'] = bmToDbm
    
        allometry_f['bmToStems'] = bmToStems
        allometry_f['yiToVol'] = yiToVol
        allometry_f['yiToBm'] = yiToBm
        allometry_f['volToLogs'] = volToLogs
        allometry_f['volToPulp'] = volToPulp
        
        allometry_f['bmToDbm'] = bmToDbm
        allometry_f['bmToFinerootLitter'] = bmToFinerootLitter
        allometry_f['bmToWoodyLitter'] = bmToWoodyLitter
        allometry_f['bmToMortalityFineRoot'] =  bmToMortalityFineRoot
        allometry_f['bmToMortalityWoody'] = bmToMortalityWoody
        allometry_f['bmToMortalityLeaves'] = bmToMortalityLeaves
    
        allometry_f['bmToNdemand'] = bmToNdemand
        allometry_f['bmToPdemand'] = bmToPdemand
        allometry_f['bmToKdemand'] = bmToKdemand
        allometry_f['bmToNFineRootLitter'] = bmToNFineRootLitter
        allometry_f['bmToPFineRootLitter'] = bmToPFineRootLitter
        allometry_f['bmToKFineRootLitter'] = bmToKFineRootLitter
    
        allometry_f['bmToNWoodyLitter'] = bmToNWoodyLitter
        allometry_f['bmToPWoodyLitter'] = bmToPWoodyLitter
        allometry_f['bmToKWoodyLitter'] = bmToKWoodyLitter
    
        allometry_f['bmToNMortalityLeaves'] = bmToNMortalityLeaves
        allometry_f['bmToPMortalityLeaves'] = bmToPMortalityLeaves
        allometry_f['bmToKMortalityLeaves'] = bmToKMortalityLeaves
    
        allometry_f['bmToNMortalityFineRoot'] = bmToNMortalityFineRoot
        allometry_f['bmToPMortalityFineRoot'] = bmToPMortalityFineRoot
        allometry_f['bmToKMortalityFineRoot'] = bmToKMortalityFineRoot
        
        allometry_f['bmToWoodyLoggingResidues'] = bmToWoodyLoggingResidues
        allometry_f['bmToNWoodyLoggingResidues'] = bmToNWoodyLoggingResidues
        allometry_f['bmToPWoodyLoggingResidues'] = bmToPWoodyLoggingResidues
        allometry_f['bmToKWoodyLoggingResidues'] = bmToKWoodyLoggingResidues
        
        allometry_f['bmToFineRoots'] = bmToFineRoots
        allometry_f['bmToNFineRoots'] = bmToNFineRoots
        allometry_f['bmToPFineRoots'] = bmToPFineRoots
        allometry_f['bmToKFineRoots'] = bmToKFineRoots
        
        allometry_f['bmToNLeafDemand'] = bmToNLeafDemand
        allometry_f['bmToPLeafDemand'] = bmToPLeafDemand
        allometry_f['bmToKLeafDemand'] = bmToKLeafDemand
        
        
        
        self.allometry_f = allometry_f
        self.sp=sp
        self.df = df
