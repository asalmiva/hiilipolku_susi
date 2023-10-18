# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:50:19 2023

@author: 03081194
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 10:13:04 2022

@author: alauren
modified by Leena Stenberg and Aura Salmivaara
"""
from netCDF4 import Dataset 
from datetime import datetime
import numpy as np
import datetime
import pandas as pd
import os
import shutil
from glob import glob
from numpy.ma import masked_array
import pandas as pd

#yksittaisen ncf tarkastelu
#scenarios=['Kuonanjoki_BAU_A', 'Kuonanjoki_BAU_B', 'Kuonanjoki_BIO_A', 'Kuonanjoki_BIO_B', 'Kuonanjoki_HII_A', 'Kuonanjoki_HII_B']
#scenarios=['Sorvasranta_BAU_A', 'Sorvasranta_BAU_B', 'Sorvasranta_BIO_A', 'Sorvasranta_BIO_B', 'Sorvasranta_HII_A', 'Sorvasranta_HII_B']
scenarios=['Halvanjoki_BAU_A', 'Halvanjoki_BAU_B', 'Halvanjoki_BIO_A', 'Halvanjoki_BIO_B', 'Halvanjoki_HII_A', 'Halvanjoki_HII_B']

#scenario=scenarios[0]#
#mottifol=r'/scratch/project_2002470/HIILIPOLKU_data/'
resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'
#scenario=scenarios[0]
missing_stands_list=[]
newrunslist=[]
def check_resfiles(scenarios):
    for scenario in scenarios:                    
        mottifiles=glob(mottifol+scenario[:-6]+'/'+scenario+'/*.xls')
        files=glob(resfol+scenario+'/*.nc')
        resstands=[i[(len(resfol)+len(scenario)+1):-6] for i in files]
        resstands = [*set(resstands)]
        resstands.sort()
        stands=[i[(len(mottifol)+len(scenario[:-6])+1+len(scenario)+1):-7] for i in mottifiles]
        stands = [*set(stands)]
        stands.sort()
        if len(resstands)==len(stands):
            print(scenario, " ok")
        else:
            missing_stands=list(set(stands) - set(resstands))   #list(set(common_ms).intersection(dtw_list))
            missing_stands_list.append(missing_stands)
            missing_stands_list.sort()
            print(scenario, len(stands)-len(resstands),"missing stands: ", missing_stands) 
            flines=[]
            with open(resfol+scenario+'/'+scenario+'.txt','rt') as myfile:
                for myline in myfile:
                    flines.append(myline.rstrip('\n'))
            nmlist=["i="+str(i) for i in range(len(stands))]
            nmfound=[]
            for line in flines:
                for i in nmlist:
                    if line.find(i)!=-1:
                        nmfound.append(i)
            nm_missing=list(set(nmlist) - set(nmfound))
            nm_missing.sort()
            print(nm_missing)
            np.unique([nms[3:] for nms in nm_missing])
            newruns=[int(nm) for nm in np.unique([nms[3:] for nms in nm_missing])]
            print(','.join(str(i) for i in newruns)) #this to run_hiilipolku_csc_kuonan_baua_repair.sh file and sbatch ***.sh file 
            newrunslist.append(scenario +" "+ ','.join(str(i) for i in newruns))
    return missing_stands_list, newrunslist 

def atomic_flatten(iterable, flattened):
    try:
        iter(iterable)
        if type(iterable) not in [str, bytes]:
            for item in iterable:
                atomic_flatten(item, flattened)
        else:
            flattened.append(iterable)
    except:
        flattened.append(iterable)
    return flattened

#sorvar_ms,sorvas_newruns=check_resfiles(scenarios)
#sorvas_ms=atomic_flatten(sorvar_ms,[])
#sorvas_ms.sort()
#sorvas_ms=[*set(sorvas_ms)]

#
#
#files=glob(resfol+scenario+'/*.nc')
#stands=[i[(len(resfol)+len(scenario)+1):-6] for i in files]
#stands = [*set(stands)]
#stands.sort()

#if any(np.isnan((np.unique(wt.data)))):
#    print(f, ev)


#ev='n0'    
#f='43132938'
#f=sorvs_prob[6]

#ff=resfol+scenario+"/"+f+"_"+ev+".nc"
#ncf=Dataset(ff, mode='r')
#dscen=1
#wt = ncf['strip']['dwtyr'][dscen,:, :]
#wt
#ncf['esom']['Mass']['P1'][dscen]
#indices=np.where(~ncf['esom']['Mass']['P1'][dscen].mask) #ei aina toimi!!!
#lastind=max(indices[0])
#lastind
#ncf['esom']['Mass']['P1'][dscen].data

#ncf['balance']['N']['to_water'][1]
#indices=np.where(~ncf['balance']['N']['to_water'][1].mask)
#lastind=max(indices[0])
#lastind
#ncf['balance']['N']['to_water'][1].data
#f


#if type(ncf['balance']['N']['to_water'][1].mask)==np.ndarray:
#    print("loytyy mask info, voi kayttaa indeksointia tulosten koonnissa")

#type(ncf['esom']['Mass']['P1'][dscen].mask)
#if ~ncf['esom']['Mass']['P1'][dscen].mask:
#    print("array maks is false , cannot use indeksointi tulosten koonnistta, eka rivia ei voi laske muut lasketaan")

#

#cols = np.shape(wt)[0]
#esoms =['L0L', 'L0W', 'LL', 'LW', 'FL', 'FW', 'H', 'P1', 'P2', 'P3'] 
#for sto in esoms[7:9]:
#    inipeat += ncf['esom']['Mass'][sto][dscen,0, :]/10000.
#inipeatb = 120/190*ncf['esom']['Mass'][esoms[9]][dscen,0, :]/10000. #syvyys skaalaus eli 2.5m turvekerroksesta vain 1.8m tms mukaan, tai bottomista 60cm-180cm eli 120 koko 190cm sta
#for sto in esoms[7:9]:
#    indices=np.where(~ncf['esom']['Mass'][sto][dscen].mask)
#    lastind=max(indices[0])
#    endpeat += ncf['esom']['Mass'][sto][dscen,lastind, :]/10000. #pitais ottaa vika jossa on arvoja 
#indices=np.where(~ncf['esom']['Mass'][esoms[9]][dscen].mask)
#lastind=max(indices[0])
#endpeatb = 120/190*ncf['esom']['Mass'][esoms[9]][dscen,lastind, :]/10000. #pitais ottaa vika jossa on arvoja 
#endpeat=endpeat+endpeatb
#inimor =  np.zeros(cols)
#for sto in esoms[:7]:
#    inimor += ncf['esom']['Mass'][sto][dscen,0, :]/10000.
#endmor =  np.zeros(cols)
#for sto in esoms[:7]:
#    indices=np.where(~ncf['esom']['Mass'][sto][dscen].mask)
#    lastind=max(indices[0])
#    endmor += ncf['esom']['Mass'][sto][dscen,lastind, :]/10000.    
#

#indices=np.where(~ncf['esom']['Mass']['P1'][dscen].mask) #ei aina toimi!!!
#lastind=max(indices[0])
#lastind
#ncf['esom']['Mass']['P1'][dscen]#

#ncf['balance']['N']['to_water'][1]
#indices=np.where(~ncf['balance']['N']['to_water'][1].mask)
#lastind=max(indices[0])
#lastind#
#
#if type(ncf['balance']['N']['to_water'][1].mask)==np.ndarray:
#    print("loytyy mask info, voi kayttaa indeksointia tulosten koonnissa")
#
#type(ncf['esom']['Mass']['P1'][dscen].mask)
#if ~ncf['esom']['Mass']['P1'][dscen].mask:
#    print("array maks is false , cannot use indeksointi tulosten koonnistta, eka rivia ei voi laske muut lasketaan")

def get_ev_dscens(stand, scen):
    tapahtumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_tapahtumat.csv', encoding='latin1', sep=';', index_col=False)
    tdata = tapahtumat[(tapahtumat['kuvio']==stand)]
    tdata = tdata.reset_index(drop=True)       
    kunnoj= tdata[tdata['_2002']>-1]['_2002'].values
    ensih = tdata[tdata['_2003']>-1]['_2003'].values # ensiharvennus
    harv = tdata[tdata['_2004']>-1]['_2004'].values # harvennus
    paateh = tdata[tdata['_2005']>-1]['_2005'].values # pÃ¤Ã¤tehakkuu
    vuodet = np.unique(np.sort(np.concatenate((kunnoj, ensih, harv, paateh))))
    # Otetaan perkkisten vuosien poistettava vuosi talteen
    tuplat = pd.Series([], dtype=int)
    tup = 0
    # Tsekataan, onko tapahtumissa perak vuosia ja poistetaan niista jalkim
    for v in range(0,len(vuodet)-1):
       if vuodet[v+1]-vuodet[v]==1:
          tuplat[tup]=vuodet[v+1]
          tup = tup+1
          vuodet = np.delete(vuodet, [v+1], None)
       if v+1>=len(vuodet)-1:
          break
    n_mottifiles = len(vuodet) + 1
    #get kunnojitus vuosi, matsaa mink nx tiedoton tapahtumna kanssa? kanssa, get nx tiedosto
    #eka dsc2n 1 (50cm, sit kunojitus dscen 3), jos alle 10 v simuloinnin loppuu niin sit dscen 3 loppun asti muuten dscen 2 seuraavas n tiedostossa
    vuodet = vuodet[vuodet>5] #jos ongelmia eli ensiharvennus 0 vuotena tai mika tahans tapahtuma -sit heti dscen olis 90cm
    ind=np.where(vuodet[vuodet==kunnoj])
    dscen={'n0':1,'n1':1,'n2':1,'n3':1,'n4':1, 'n5':1, 'n6':1, 'n7':1,'n8':1,'n9':1}      
    if ind[0].size>0:
       eventid='n'+str(ind[0][0]+1)
       dscen[eventid]=3
       for i in range(2,9):
           dscen['n'+str(ind[0][0]+i)]=3
    return dscen

def peat_ini_end(ff, dscen):
    
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    
    wt = np.mean(ncf['strip']['dwtyr'][dscen,:, :], axis = 0)
    cols = np.shape(wt)[0]
    
    esoms =['L0L', 'L0W', 'LL', 'LW', 'FL', 'FW', 'H', 'P1', 'P2', 'P3'] 
    inipeat =  np.zeros(cols)
    for sto in esoms[7:]:
        inipeat += ncf['esom']['Mass'][sto][dscen,0, :]/10000.
    endpeat =  np.zeros(cols)
    for sto in esoms[7:]:
        indices=np.where(~ncf['esom']['Mass'][sto][dscen].mask)
        lastind=max(indices[0])
        endpeat += ncf['esom']['Mass'][sto][dscen,lastind, :]/10000. #pitais ottaa vika jossa on arvoja 
    inimor =  np.zeros(cols)
    for sto in esoms[:7]:
        inimor += ncf['esom']['Mass'][sto][dscen,0, :]/10000.
    endmor =  np.zeros(cols)
    for sto in esoms[:7]:
        indices=np.where(~ncf['esom']['Mass'][sto][dscen].mask)
        lastind=max(indices[0])
        endmor += ncf['esom']['Mass'][sto][dscen,lastind, :]/10000.
    
    return [np.mean(inipeat), np.mean(inimor), np.mean(endpeat), np.mean(endmor)]

def peat_ini_end_scaled(ff, dscen):
    
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    
    wt = np.mean(ncf['strip']['dwtyr'][dscen,:, :], axis = 0)
    cols = np.shape(wt)[0]
    
    esoms =['L0L', 'L0W', 'LL', 'LW', 'FL', 'FW', 'H', 'P1', 'P2', 'P3'] 
    inipeat =  np.zeros(cols)
    inipeatb =  np.zeros(cols)
    for sto in esoms[7:9]:
        inipeat += ncf['esom']['Mass'][sto][dscen,0, :]/10000.
    inipeatb = 120/190*ncf['esom']['Mass'][esoms[9]][dscen,0, :]/10000. #syvyys skaalaus eli 2.5m turvekerroksesta vain 1.8m tms mukaan, tai bottomista 60cm-180cm eli 120 koko 190cm sta
    inipeat=inipeat+inipeatb
    endpeat =  np.zeros(cols)
    endpeatb =  np.zeros(cols)
    for sto in esoms[7:9]:
        if type(ncf['esom']['Mass'][sto][dscen].mask)==np.ndarray:
            indices=np.where(~ncf['esom']['Mass'][sto][dscen].mask)
            lastind=max(indices[0])
            endpeat += ncf['esom']['Mass'][sto][dscen,lastind, :]/10000. #pitais ottaa vika jossa on arvoja 
        elif ~ncf['esom']['Mass'][sto][dscen].mask:
            print("esom P1 mask  is false")
            print(ncf['esom']['Mass'][sto][dscen,:, :]/10000)
            #indices=np.where(~ncf['esom']['Mass'][sto][dscen].mask)
            #lastind=max(indices[0])
            endpeat += ncf['esom']['Mass'][sto][dscen,-1, :]/10000. #pitais ottaa vika jossa on arvoja 
    if type(ncf['esom']['Mass'][esoms[9]][dscen].mask)==np.ndarray:
        indices=np.where(~ncf['esom']['Mass'][esoms[9]][dscen].mask)
        lastind=max(indices[0])
        endpeatb = 120/190*ncf['esom']['Mass'][esoms[9]][dscen,lastind, :]/10000. #pitais ottaa vika jossa on arvoja 
    elif ~ncf['esom']['Mass'][sto][dscen].mask:
        print("esom P1 mask  is false")
        print(ncf['esom']['Mass'][sto][dscen,:, :]/10000)
        print(ncf['esom']['Mass'][sto][dscen,-1, :]/10000)
        #indices=np.where(~ncf['esom']['Mass'][sto][dscen].mask)
        #lastind=max(indices[0])
        endpeatb = 120/190*ncf['esom']['Mass'][esoms[9]][dscen,-1, :]/10000. #pitais ottaa vika jossa on arvoja 
    endpeat=endpeat+endpeatb
    #indices=np.where(~ncf['esom']['Mass'][esoms[9]][dscen].mask)
    #lastind=max(indices[0])
    #endpeatb = 120/190*ncf['esom']['Mass'][esoms[9]][dscen,lastind, :]/10000. #pitais ottaa vika jossa on arvoja 
    #endpeat=endpeat+endpeatb
    inimor =  np.zeros(cols)
    for sto in esoms[:7]:
        inimor += ncf['esom']['Mass'][sto][dscen,0, :]/10000.
    endmor =  np.zeros(cols)
    for sto in esoms[:7]:
        if type(ncf['esom']['Mass'][sto][dscen].mask)==np.ndarray:
            indices=np.where(~ncf['esom']['Mass'][sto][dscen].mask)
            lastind=max(indices[0])
            endmor += ncf['esom']['Mass'][sto][dscen,lastind, :]/10000. #pitais ottaa vika jossa on arvoja 
        elif ~ncf['esom']['Mass'][sto][dscen].mask:
            print(sto, "esom  mask  is false")
            print(ncf['esom']['Mass'][sto][dscen,:, :]/10000)
            #indices=np.where(~ncf['esom']['Mass'][sto][dscen].mask)
            #lastind=max(indices[0])
            endmor += ncf['esom']['Mass'][sto][dscen,-1, :]/10000. #pitais ottaa vika jossa on arvoja 
        #indices=np.where(~ncf['esom']['Mass'][sto][dscen].mask)
        #lastind=max(indices[0])
        #endmor += ncf['esom']['Mass'][sto][dscen,lastind, :]/10000.
    
    return [np.mean(inipeat), np.mean(inimor), np.mean(endpeat), np.mean(endmor)]


def make_weather(ID, start_date,end_date, simyears, sourcefile=None):
    ID=0
    
    fmi=pd.read_csv(sourcefile, sep=';', header='infer', 
                    usecols=['OmaTunniste','Kunta','aika','longitude','latitude','t_mean','t_max','t_min','rainfall',\
                             'radiation','hpa'], encoding= 'ISO-8859-1')
    
    fmi_orig=fmi.copy()   
    time=pd.to_datetime(fmi['aika'],format='%Y%m%d')
    
    fmi.index=time
    fmi1=fmi[(fmi.index >= start_date) & (fmi.index <= end_date)]
    
    dates = pd.date_range(start_date, end_date).tolist()
    
    if len(dates) != len(fmi1):
        print(str(len(dates) - len(fmi1)) + ' days missing from forcing file, interpolated') 
            
    forcing = pd.DataFrame(index=dates, columns=[])
    forcing = forcing.merge(fmi1, how='outer', left_index=True, right_index=True)
    
    forcing['rainfall'] = forcing['rainfall'].fillna(0)
    forcing = forcing.fillna(method='ffill')  
    newenddateyear=end_date.year+simyears
    newend_date = datetime.datetime(newenddateyear,12,31) # for forc only
    climstart_date = datetime.datetime(end_date.year-simyears+1,1,1) # for forc only
    newstart_date = datetime.datetime(end_date.year+1,1,1) # for forc only
    fmi2=forcing[(forcing.index >= climstart_date ) & (forcing.index <= end_date)]
    dates_new = pd.date_range(newstart_date, newend_date).tolist()
    fmi2=fmi2.reset_index()
    
    forcing2 = pd.DataFrame(dates_new, columns=['indexn'])
    forcing2 = pd.concat([forcing2,fmi2], axis="columns")
    timen=pd.to_datetime(forcing2['indexn'],format='%Y%m%d')
    forcing2.index=timen
    forcing2=forcing2.drop(columns=['index', 'indexn']) 
    forcing_f=pd.concat([forcing, forcing2],axis="rows")
    #forcing2['Prec'] = forcing2['Prec'].fillna(0)
    #forcing2 = forcing2.fillna(method='ffill')
    forcing_f=forcing_f.reset_index()
    for i in range(len(forcing_f)):
        forcing_f.loc[i,'aika']=forcing_f['index'][i].strftime("%Y%m%d")
    forcing_f['aika']=forcing_f['aika'].astype(int) 
    
    forcing_f=forcing_f.drop(columns=['index']) 
    forcing_f.to_csv(sourcefile[:-4]+'cont.csv', sep=';', encoding= 'ISO-8859-1', index=False)
    return forcing_f 

def get_hiili(scenario):
    resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'
    files=glob(resfol+scenario+'/*.nc')
    stands=[i[(len(resfol)+len(scenario)+1):-6] for i in files]
    stands = [*set(stands)]
    stands.sort()
    storages=[]
    endpeats=[]
    endmors=[]
    for f in stands:
        dscen=get_ev_dscens(f, scenario)
        events=[i[(len(resfol)+len(scenario)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scenario+'/'+f)]
        events.sort()
        for ev in events:
            sto_values=peat_ini_end(resfol+scenario+"/"+f+"_"+ev+".nc", dscen[ev])
            storages.append(sto_values)
        endpeat= storages[0][2]
        endmor= storages[0][3]
        for i in range(1,len(events)):
            endpeat = endpeat-(storages[i][0]-storages[i][2])
            endmor = endmor-(storages[i][1]-storages[i][3])     
        endpeats.append(endpeat)    
        endmors.append(endmor)
        #print(stands.index(f),"of",len(stands), "done")    
    return [endpeats, endmors]

def get_hiili_scaled(scenario):
    resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'
    files=glob(resfol+scenario+'/*.nc')
    stands=[i[(len(resfol)+len(scenario)+1):-6] for i in files]
    stands = [*set(stands)]
    stands.sort()
    storages=[]
    endpeats=[]
    endmors=[]
    for f in stands:
        dscen=get_ev_dscens(f, scenario)
        events=[i[(len(resfol)+len(scenario)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scenario+'/'+f)]
        events.sort()
        for ev in events:
            sto_values=peat_ini_end_scaled(resfol+scenario+"/"+f+"_"+ev+".nc", dscen[ev])
            storages.append(sto_values)
        endpeat= storages[0][2]
        endmor= storages[0][3]
        for i in range(1,len(events)):
            endpeat = endpeat-(storages[i][0]-storages[i][2])
            endmor = endmor-(storages[i][1]-storages[i][3])     
        endpeats.append(endpeat)    
        endmors.append(endmor)
        print(stands.index(f),"of",len(stands), "done")    
    return [endpeats, endmors]

def get_hiili_scaled_exp(scenario):
    resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'
    files=glob(resfol+scenario+'/*.nc')
    stands=[i[(len(resfol)+len(scenario)+1):-6] for i in files]
    stands = [*set(stands)]
    stands.remove('43069586')
    stands.sort()
    storages=[]
    endpeats=[]
    endmors=[]
    for f in stands:
        dscen=get_ev_dscens(f, scenario)
        events=[i[(len(resfol)+len(scenario)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scenario+'/'+f)]
        events.sort()
        for ev in events:
            sto_values=peat_ini_end_scaled(resfol+scenario+"/"+f+"_"+ev+".nc", dscen[ev])
            storages.append(sto_values)
        endpeat= storages[0][2]
        endmor= storages[0][3]
        for i in range(1,len(events)):
            endpeat = endpeat-(storages[i][0]-storages[i][2])
            endmor = endmor-(storages[i][1]-storages[i][3])     
        endpeats.append(endpeat)    
        endmors.append(endmor)
        print(stands.index(f),"of",len(stands), "done")    
    return [endpeats, endmors]

def get_hiili_scaled_exp2(scenario):
    resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'
    files=glob(resfol+scenario+'/*.nc')
    stands=[i[(len(resfol)+len(scenario)+1):-6] for i in files]
    stands = [*set(stands)]
    #stands.remove('43069586')
    stands.remove('35394412')
    stands.sort()
    storages=[]
    endpeats=[]
    endmors=[]
    for f in stands:
        dscen=get_ev_dscens(f, scenario)
        events=[i[(len(resfol)+len(scenario)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scenario+'/'+f)]
        events.sort()
        for ev in events:
            sto_values=peat_ini_end_scaled(resfol+scenario+"/"+f+"_"+ev+".nc", dscen[ev])
            storages.append(sto_values)
        endpeat= storages[0][2]
        endmor= storages[0][3]
        for i in range(1,len(events)):
            endpeat = endpeat-(storages[i][0]-storages[i][2])
            endmor = endmor-(storages[i][1]-storages[i][3])     
        endpeats.append(endpeat)    
        endmors.append(endmor)
        print(stands.index(f),"of",len(stands), "done")    
    return [endpeats, endmors]


def get_hiili_scaled_exp_sorvas(scenario):
    resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'
    files=glob(resfol+scenario+'/*.nc')
    stands=[i[(len(resfol)+len(scenario)+1):-6] for i in files]
    stands = [*set(stands)]
    #stands.remove('43069586')
    #stands.remove('35394412')
    prob_stands=['24298445','24301359','24301365','24301373','24301452','29579511','35393481','35723701','35723705','35725748','35725767','35725769','35725773','35725774','35725824','35725834','35725949','36109336','40699198','41267888','43034912','43095058','43095060','43095061','43095840','43096868','43097183','43097204']
    for p in prob_stands:
        try:
            stands.remove(p)
        except ValueError:
            pass
    stands.sort()
    storages=[]
    endpeats=[]
    endmors=[]
    for f in stands:
        dscen=get_ev_dscens(f, scenario)
        events=[i[(len(resfol)+len(scenario)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scenario+'/'+f)]
        events.sort()
        for ev in events:
            sto_values=peat_ini_end_scaled(resfol+scenario+"/"+f+"_"+ev+".nc", dscen[ev])
            storages.append(sto_values)
        endpeat= storages[0][2]
        endmor= storages[0][3]
        for i in range(1,len(events)):
            endpeat = endpeat-(storages[i][0]-storages[i][2])
            endmor = endmor-(storages[i][1]-storages[i][3])     
        endpeats.append(endpeat)    
        endmors.append(endmor)
        print(stands.index(f),"of",len(stands), "done")    
    return [endpeats, endmors]

def get_hiili_scaled_prex(scenario, problem_stand_list):
    resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'
    files=glob(resfol+scenario+'/*.nc')
    stands=[i[(len(resfol)+len(scenario)+1):-6] for i in files]
    stands = [*set(stands)]
    #stands.remove('43069586')
    #stands.remove('35394412')
    #prob_stands=['24298445','24301359','24301365','24301373','24301452','29579511','35393481','35723701','35723705','35725748','35725767','35725769','35725773','35725774','35725824','35725834','35725949','36109336','40699198','41267888','43034912','43095058','43095060','43095061','43095840','43096868','43097183','43097204']
    prob_stands=problem_stand_list
    for p in prob_stands:
        try:
            stands.remove(p)
        except ValueError:
            pass
    stands.sort()
    storages=[]
    endpeats=[]
    endmors=[]
    for f in stands:
        dscen=get_ev_dscens(f, scenario)
        events=[i[(len(resfol)+len(scenario)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scenario+'/'+f)]
        events.sort()
        for ev in events:
            sto_values=peat_ini_end_scaled(resfol+scenario+"/"+f+"_"+ev+".nc", dscen[ev])
            storages.append(sto_values)
        endpeat= storages[0][2]
        endmor= storages[0][3]
        for i in range(1,len(events)):
            endpeat = endpeat-(storages[i][0]-storages[i][2])
            endmor = endmor-(storages[i][1]-storages[i][3])     
        endpeats.append(endpeat)    
        endmors.append(endmor)
        print(stands.index(f),"of",len(stands), "done")    
    return [endpeats, endmors]



#scenario=['Sorvasranta_BAU_A', 'Sorvasranta_BAU_B']
#scenarios=['Sorvasranta_BAU_A', 'Sorvasranta_BAU_B', 'Sorvasranta_BIO_A', 'Sorvasranta_BIO_B','Sorvasranta_HII_A', 'Sorvasranta_HII_B']

resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'


#files=os.listdir(resfol+scenario[0]+'/')
#files=os.listdir(resfol+scenario[0]+'/')

def get_NP(scen):
    files=glob(resfol+scen+'/*.nc')
    stands=[i[(len(resfol)+len(scen)+1):-6] for i in files]
    stands = [*set(stands)]
    stands.sort()
    ncfs={}
    n=[]
    p=[]
    nsum=[]
    psum=[]
    areas=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', sep=';', encoding='latin1')
    for f in stands:
        #print(f, "alkaa", stands.index(f))
        ncfs={}
        events=[i[(len(resfol)+len(scen)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scen+'/'+f)]
        #print(f, events)
        standha= areas['ala'][areas['kuvio']==int(f)].values[0]
        dscen=get_ev_dscens(f, scen)
        if 'n'+str(0) in events:
            ncfs["n{0}".format(0)] = Dataset(resfol+scen+'/'+f+'_n'+str(0)+'.nc', mode='r') 
            n0ndata=np.mean(ncfs["n{0}".format(0)]['balance']['N']['to_water'][dscen['n0']][:,1:-1],axis=1)[1:]
            n0ndata=n0ndata[n0ndata.mask==False].data
            n0pdata=np.mean(ncfs["n{0}".format(0)]['balance']['P']['to_water'][dscen['n0']][:,1:-1],axis=1)[1:]
            n0pdata=n0pdata[n0pdata.mask==False].data
            ncfs["n{0}".format(0)].close
            for j in range(1,len(events)):
                ncfs["n{0}".format(j)] = Dataset(resfol+scen+'/'+f+'_n'+str(j)+'.nc', mode='r') 
                nxndata=np.mean(ncfs["n{0}".format(j)]['balance']['N']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:] 
                nxndata=nxndata[nxndata.mask==False].data
                n0ndata=np.concatenate((n0ndata,nxndata),axis=0)
                nxpdata=np.mean(ncfs["n{0}".format(j)]['balance']['P']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxpdata=nxpdata[nxpdata.mask==False].data
                n0pdata=np.concatenate((n0pdata,nxpdata),axis=0)
                ncfs["n{0}".format(j)].close        
            n0ndata_kg=standha*n0ndata
            n0pdata_kg=standha*n0pdata
            if len(n0ndata)>0:
                n.append(sum(n0ndata_kg)/len(n0ndata))
                p.append(sum(n0pdata_kg)/len(n0pdata))
            else:
                n.append(np.nan)
                p.append(np.nan)
            nsum.append(sum(n0ndata_kg))
            psum.append(sum(n0pdata_kg))
            #print(stands.index(f),"of",len(stands), "done")
        else:  
            ncfs["n{0}".format(1)] = Dataset(resfol+scen+'/'+f+'_n'+str(1)+'.nc', mode='r') 
            n0ndata=np.mean(ncfs["n{0}".format(1)]['balance']['N']['to_water'][dscen["n{0}".format(1)]][:,1:-1],axis=1)[1:] #kg/ha/yr mean over stand per year
            n0ndata=n0ndata[n0ndata.mask==False].data
            n0pdata=np.mean(ncfs["n{0}".format(1)]['balance']['P']['to_water'][dscen["n{0}".format(1)]][:,1:-1],axis=1)[1:]
            n0pdata=n0pdata[n0pdata.mask==False].data
            ncfs["n{0}".format(1)].close          
            for j in range(2,len(events)+1):
                ncfs["n{0}".format(j)] = Dataset(resfol+scen+'/'+f+'_n'+str(j)+'.nc', mode='r') 
                nxndata=np.mean(ncfs["n{0}".format(j)]['balance']['N']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxndata=nxndata[nxndata.mask==False].data
                n0ndata=np.concatenate((n0ndata,nxndata),axis=0)
                nxpdata=np.mean(ncfs["n{0}".format(j)]['balance']['P']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxpdata=nxpdata[nxpdata.mask==False].data
                n0pdata=np.concatenate((n0pdata,nxpdata),axis=0)
                ncfs["n{0}".format(j)].close        
            n0ndata_kg=standha*n0ndata
            n0pdata_kg=standha*n0pdata
            if len(n0ndata)>0:
                n.append(sum(n0ndata_kg)/len(n0ndata))
                p.append(sum(n0pdata_kg)/len(n0pdata))
            else:
                n.append(np.nan)
                p.append(np.nan)
            nsum.append(sum(n0ndata_kg))
            psum.append(sum(n0pdata_kg))
            #n.append(sum(n0ndata)/len(n0ndata))
            #p.append(sum(n0pdata)/len(n0pdata))
            #print(stands.index(f),"of",len(stands), "done")
        #print(f, "loppuu", stands.index(f))
    return {'nsum':nsum, 'psum':psum, 'n':n, 'p':p}      

def get_NP_exp(scen):
    files=glob(resfol+scen+'/*.nc')
    stands=[i[(len(resfol)+len(scen)+1):-6] for i in files]
    stands = [*set(stands)]
    stands.remove('43069586')
    #stands.remove('35394412')
    stands.sort()
    ncfs={}
    n=[]
    p=[]
    nsum=[]
    psum=[]
    areas=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', sep=';', encoding='latin1')
    for f in stands:
        #print(f, "alkaa", stands.index(f))
        ncfs={}
        events=[i[(len(resfol)+len(scen)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scen+'/'+f)]
        #print(f, events)
        standha= areas['ala'][areas['kuvio']==int(f)].values[0]
        dscen=get_ev_dscens(f, scen)
        if 'n'+str(0) in events:
            ncfs["n{0}".format(0)] = Dataset(resfol+scen+'/'+f+'_n'+str(0)+'.nc', mode='r') 
            n0ndata=np.mean(ncfs["n{0}".format(0)]['balance']['N']['to_water'][dscen['n0']][:,1:-1],axis=1)[1:]
            n0ndata=n0ndata[n0ndata.mask==False].data
            n0pdata=np.mean(ncfs["n{0}".format(0)]['balance']['P']['to_water'][dscen['n0']][:,1:-1],axis=1)[1:]
            n0pdata=n0pdata[n0pdata.mask==False].data
            ncfs["n{0}".format(0)].close
            for j in range(1,len(events)):
                ncfs["n{0}".format(j)] = Dataset(resfol+scen+'/'+f+'_n'+str(j)+'.nc', mode='r') 
                nxndata=np.mean(ncfs["n{0}".format(j)]['balance']['N']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:] 
                nxndata=nxndata[nxndata.mask==False].data
                n0ndata=np.concatenate((n0ndata,nxndata),axis=0)
                nxpdata=np.mean(ncfs["n{0}".format(j)]['balance']['P']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxpdata=nxpdata[nxpdata.mask==False].data
                n0pdata=np.concatenate((n0pdata,nxpdata),axis=0)
                ncfs["n{0}".format(j)].close        
            n0ndata_kg=standha*n0ndata
            n0pdata_kg=standha*n0pdata
            if len(n0ndata)>0:
                n.append(sum(n0ndata_kg)/len(n0ndata))
                p.append(sum(n0pdata_kg)/len(n0pdata))
            else:
                n.append(np.nan)
                p.append(np.nan)
            nsum.append(sum(n0ndata_kg))
            psum.append(sum(n0pdata_kg))
            #print(stands.index(f),"of",len(stands), "done")
        else:  
            ncfs["n{0}".format(1)] = Dataset(resfol+scen+'/'+f+'_n'+str(1)+'.nc', mode='r') 
            n0ndata=np.mean(ncfs["n{0}".format(1)]['balance']['N']['to_water'][dscen["n{0}".format(1)]][:,1:-1],axis=1)[1:] #kg/ha/yr mean over stand per year
            n0ndata=n0ndata[n0ndata.mask==False].data
            n0pdata=np.mean(ncfs["n{0}".format(1)]['balance']['P']['to_water'][dscen["n{0}".format(1)]][:,1:-1],axis=1)[1:]
            n0pdata=n0pdata[n0pdata.mask==False].data
            ncfs["n{0}".format(1)].close          
            for j in range(2,len(events)+1):
                ncfs["n{0}".format(j)] = Dataset(resfol+scen+'/'+f+'_n'+str(j)+'.nc', mode='r') 
                nxndata=np.mean(ncfs["n{0}".format(j)]['balance']['N']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxndata=nxndata[nxndata.mask==False].data
                n0ndata=np.concatenate((n0ndata,nxndata),axis=0)
                nxpdata=np.mean(ncfs["n{0}".format(j)]['balance']['P']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxpdata=nxpdata[nxpdata.mask==False].data
                n0pdata=np.concatenate((n0pdata,nxpdata),axis=0)
                ncfs["n{0}".format(j)].close        
            n0ndata_kg=standha*n0ndata
            n0pdata_kg=standha*n0pdata
            if len(n0ndata)>0:
                n.append(sum(n0ndata_kg)/len(n0ndata))
                p.append(sum(n0pdata_kg)/len(n0pdata))
            else:
                n.append(np.nan)
                p.append(np.nan)
            nsum.append(sum(n0ndata_kg))
            psum.append(sum(n0pdata_kg))
            #n.append(sum(n0ndata)/len(n0ndata))
            #p.append(sum(n0pdata)/len(n0pdata))
            #print(stands.index(f),"of",len(stands), "done")
        #print(f, "loppuu", stands.index(f))
    return {'nsum':nsum, 'psum':psum, 'n':n, 'p':p}      

def get_NP_exp2(scen):
    files=glob(resfol+scen+'/*.nc')
    stands=[i[(len(resfol)+len(scen)+1):-6] for i in files]
    stands = [*set(stands)]
    #stands.remove('43069586')
    stands.remove('35394412')
    stands.sort()
    ncfs={}
    n=[]
    p=[]
    nsum=[]
    psum=[]
    areas=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', sep=';', encoding='latin1')
    for f in stands:
        #print(f, "alkaa", stands.index(f))
        ncfs={}
        events=[i[(len(resfol)+len(scen)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scen+'/'+f)]
        #print(f, events)
        standha= areas['ala'][areas['kuvio']==int(f)].values[0]
        dscen=get_ev_dscens(f, scen)
        if 'n'+str(0) in events:
            ncfs["n{0}".format(0)] = Dataset(resfol+scen+'/'+f+'_n'+str(0)+'.nc', mode='r') 
            n0ndata=np.mean(ncfs["n{0}".format(0)]['balance']['N']['to_water'][dscen['n0']][:,1:-1],axis=1)[1:]
            n0ndata=n0ndata[n0ndata.mask==False].data
            n0pdata=np.mean(ncfs["n{0}".format(0)]['balance']['P']['to_water'][dscen['n0']][:,1:-1],axis=1)[1:]
            n0pdata=n0pdata[n0pdata.mask==False].data
            ncfs["n{0}".format(0)].close
            for j in range(1,len(events)):
                ncfs["n{0}".format(j)] = Dataset(resfol+scen+'/'+f+'_n'+str(j)+'.nc', mode='r') 
                nxndata=np.mean(ncfs["n{0}".format(j)]['balance']['N']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:] 
                nxndata=nxndata[nxndata.mask==False].data
                n0ndata=np.concatenate((n0ndata,nxndata),axis=0)
                nxpdata=np.mean(ncfs["n{0}".format(j)]['balance']['P']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxpdata=nxpdata[nxpdata.mask==False].data
                n0pdata=np.concatenate((n0pdata,nxpdata),axis=0)
                ncfs["n{0}".format(j)].close        
            n0ndata_kg=standha*n0ndata
            n0pdata_kg=standha*n0pdata
            if len(n0ndata)>0:
                n.append(sum(n0ndata_kg)/len(n0ndata))
                p.append(sum(n0pdata_kg)/len(n0pdata))
            else:
                n.append(np.nan)
                p.append(np.nan)
            nsum.append(sum(n0ndata_kg))
            psum.append(sum(n0pdata_kg))
            #print(stands.index(f),"of",len(stands), "done")
        else:  
            ncfs["n{0}".format(1)] = Dataset(resfol+scen+'/'+f+'_n'+str(1)+'.nc', mode='r') 
            n0ndata=np.mean(ncfs["n{0}".format(1)]['balance']['N']['to_water'][dscen["n{0}".format(1)]][:,1:-1],axis=1)[1:] #kg/ha/yr mean over stand per year
            n0ndata=n0ndata[n0ndata.mask==False].data
            n0pdata=np.mean(ncfs["n{0}".format(1)]['balance']['P']['to_water'][dscen["n{0}".format(1)]][:,1:-1],axis=1)[1:]
            n0pdata=n0pdata[n0pdata.mask==False].data
            ncfs["n{0}".format(1)].close          
            for j in range(2,len(events)+1):
                ncfs["n{0}".format(j)] = Dataset(resfol+scen+'/'+f+'_n'+str(j)+'.nc', mode='r') 
                nxndata=np.mean(ncfs["n{0}".format(j)]['balance']['N']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxndata=nxndata[nxndata.mask==False].data
                n0ndata=np.concatenate((n0ndata,nxndata),axis=0)
                nxpdata=np.mean(ncfs["n{0}".format(j)]['balance']['P']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxpdata=nxpdata[nxpdata.mask==False].data
                n0pdata=np.concatenate((n0pdata,nxpdata),axis=0)
                ncfs["n{0}".format(j)].close        
            n0ndata_kg=standha*n0ndata
            n0pdata_kg=standha*n0pdata
            if len(n0ndata)>0:
                n.append(sum(n0ndata_kg)/len(n0ndata))
                p.append(sum(n0pdata_kg)/len(n0pdata))
            else:
                n.append(np.nan)
                p.append(np.nan)
            nsum.append(sum(n0ndata_kg))
            psum.append(sum(n0pdata_kg))
            #n.append(sum(n0ndata)/len(n0ndata))
            #p.append(sum(n0pdata)/len(n0pdata))
            #print(stands.index(f),"of",len(stands), "done")
        #print(f, "loppuu", stands.index(f))
    return {'nsum':nsum, 'psum':psum, 'n':n, 'p':p}      

def get_NP_sorvas(scen):
    files=glob(resfol+scen+'/*.nc')
    stands=[i[(len(resfol)+len(scen)+1):-6] for i in files]
    stands = [*set(stands)]
    #stands.remove('43069586')
    prob_stands=['24298445','24301359','24301365','24301373','24301452','29579511','35393481','35723701','35723705','35725748','35725767','35725769','35725773','35725774','35725824','35725834','35725949','36109336','40699198','41267888','43034912','43095058','43095060','43095061','43095840','43096868','43097183','43097204']
    for p in prob_stands:
        try:
            stands.remove(p)
        except ValueError:
            pass
    #stands.remove('24298445','24301359','24301365','24301373','24301452','29579511','35393481','35723701','35723705','35725748','35725767','35725769','35725773','35725774','35725824','35725834','35725949','36109336','40699198','41267888','43034912','43095058','43095060','43095061','43095840','43096868','43097183','43097204')
    stands.sort()
    ncfs={}
    n=[]
    p=[]
    nsum=[]
    psum=[]
    areas=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', sep=';', encoding='latin1')
    for f in stands:
        #print(f, "alkaa", stands.index(f))
        ncfs={}
        events=[i[(len(resfol)+len(scen)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scen+'/'+f)]
        #print(f, events)
        standha= areas['ala'][areas['kuvio']==int(f)].values[0]
        dscen=get_ev_dscens(f, scen)
        if 'n'+str(0) in events:
            ncfs["n{0}".format(0)] = Dataset(resfol+scen+'/'+f+'_n'+str(0)+'.nc', mode='r') 
            n0ndata=np.mean(ncfs["n{0}".format(0)]['balance']['N']['to_water'][dscen['n0']][:,1:-1],axis=1)[1:]
            n0ndata=n0ndata[n0ndata.mask==False].data
            n0pdata=np.mean(ncfs["n{0}".format(0)]['balance']['P']['to_water'][dscen['n0']][:,1:-1],axis=1)[1:]
            n0pdata=n0pdata[n0pdata.mask==False].data
            ncfs["n{0}".format(0)].close
            for j in range(1,len(events)):
                ncfs["n{0}".format(j)] = Dataset(resfol+scen+'/'+f+'_n'+str(j)+'.nc', mode='r') 
                nxndata=np.mean(ncfs["n{0}".format(j)]['balance']['N']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:] 
                nxndata=nxndata[nxndata.mask==False].data
                n0ndata=np.concatenate((n0ndata,nxndata),axis=0)
                nxpdata=np.mean(ncfs["n{0}".format(j)]['balance']['P']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxpdata=nxpdata[nxpdata.mask==False].data
                n0pdata=np.concatenate((n0pdata,nxpdata),axis=0)
                ncfs["n{0}".format(j)].close        
            n0ndata_kg=standha*n0ndata
            n0pdata_kg=standha*n0pdata
            if len(n0ndata)>0:
                n.append(sum(n0ndata_kg)/len(n0ndata))
                p.append(sum(n0pdata_kg)/len(n0pdata))
            else:
                n.append(np.nan)
                p.append(np.nan)
            nsum.append(sum(n0ndata_kg))
            psum.append(sum(n0pdata_kg))
            #print(stands.index(f),"of",len(stands), "done")
        else:  
            ncfs["n{0}".format(1)] = Dataset(resfol+scen+'/'+f+'_n'+str(1)+'.nc', mode='r') 
            n0ndata=np.mean(ncfs["n{0}".format(1)]['balance']['N']['to_water'][dscen["n{0}".format(1)]][:,1:-1],axis=1)[1:] #kg/ha/yr mean over stand per year
            n0ndata=n0ndata[n0ndata.mask==False].data
            n0pdata=np.mean(ncfs["n{0}".format(1)]['balance']['P']['to_water'][dscen["n{0}".format(1)]][:,1:-1],axis=1)[1:]
            n0pdata=n0pdata[n0pdata.mask==False].data
            ncfs["n{0}".format(1)].close          
            for j in range(2,len(events)+1):
                ncfs["n{0}".format(j)] = Dataset(resfol+scen+'/'+f+'_n'+str(j)+'.nc', mode='r') 
                nxndata=np.mean(ncfs["n{0}".format(j)]['balance']['N']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxndata=nxndata[nxndata.mask==False].data
                n0ndata=np.concatenate((n0ndata,nxndata),axis=0)
                nxpdata=np.mean(ncfs["n{0}".format(j)]['balance']['P']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxpdata=nxpdata[nxpdata.mask==False].data
                n0pdata=np.concatenate((n0pdata,nxpdata),axis=0)
                ncfs["n{0}".format(j)].close        
            n0ndata_kg=standha*n0ndata
            n0pdata_kg=standha*n0pdata
            if len(n0ndata)>0:
                n.append(sum(n0ndata_kg)/len(n0ndata))
                p.append(sum(n0pdata_kg)/len(n0pdata))
            else:
                n.append(np.nan)
                p.append(np.nan)
            nsum.append(sum(n0ndata_kg))
            psum.append(sum(n0pdata_kg))
            #n.append(sum(n0ndata)/len(n0ndata))
            #p.append(sum(n0pdata)/len(n0pdata))
            #print(stands.index(f),"of",len(stands), "done")
        #print(f, "loppuu", stands.index(f))
    return {'nsum':nsum, 'psum':psum, 'n':n, 'p':p}      

def get_NP_prex(scen,problem_stand_list):
    files=glob(resfol+scen+'/*.nc')
    stands=[i[(len(resfol)+len(scen)+1):-6] for i in files]
    stands = [*set(stands)]
    #stands.remove('43069586')
    #prob_stands=['24298445','24301359','24301365','24301373','24301452','29579511','35393481','35723701','35723705','35725748','35725767','35725769','35725773','35725774','35725824','35725834','35725949','36109336','40699198','41267888','43034912','43095058','43095060','43095061','43095840','43096868','43097183','43097204']
    prob_stands=problem_stand_list
    for p in prob_stands:
        try:
            stands.remove(p)
        except ValueError:
            pass
    #stands.remove('24298445','24301359','24301365','24301373','24301452','29579511','35393481','35723701','35723705','35725748','35725767','35725769','35725773','35725774','35725824','35725834','35725949','36109336','40699198','41267888','43034912','43095058','43095060','43095061','43095840','43096868','43097183','43097204')
    stands.sort()
    ncfs={}
    n=[]
    p=[]
    nsum=[]
    psum=[]
    areas=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', sep=';', encoding='latin1')
    for f in stands:
        #print(f, "alkaa", stands.index(f))
        ncfs={}
        events=[i[(len(resfol)+len(scen)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scen+'/'+f)]
        #print(f, events)
        standha= areas['ala'][areas['kuvio']==int(f)].values[0]
        dscen=get_ev_dscens(f, scen)
        if 'n'+str(0) in events:
            ncfs["n{0}".format(0)] = Dataset(resfol+scen+'/'+f+'_n'+str(0)+'.nc', mode='r') 
            n0ndata=np.mean(ncfs["n{0}".format(0)]['balance']['N']['to_water'][dscen['n0']][:,1:-1],axis=1)[1:]
            n0ndata=n0ndata[n0ndata.mask==False].data
            n0pdata=np.mean(ncfs["n{0}".format(0)]['balance']['P']['to_water'][dscen['n0']][:,1:-1],axis=1)[1:]
            n0pdata=n0pdata[n0pdata.mask==False].data
            ncfs["n{0}".format(0)].close
            for j in range(1,len(events)):
                ncfs["n{0}".format(j)] = Dataset(resfol+scen+'/'+f+'_n'+str(j)+'.nc', mode='r') 
                nxndata=np.mean(ncfs["n{0}".format(j)]['balance']['N']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:] 
                nxndata=nxndata[nxndata.mask==False].data
                n0ndata=np.concatenate((n0ndata,nxndata),axis=0)
                nxpdata=np.mean(ncfs["n{0}".format(j)]['balance']['P']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxpdata=nxpdata[nxpdata.mask==False].data
                n0pdata=np.concatenate((n0pdata,nxpdata),axis=0)
                ncfs["n{0}".format(j)].close        
            n0ndata_kg=standha*n0ndata
            n0pdata_kg=standha*n0pdata
            if len(n0ndata)>0:
                n.append(sum(n0ndata_kg)/len(n0ndata))
                p.append(sum(n0pdata_kg)/len(n0pdata))
            else:
                n.append(np.nan)
                p.append(np.nan)
            nsum.append(sum(n0ndata_kg))
            psum.append(sum(n0pdata_kg))
            #print(stands.index(f),"of",len(stands), "done")
        else:  
            ncfs["n{0}".format(1)] = Dataset(resfol+scen+'/'+f+'_n'+str(1)+'.nc', mode='r') 
            n0ndata=np.mean(ncfs["n{0}".format(1)]['balance']['N']['to_water'][dscen["n{0}".format(1)]][:,1:-1],axis=1)[1:] #kg/ha/yr mean over stand per year
            n0ndata=n0ndata[n0ndata.mask==False].data
            n0pdata=np.mean(ncfs["n{0}".format(1)]['balance']['P']['to_water'][dscen["n{0}".format(1)]][:,1:-1],axis=1)[1:]
            n0pdata=n0pdata[n0pdata.mask==False].data
            ncfs["n{0}".format(1)].close          
            for j in range(2,len(events)+1):
                ncfs["n{0}".format(j)] = Dataset(resfol+scen+'/'+f+'_n'+str(j)+'.nc', mode='r') 
                nxndata=np.mean(ncfs["n{0}".format(j)]['balance']['N']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxndata=nxndata[nxndata.mask==False].data
                n0ndata=np.concatenate((n0ndata,nxndata),axis=0)
                nxpdata=np.mean(ncfs["n{0}".format(j)]['balance']['P']['to_water'][dscen["n{0}".format(j)]][:,1:-1],axis=1)[1:]
                nxpdata=nxpdata[nxpdata.mask==False].data
                n0pdata=np.concatenate((n0pdata,nxpdata),axis=0)
                ncfs["n{0}".format(j)].close        
            n0ndata_kg=standha*n0ndata
            n0pdata_kg=standha*n0pdata
            if len(n0ndata)>0:
                n.append(sum(n0ndata_kg)/len(n0ndata))
                p.append(sum(n0pdata_kg)/len(n0pdata))
            else:
                n.append(np.nan)
                p.append(np.nan)
            nsum.append(sum(n0ndata_kg))
            psum.append(sum(n0pdata_kg))
            #n.append(sum(n0ndata)/len(n0ndata))
            #p.append(sum(n0pdata)/len(n0pdata))
            #print(stands.index(f),"of",len(stands), "done")
        #print(f, "loppuu", stands.index(f))
    return {'nsum':nsum, 'psum':psum, 'n':n, 'p':p}      

def get_areas(scen):
    files=glob(resfol+scen+'/*.nc')
    stands=[i[(len(resfol)+len(scen)+1):-6] for i in files]
    stands = [*set(stands)]
    stands.sort()
    areas=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', sep=';', encoding='latin1')
    standareas=[]
    for f in stands:
        standha= areas['ala'][areas['kuvio']==int(f)].values[0]
        standareas.append(standha)
    return standareas

def get_areas_exp(scen):
    files=glob(resfol+scen+'/*.nc')
    stands=[i[(len(resfol)+len(scen)+1):-6] for i in files]
    stands = [*set(stands)]
    stands.remove('43069586')
    stands.sort()
    areas=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', sep=';', encoding='latin1')
    standareas=[]
    for f in stands:
        standha= areas['ala'][areas['kuvio']==int(f)].values[0]
        standareas.append(standha)
    return standareas

def get_areas_exp2(scen):
    files=glob(resfol+scen+'/*.nc')
    stands=[i[(len(resfol)+len(scen)+1):-6] for i in files]
    stands = [*set(stands)]
    #stands.remove('43069586')
    stands.remove('35394412')
    stands.sort()
    areas=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', sep=';', encoding='latin1')
    standareas=[]
    for f in stands:
        standha= areas['ala'][areas['kuvio']==int(f)].values[0]
        standareas.append(standha)
    return standareas

def get_areas_sorvas(scen):
    files=glob(resfol+scen+'/*.nc')
    stands=[i[(len(resfol)+len(scen)+1):-6] for i in files]
    stands = [*set(stands)]
    #stands.remove('43069586')
    prob_stands=['24298445','24301359','24301365','24301373','24301452','29579511','35393481','35723701','35723705','35725748','35725767','35725769','35725773','35725774','35725824','35725834','35725949','36109336','40699198','41267888','43034912','43095058','43095060','43095061','43095840','43096868','43097183','43097204']
    for p in prob_stands:
        try:
            stands.remove(p)
        except ValueError:
            pass
    #stands.remove('24298445','24301359','24301365','24301373','24301452','29579511','35393481','35723701','35723705','35725748','35725767','35725769','35725773','35725774','35725824','35725834','35725949','36109336','40699198','41267888','43034912','43095058','43095060','43095061','43095840','43096868','43097183','43097204')
    stands.sort()
    areas=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', sep=';', encoding='latin1')
    standareas=[]
    for f in stands:
        standha= areas['ala'][areas['kuvio']==int(f)].values[0]
        standareas.append(standha)
    return standareas

def get_areas_prex(scen, problem_stand_list):
    files=glob(resfol+scen+'/*.nc')
    stands=[i[(len(resfol)+len(scen)+1):-6] for i in files]
    stands = [*set(stands)]
    #stands.remove('43069586')
    #prob_stands=['24298445','24301359','24301365','24301373','24301452','29579511','35393481','35723701','35723705','35725748','35725767','35725769','35725773','35725774','35725824','35725834','35725949','36109336','40699198','41267888','43034912','43095058','43095060','43095061','43095840','43096868','43097183','43097204']
    prob_stands=problem_stand_list
    for p in prob_stands:
        try:
            stands.remove(p)
        except ValueError:
            pass
    #stands.remove('24298445','24301359','24301365','24301373','24301452','29579511','35393481','35723701','35723705','35725748','35725767','35725769','35725773','35725774','35725824','35725834','35725949','36109336','40699198','41267888','43034912','43095058','43095060','43095061','43095840','43096868','43097183','43097204')
    stands.sort()
    areas=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', sep=';', encoding='latin1')
    standareas=[]
    for f in stands:
        standha= areas['ala'][areas['kuvio']==int(f)].values[0]
        standareas.append(standha)
    return standareas


def checkres(scenario):
    resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'
    files=glob(resfol+scenario+'/*.nc')
    stands=[i[(len(resfol)+len(scenario)+1):-6] for i in files]
    stands = [*set(stands)]
    #stands.remove('43069586')
    stands.sort()
    problems=[]
    problem_stands=[]
    for f in stands:
        print(f, "alkaa")
        dscen=get_ev_dscens(f, scenario)
        events=[i[(len(resfol)+len(scenario)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scenario+'/'+f)]
        events.sort()
        for ev in events:
            ff=resfol+scenario+"/"+f+"_"+ev+".nc"
            ncf=Dataset(ff, mode='r')
            wt = ncf['strip']['dwtyr'][dscen[ev],:, :]
            if any(np.isnan((np.unique(wt.data)))):
                print(f, ev)
                problems.append(f+" "+ev+" dwtyr nan")
                problem_stands.append(f)
            if type(ncf['balance']['N']['to_water'][dscen[ev]].mask)==np.ndarray:
                print(f, ev," loytyy mask info, voi kayttaa indeksointia tulosten koonnissa")
                #problems.append(f+" "+ev+" loytyy mask info, voi kayttaa indeksointia ")
            else:
                print(f,ev,  "towater ei voikayttaa indeksointia tulosten koonnissa")
                #problems.append(f+" "+ev+" towater ei voi kayttaa indeksointia ")
            if type(ncf['esom']['Mass']['P1'][dscen[ev]].mask)==np.ndarray:
                print(f, "loytyy mask info, voi kayttaa indeksointia tulosten koonnissa")
                #problems.append(f+" "+ev+" loytyy mask info, voi kayttaa indeksointia ")
            elif ~ncf['esom']['Mass']['P1'][dscen[ev]].mask:
                print(f, "esom P1 mask  is false")
                #problems.append(f+" "+ev+" esom P1 mask is false, ei indeksointi")
            p1=ncf['esom']['Mass']['P1'][dscen[ev]]
            if any(np.isnan((np.unique(p1.data)))):
                print(f, ev)
                problems.append(f+" "+ev+" esom p1 nan")
                problem_stands.append(f)
            bal1=ncf['balance']['N']['to_water'][dscen[ev]]
            if any(np.isnan((np.unique(bal1.data)))):
                print(f, ev)
                problems.append(f+" "+ev+" balance N nan")
                problem_stands.append(f)
        print(f, "loppuu", stands.index(f), "of ", len(stands))
    problem_stands=[*set(problem_stands)]
    problem_stands.sort()
    with open('/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'+scenario+'_res_prob.txt', "a") as myfile:
        for item in problems:
            myfile.write("%s\n" % item)
    return problem_stands 

#scenarios=['Kuonanjoki_BAU_A', 'Kuonanjoki_BAU_B', 'Kuonanjoki_BIO_A', 'Kuonanjoki_BIO_B', 'Kuonanjoki_HII_A', 'Kuonanjoki_HII_B']

#k_baua_prs=checkres(scenarios[0])
#k_baub_prs=checkres(scenarios[1])
#k_bioa_prs=checkres(scenarios[2])
#k_biob_prs=checkres(scenarios[3])
#k_hiia_prs=checkres(scenarios[4])
#k_hiib_prs=checkres(scenarios[5])

scenarios=['Sorvasranta_BAU_A', 'Sorvasranta_BAU_B', 'Sorvasranta_BIO_A', 'Sorvasranta_BIO_B', 'Sorvasranta_HII_A', 'Sorvasranta_HII_B']

s_baua_prs=checkres(scenarios[0])
s_baub_prs=checkres(scenarios[1])
s_bioa_prs=checkres(scenarios[2])
s_biob_prs=checkres(scenarios[3])
s_hiia_prs=checkres(scenarios[4])
s_hiib_prs=checkres(scenarios[5])

s_all_pr=s_baua_prs+s_baub_prs+s_bioa_prs+s_biob_prs+s_hiia_prs+s_hiib_prs
s_all_pr = [*set(s_all_pr)]
s_all_pr.sort()

#scenarios=['Halvanjoki_BAU_A', 'Halvanjoki_BAU_B', 'Halvanjoki_BIO_A', 'Halvanjoki_BIO_B', 'Halvanjoki_HII_A', 'Halvanjoki_HII_B']

#h_baua_prs=checkres(scenarios[0])
#h_baub_prs=checkres(scenarios[1])
#h_bioa_prs=checkres(scenarios[2])
#h_biob_prs=checkres(scenarios[3])
#h_hiia_prs=checkres(scenarios[4])
#h_hiib_prs=checkres(scenarios[5])


def checkres_sorvas(scenario):
    resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'
    files=glob(resfol+scenario+'/*.nc')
    stands=[i[(len(resfol)+len(scenario)+1):-6] for i in files]
    stands = [*set(stands)]
    #stands.remove('43069586')
    stands.sort()
    problems=[]
    for f in stands:
        print(f, "alkaa")
        dscen=get_ev_dscens(f, scenario)
        events=[i[(len(resfol)+len(scenario)+1+len(f)+1):-3] for i in files if i.startswith(resfol+scenario+'/'+f)]
        events.sort()
        for ev in events:
            ff=resfol+scenario+"/"+f+"_"+ev+".nc"
            ncf=Dataset(ff, mode='r')
            wt = ncf['strip']['dwtyr'][dscen[ev],:, :]
            if any(np.isnan((np.unique(wt.data)))):
                print(f, ev)
                problems.append(f+" "+ev+" dwtyr nan")
            if type(ncf['balance']['N']['to_water'][dscen[ev]].mask)==np.ndarray:
                print(f, "loytyy mask info, voi kayttaa indeksointia tulosten koonnissa")
                #problems.append(f+" "+ev+" loytyy mask info, voi kayttaa indeksointia ")
            else:
                print(f, "towater ei voikayttaa indeksointia tulosten koonnissa")
                #problems.append(f+" "+ev+" towater ei voi kayttaa indeksointia ")
            if type(ncf['esom']['Mass']['P1'][dscen[ev]].mask)==np.ndarray:
                print(f, "loytyy mask info, voi kayttaa indeksointia tulosten koonnissa")
                #problems.append(f+" "+ev+" loytyy mask info, voi kayttaa indeksointia ")
            elif ~ncf['esom']['Mass']['P1'][dscen[ev]].mask:
                print(f, "esom P1 mask  is false")
                #problems.append(f+" "+ev+" esom P1 mask is false, ei indeksointi")
            p1=ncf['esom']['Mass']['P1'][dscen[ev]]
            if any(np.isnan((np.unique(p1.data)))):
                print(f, ev)
                problems.append(f+" "+ev+" esom p1 nan")
            bal1=ncf['balance']['N']['to_water'][dscen[ev]]
            if any(np.isnan((np.unique(bal1.data)))):
                print(f, ev)
                problems.append(f+" "+ev+" balance N nan")
        print(f, "loppuu", stands.index(f), "of ", len(stands))
    with open('/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'+scenario+'resok.txt', "a") as myfile:
        for item in problems:
            myfile.write("%s\n" % item)        



#puuston nykyinen hiilivarasto
##motista puusto kaikille kuvioille biomassojen mukaan
##tuloksesta eka rivi, skaalaus alku tilavuuden mukaan?
#tree_c0_data = pd.read_excel(r'/scratch/project_2002470/SUSI_HIILIPOLKU/Sorvasranta_BAU_A/motti/Sorvasranta_BAU_A.xlsx', sheet_name='puustot') 
#tree_c0_data.columns()
#tree_c0_data[bm_sum
#['Runkopuukg','Hukkapuukg','Oksatkg','Neulasetkg','Kannot_juuretkg','h_juuretkg','_2mm_juuretkg']
#maaperan hiilivarasto
##mineraalimaille YASSO

##turvemaille SUSI

#Turvekerroksen massat nc-tiedostosta. Kaivelin figures.py:ta tarvittavat koodit, joista tein def peat_ini_end(ff, scen), jolla saa maaperan alkutilanteen (inipeat+inimor). 
#Kuvassa nuo oli nimetty ”organic soil mass, kg m-2”, joten oletan, et luvut  kertoa 0.5:lla, et saadaan hiili.
#Simuloinneissa nyt oletuksena 2.5 m turvekerros.
#Periaatteessa tuosta saa myos lopputilanteen. Jonkinlainen ongelma on kuitenkin se, et hakkuun tullessa simulointi alkaa uudelleen ja turvekerroskin palautuu takaisin entiseen massaansa… Yhden kuvion lopputilanteen voisi kuitenkin laskea jotenkin nain: #

#endpeat = endpeat_n0 - (inipeat_n1 – endpeat_n1) - (inipeat_n2 – endpeat_n2) - … 
#endmor = endmor_n0 - (inimor_n1 – endmor_n1) - (inimor_n2 – endmor_n2) - …

#n viittaa tuossa hakkuiden paloittelemiin simulointeihin.#

#Parempi ratkaisu olisi tietysti jatkaa simulointia aina sellaisesta turvekerroksesta, joka edellisessa n:ssa jai jaljelle. 
#Tai toteuttaa hakkuiden simulointi niin, ettei simulointia tarvitse paloitella. 
#Kumpikin  tuntuu sen verran tyolailta ratkaisuilta


#scen=scenarios[0]

    

#start_yr_ini=1991
##start_date=datetime.datetime(191,1,1)
#totSimYears=30
#end_yr=start_yr_ini + totSimYears 
##end_date=datetime.datetime(2020,12,31)
##simyears=20
##sourcefile=r'/scratch/project_2002470/SUSI_HIILIPOLKU/inputs/sorvasranta_weather.csv'
##forc=make_weather(0, start_date, end_date, 20, sourcefile=sourcefile)

##sourcefile=r'/scratch/project_2002470/SUSI_HIILIPOLKU/inputs/kuonanjoki_weather.csv'
##forc=make_weather(0, start_date, end_date, 20, sourcefile=sourcefile)

##sourcefile=r'/scratch/project_2002470/SUSI_HIILIPOLKU/inputs/halvansuo_weather.csv'
##forc=make_weather(0, start_date, end_date, 20, sourcefile=sourcefile)



#indices = np.where(flag & ~flag.mask)
#print(data[indices])


#baua_sto=get_hiili(scenarios[0]); baub_sto=get_hiili(scenarios[1]); bioa_sto=get_hiili(scenarios[2]); biob_sto=get_hiili(scenarios[3]); hiia_sto=get_hiili(scenarios[4]); hiib_sto=get_hiili(scenarios[5]);



scenarios=['Sorvasranta_BAU_A', 'Sorvasranta_BAU_B', 'Sorvasranta_BIO_A', 'Sorvasranta_BIO_B', 'Sorvasranta_HII_A', 'Sorvasranta_HII_B']


#standareas_baua=get_areas(scenarios[0])
#standareas_baub=get_areas(scenarios[1]);
#standareas_bioa=get_areas(scenarios[2]);
#standareas_biob=get_areas(scenarios[3]);
#standareas_hiia=get_areas(scenarios[4]);
#standareas_hiib=get_areas(scenarios[5])

#standareas_baua=get_areas_prex(scenarios[0],s_baua_prs);
#standareas_baub=get_areas_prex(scenarios[1],s_baub_prs);
#standareas_bioa=get_areas_prex(scenarios[2],s_bioa_prs);
#standareas_biob=get_areas_prex(scenarios[3],s_biob_prs);
#standareas_hiia=get_areas_prex(scenarios[4],s_hiia_prs);
#standareas_hiib=get_areas_prex(scenarios[5],s_hiib_prs);

standareas_baua=get_areas_prex(scenarios[0],s_all_pr);
standareas_baub=get_areas_prex(scenarios[1],s_all_pr);
standareas_bioa=get_areas_prex(scenarios[2],s_all_pr);
standareas_biob=get_areas_prex(scenarios[3],s_all_pr);
standareas_hiia=get_areas_prex(scenarios[4],s_all_pr);
standareas_hiib=get_areas_prex(scenarios[5],s_all_pr);



#baua_sto_scaled=get_hiili_scaled(scenarios[0]); 
#baub_sto_scaled=get_hiili_scaled(scenarios[1]); 
#bioa_sto_scaled=get_hiili_scaled(scenarios[2]);
#biob_sto_scaled=get_hiili_scaled(scenarios[3]); 
#hiia_sto_scaled=get_hiili_scaled(scenarios[4]); 
#hiib_sto_scaled=get_hiili_scaled(scenarios[5]); 

#baua_sto_scaled=get_hiili_scaled_prex(scenarios[0],s_baua_prs); 
#baub_sto_scaled=get_hiili_scaled_prex(scenarios[1],s_baua_prs); 
#bioa_sto_scaled=get_hiili_scaled_prex(scenarios[2],s_baua_prs);
#biob_sto_scaled=get_hiili_scaled_prex(scenarios[3],s_baua_prs); 
#hiia_sto_scaled=get_hiili_scaled_prex(scenarios[4],s_baua_prs); 
#hiib_sto_scaled=get_hiili_scaled_prex(scenarios[5],s_baua_prs); 

baua_sto_scaled=get_hiili_scaled_prex(scenarios[0],s_all_pr); 
baub_sto_scaled=get_hiili_scaled_prex(scenarios[1],s_all_pr); 
bioa_sto_scaled=get_hiili_scaled_prex(scenarios[2],s_all_pr);
biob_sto_scaled=get_hiili_scaled_prex(scenarios[3],s_all_pr); 
hiia_sto_scaled=get_hiili_scaled_prex(scenarios[4],s_all_pr); 
hiib_sto_scaled=get_hiili_scaled_prex(scenarios[5],s_all_pr); 

#standareas_baua=get_areas_prex(scenarios[0],s_baua_prs)
#standareas_baub=get_areas_sorvas(scenarios[1]);
#standareas_bioa=get_areas_sorvas(scenarios[2]);
#standareas_biob=get_areas_sorvas(scenarios[3]);
#standareas_hiia=get_areas_sorvas(scenarios[4]);
#standareas_hiib=get_areas_sorvas(scenarios[5])

#baua_sto_scaled=get_hiili_scaled_exp_sorvas(scenarios[0]); 
#baub_sto_scaled=get_hiili_scaled_exp_sorvas(scenarios[1]); 
#bioa_sto_scaled=get_hiili_scaled_exp_sorvas(scenarios[2]);
#biob_sto_scaled=get_hiili_scaled_exp_sorvas(scenarios[3]); 
#hiia_sto_scaled=get_hiili_scaled_exp_sorvas(scenarios[4]); 
#hiib_sto_scaled=get_hiili_scaled_exp_sorvas(scenarios[5]); 




baua_c_sto_peat_scaled=[]
for i in range(len(baua_sto_scaled[0])):
    baua_c_sto_peat_scaled.append(baua_sto_scaled[0][i]*standareas_baua[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baua_co2eq=(sum(baua_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baua)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

baub_c_sto_peat_scaled=[]
for i in range(len(baub_sto_scaled[0])):
    baub_c_sto_peat_scaled.append(baub_sto_scaled[0][i]*standareas_baub[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baub_co2eq=(sum(baub_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baub)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

bioa_c_sto_peat_scaled=[]
for i in range(len(bioa_sto_scaled[0])):
    bioa_c_sto_peat_scaled.append(bioa_sto_scaled[0][i]*standareas_bioa[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

bioa_co2eq=(sum(bioa_c_sto_peat_scaled)*0.5/1000/(sum(standareas_bioa)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

biob_c_sto_peat_scaled=[]
for i in range(len(biob_sto_scaled[0])):
    biob_c_sto_peat_scaled.append(biob_sto_scaled[0][i]*standareas_biob[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

biob_co2eq=(sum(biob_c_sto_peat_scaled)*0.5/1000/(sum(standareas_biob)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiia_c_sto_peat_scaled=[]
for i in range(len(hiia_sto_scaled[0])):
    hiia_c_sto_peat_scaled.append(hiia_sto_scaled[0][i]*standareas_hiia[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiia_co2eq=(sum(hiia_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiia)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiib_c_sto_peat_scaled=[]
for i in range(len(hiib_sto_scaled[0])):
    hiib_c_sto_peat_scaled.append(hiib_sto_scaled[0][i]*standareas_hiib[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiib_co2eq=(sum(hiib_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiib)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

print(scenarios[0], "peat co2eq", baua_co2eq)
print(scenarios[1], "peat co2eq", baub_co2eq)
print(scenarios[2], "peat co2eq", bioa_co2eq)
print(scenarios[3], "peat co2eq", biob_co2eq)
print(scenarios[4], "peat co2eq", hiia_co2eq)
print(scenarios[5], "peat co2eq", hiib_co2eq)

baua_c_sto_peat_scaled=[] #now mor
for i in range(len(baua_sto_scaled[1])):
    baua_c_sto_peat_scaled.append(baua_sto_scaled[1][i]*standareas_baua[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baua_co2eq=(sum(baua_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baua)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

baub_c_sto_peat_scaled=[]
for i in range(len(baub_sto_scaled[1])):
    baub_c_sto_peat_scaled.append(baub_sto_scaled[1][i]*standareas_baub[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baub_co2eq=(sum(baub_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baub)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

bioa_c_sto_peat_scaled=[]
for i in range(len(bioa_sto_scaled[1])):
    bioa_c_sto_peat_scaled.append(bioa_sto_scaled[1][i]*standareas_bioa[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

bioa_co2eq=(sum(bioa_c_sto_peat_scaled)*0.5/1000/(sum(standareas_bioa)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

biob_c_sto_peat_scaled=[]
for i in range(len(biob_sto_scaled[1])):
    biob_c_sto_peat_scaled.append(biob_sto_scaled[1][i]*standareas_biob[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

biob_co2eq=(sum(biob_c_sto_peat_scaled)*0.5/1000/(sum(standareas_biob)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiia_c_sto_peat_scaled=[]
for i in range(len(hiia_sto_scaled[1])):
    hiia_c_sto_peat_scaled.append(hiia_sto_scaled[1][i]*standareas_hiia[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiia_co2eq=(sum(hiia_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiia)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiib_c_sto_peat_scaled=[]
for i in range(len(hiib_sto_scaled[1])):
    hiib_c_sto_peat_scaled.append(hiib_sto_scaled[1][i]*standareas_hiib[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiib_co2eq=(sum(hiib_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiib)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

print(scenarios[0],"mor co2eq", baua_co2eq)
print(scenarios[1],"mor co2eq", baub_co2eq)
print(scenarios[2],"mor co2eq", bioa_co2eq)
print(scenarios[3],"mor co2eq", biob_co2eq)
print(scenarios[4],"mor co2eq", hiia_co2eq)
print(scenarios[5],"mor co2eq", hiib_co2eq)


#baua=get_NP_sorvas(scenarios[0]); 
#baub=get_NP_sorvas(scenarios[1]); 
#bioa=get_NP_sorvas(scenarios[2]);
#biob=get_NP_sorvas(scenarios[3]); 
#hiia=get_NP_sorvas(scenarios[4]);
#hiib=get_NP_sorvas(scenarios[5]);
#baua=get_NP(scenarios[0]); 
#baub=get_NP(scenarios[1]); 
#bioa=get_NP(scenarios[2]);
#biob=get_NP(scenarios[3]); 
#hiia=get_NP(scenarios[4]);
#hiib=get_NP(scenarios[5]);

#baua=get_NP_prex(scenarios[0],s_baua_prs); 
baua=get_NP_prex(scenarios[0],s_all_pr); 
baub=get_NP_prex(scenarios[1],s_all_pr); 
bioa=get_NP_prex(scenarios[2],s_all_pr);
biob=get_NP_prex(scenarios[3],s_all_pr); 
hiia=get_NP_prex(scenarios[4],s_all_pr);
hiib=get_NP_prex(scenarios[5],s_all_pr);

print(scenarios[0], "nsum", sum(baua['nsum']))
print(scenarios[1], "nsum", sum(baub['nsum']))
print(scenarios[2], "nsum", sum(bioa['nsum']))
print(scenarios[3], "nsum", sum(biob['nsum']))
print(scenarios[4], "nsum", sum(hiia['nsum']))
print(scenarios[5], "nsum", sum(hiib['nsum']))

print(scenarios[0], "psum", sum(baua['psum']))
print(scenarios[1], "psum", sum(baub['psum']))
print(scenarios[2], "psum", sum(bioa['psum']))
print(scenarios[3], "psum", sum(biob['psum']))
print(scenarios[4], "psum", sum(hiia['psum']))
print(scenarios[5], "psum", sum(hiib['psum']))


print(scenarios[0],"n kg/ha/yr", sum(baua['nsum'])/sum(standareas_baua)/50)
print(scenarios[1],"n kg/ha/yr", sum(baub['nsum'])/sum(standareas_baub)/50)
print(scenarios[2],"n kg/ha/yr", sum(bioa['nsum'])/sum(standareas_bioa)/50)
print(scenarios[3],"n kg/ha/yr", sum(biob['nsum'])/sum(standareas_biob)/50)
print(scenarios[4],"n kg/ha/yr", sum(hiia['nsum'])/sum(standareas_hiia)/50)
print(scenarios[5],"n kg/ha/yr", sum(hiib['nsum'])/sum(standareas_hiib)/50)

print(scenarios[0],"p kg/ha/yr", sum(baua['psum'])/sum(standareas_baua)/50)
print(scenarios[1],"p kg/ha/yr", sum(baub['psum'])/sum(standareas_baub)/50)
print(scenarios[2],"p kg/ha/yr", sum(bioa['psum'])/sum(standareas_bioa)/50)
print(scenarios[3],"p kg/ha/yr", sum(biob['psum'])/sum(standareas_biob)/50)
print(scenarios[4],"p kg/ha/yr", sum(hiia['psum'])/sum(standareas_hiia)/50)
print(scenarios[5],"p kg/ha/yr", sum(hiib['psum'])/sum(standareas_hiib)/50)



#print("check res", scenarios[0])
#checkres_sorvas(scenarios[0])
#print("check res", scenarios[1])
#checkres_sorvas(scenarios[1])
#print("check res", scenarios[2])
#checkres_sorvas(scenarios[2])
#print("check res", scenarios[3])
#checkres_sorvas(scenarios[3])
#print("check res", scenarios[4])
#checkres_sorvas(scenarios[4])
#print("check res", scenarios[5])
#checkres_sorvas(scenarios[5])

"""
print("kuonanjoki")
scenarios=['Kuonanjoki_BAU_A', 'Kuonanjoki_BAU_B', 'Kuonanjoki_BIO_A', 'Kuonanjoki_BIO_B', 'Kuonanjoki_HII_A', 'Kuonanjoki_HII_B']

standareas_baua=get_areas(scenarios[0])
standareas_baub=get_areas(scenarios[1]);standareas_bioa=get_areas(scenarios[2]);standareas_biob=get_areas(scenarios[3]);standareas_hiia=get_areas(scenarios[4]);standareas_hiib=get_areas(scenarios[5])

#standareas_baua=get_areas_exp(scenarios[0])
#standareas_baub=get_areas_exp(scenarios[1]);standareas_bioa=get_areas_exp(scenarios[2]);standareas_biob=get_areas_exp(scenarios[3]);standareas_hiia=get_areas_exp(scenarios[4]);standareas_hiib=get_areas_exp(scenarios[5])

print(scenarios)

baua_sto_scaled=get_hiili_scaled(scenarios[0]); 
#43069586 kuvio poistettava ku probleemei
#baua_sto_scaled=get_hiili_scaled_exp(scenarios[0]); 
#baub_sto_scaled=get_hiili_scaled_exp(scenarios[1]); 
#bioa_sto_scaled=get_hiili_scaled_exp(scenarios[2]); 
#biob_sto_scaled=get_hiili_scaled_exp(scenarios[3]); 
#hiia_sto_scaled=get_hiili_scaled_exp(scenarios[4]); 
#hiib_sto_scaled=get_hiili_scaled_exp(scenarios[5]); 

baub_sto_scaled=get_hiili_scaled(scenarios[1]); 
bioa_sto_scaled=get_hiili_scaled(scenarios[2]); 
biob_sto_scaled=get_hiili_scaled(scenarios[3]); 
hiia_sto_scaled=get_hiili_scaled(scenarios[4]); 
hiib_sto_scaled=get_hiili_scaled(scenarios[5]); 


baua_c_sto_peat_scaled=[]
for i in range(len(baua_sto_scaled[0])):
    baua_c_sto_peat_scaled.append(baua_sto_scaled[0][i]*standareas_baua[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baua_co2eq=(sum(baua_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baua)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

baub_c_sto_peat_scaled=[]
for i in range(len(baub_sto_scaled[0])):
    baub_c_sto_peat_scaled.append(baub_sto_scaled[0][i]*standareas_baub[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baub_co2eq=(sum(baub_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baub)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

bioa_c_sto_peat_scaled=[]
for i in range(len(bioa_sto_scaled[0])):
    bioa_c_sto_peat_scaled.append(bioa_sto_scaled[0][i]*standareas_bioa[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

bioa_co2eq=(sum(bioa_c_sto_peat_scaled)*0.5/1000/(sum(standareas_bioa)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

biob_c_sto_peat_scaled=[]
for i in range(len(biob_sto_scaled[0])):
    biob_c_sto_peat_scaled.append(biob_sto_scaled[0][i]*standareas_biob[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

biob_co2eq=(sum(biob_c_sto_peat_scaled)*0.5/1000/(sum(standareas_biob)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiia_c_sto_peat_scaled=[]
for i in range(len(hiia_sto_scaled[0])):
    hiia_c_sto_peat_scaled.append(hiia_sto_scaled[0][i]*standareas_hiia[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiia_co2eq=(sum(hiia_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiia)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiib_c_sto_peat_scaled=[]
for i in range(len(hiib_sto_scaled[0])):
    hiib_c_sto_peat_scaled.append(hiib_sto_scaled[0][i]*standareas_hiib[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiib_co2eq=(sum(hiib_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiib)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

print(scenarios[0], "peat co2eq", baua_co2eq)
print(scenarios[1], "peat co2eq", baub_co2eq)
print(scenarios[2], "peat co2eq", bioa_co2eq)
print(scenarios[3], "peat co2eq", biob_co2eq)
print(scenarios[4], "peat co2eq", hiia_co2eq)
print(scenarios[5], "peat co2eq", hiib_co2eq)

baua_c_sto_peat_scaled=[] #now mor
for i in range(len(baua_sto_scaled[1])):
    baua_c_sto_peat_scaled.append(baua_sto_scaled[1][i]*standareas_baua[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baua_co2eq=(sum(baua_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baua)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

baub_c_sto_peat_scaled=[]
for i in range(len(baub_sto_scaled[1])):
    baub_c_sto_peat_scaled.append(baub_sto_scaled[1][i]*standareas_baub[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baub_co2eq=(sum(baub_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baub)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

bioa_c_sto_peat_scaled=[]
for i in range(len(bioa_sto_scaled[1])):
    bioa_c_sto_peat_scaled.append(bioa_sto_scaled[1][i]*standareas_bioa[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

bioa_co2eq=(sum(bioa_c_sto_peat_scaled)*0.5/1000/(sum(standareas_bioa)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

biob_c_sto_peat_scaled=[]
for i in range(len(biob_sto_scaled[1])):
    biob_c_sto_peat_scaled.append(biob_sto_scaled[1][i]*standareas_biob[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

biob_co2eq=(sum(biob_c_sto_peat_scaled)*0.5/1000/(sum(standareas_biob)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiia_c_sto_peat_scaled=[]
for i in range(len(hiia_sto_scaled[1])):
    hiia_c_sto_peat_scaled.append(hiia_sto_scaled[1][i]*standareas_hiia[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiia_co2eq=(sum(hiia_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiia)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiib_c_sto_peat_scaled=[]
for i in range(len(hiib_sto_scaled[1])):
    hiib_c_sto_peat_scaled.append(hiib_sto_scaled[1][i]*standareas_hiib[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiib_co2eq=(sum(hiib_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiib)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

print(scenarios[0],"mor co2eq", baua_co2eq)
print(scenarios[1],"mor co2eq", baub_co2eq)
print(scenarios[2],"mor co2eq", bioa_co2eq)
print(scenarios[3],"mor co2eq", biob_co2eq)
print(scenarios[4],"mor co2eq", hiia_co2eq)
print(scenarios[5],"mor co2eq", hiib_co2eq)




baua=get_NP_exp(scenarios[0]); 
baub=get_NP_exp(scenarios[1]); 
bioa=get_NP_exp(scenarios[2]);
biob=get_NP_exp(scenarios[3]); 
hiia=get_NP_exp(scenarios[4]);
hiib=get_NP_exp(scenarios[5]);

print(scenarios[0], "nsum", sum(baua['nsum']))
print(scenarios[1], "nsum", sum(baub['nsum']))
print(scenarios[2], "nsum", sum(bioa['nsum']))
print(scenarios[3], "nsum", sum(biob['nsum']))
print(scenarios[4], "nsum", sum(hiia['nsum']))
print(scenarios[5], "nsum", sum(hiib['nsum']))

print(scenarios[0], "psum", sum(baua['psum']))
print(scenarios[1], "psum", sum(baub['psum']))
print(scenarios[2], "psum", sum(bioa['psum']))
print(scenarios[3], "psum", sum(biob['psum']))
print(scenarios[4], "psum", sum(hiia['psum']))
print(scenarios[5], "psum", sum(hiib['psum']))


print(scenarios[0],"n kg/ha/yr", sum(baua['nsum'])/sum(standareas_baua)/50)
print(scenarios[1],"n kg/ha/yr", sum(baub['nsum'])/sum(standareas_baub)/50)
print(scenarios[2],"n kg/ha/yr", sum(bioa['nsum'])/sum(standareas_bioa)/50)
print(scenarios[3],"n kg/ha/yr", sum(biob['nsum'])/sum(standareas_biob)/50)
print(scenarios[4],"n kg/ha/yr", sum(hiia['nsum'])/sum(standareas_hiia)/50)
print(scenarios[5],"n kg/ha/yr", sum(hiib['nsum'])/sum(standareas_hiib)/50)

print(scenarios[0],"p kg/ha/yr", sum(baua['psum'])/sum(standareas_baua)/50)
print(scenarios[1],"p kg/ha/yr", sum(baub['psum'])/sum(standareas_baub)/50)
print(scenarios[2],"p kg/ha/yr", sum(bioa['psum'])/sum(standareas_bioa)/50)
print(scenarios[3],"p kg/ha/yr", sum(biob['psum'])/sum(standareas_biob)/50)
print(scenarios[4],"p kg/ha/yr", sum(hiia['psum'])/sum(standareas_hiia)/50)
print(scenarios[5],"p kg/ha/yr", sum(hiib['psum'])/sum(standareas_hiib)/50)

print("Halvanjoki")

scenarios=['Halvanjoki_BAU_A', 'Halvanjoki_BAU_B', 'Halvanjoki_BIO_A', 'Halvanjoki_BIO_B', 'Halvanjoki_HII_A', 'Halvanjoki_HII_B']
#standareas_baua=get_areas_exp2(scenarios[0])
#standareas_baub=get_areas_exp2(scenarios[1]);
#standareas_bioa=get_areas(scenarios[2]);
#standareas_biob=get_areas_exp2(scenarios[3]);
#standareas_hiia=get_areas_exp2(scenarios[4]);
#standareas_hiib=get_areas_exp2(scenarios[5])

standareas_baua=get_areas(scenarios[0]);
standareas_baub=get_areas(scenarios[1]);
standareas_bioa=get_areas(scenarios[2]);
standareas_biob=get_areas(scenarios[3]);
standareas_hiia=get_areas(scenarios[4]);
standareas_hiib=get_areas(scenarios[5]);


print(scenarios)

#baua_sto_scaled=get_hiili_scaled_exp2(scenarios[0]); 
#baub_sto_scaled=get_hiili_scaled_exp2(scenarios[1]); 
#bioa_sto_scaled=get_hiili_scaled(scenarios[2]); 
#biob_sto_scaled=get_hiili_scaled_exp2(scenarios[3]); 
#hiia_sto_scaled=get_hiili_scaled_exp2(scenarios[4]); 
#hiib_sto_scaled=get_hiili_scaled_exp2(scenarios[5]); 
baua_sto_scaled=get_hiili_scaled(scenarios[0]); 
baub_sto_scaled=get_hiili_scaled(scenarios[1]); 
bioa_sto_scaled=get_hiili_scaled(scenarios[2]); 
biob_sto_scaled=get_hiili_scaled(scenarios[3]); 
hiia_sto_scaled=get_hiili_scaled(scenarios[4]); 
hiib_sto_scaled=get_hiili_scaled(scenarios[5]); 

baua_c_sto_peat_scaled=[]
for i in range(len(baua_sto_scaled[0])):
    baua_c_sto_peat_scaled.append(baua_sto_scaled[0][i]*standareas_baua[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baua_co2eq=(sum(baua_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baua)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

baub_c_sto_peat_scaled=[]
for i in range(len(baub_sto_scaled[0])):
    baub_c_sto_peat_scaled.append(baub_sto_scaled[0][i]*standareas_baub[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baub_co2eq=(sum(baub_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baub)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

bioa_c_sto_peat_scaled=[]
for i in range(len(bioa_sto_scaled[0])):
    bioa_c_sto_peat_scaled.append(bioa_sto_scaled[0][i]*standareas_bioa[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

bioa_co2eq=(sum(bioa_c_sto_peat_scaled)*0.5/1000/(sum(standareas_bioa)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

biob_c_sto_peat_scaled=[]
for i in range(len(biob_sto_scaled[0])):
    biob_c_sto_peat_scaled.append(biob_sto_scaled[0][i]*standareas_biob[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

biob_co2eq=(sum(biob_c_sto_peat_scaled)*0.5/1000/(sum(standareas_biob)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiia_c_sto_peat_scaled=[]
for i in range(len(hiia_sto_scaled[0])):
    hiia_c_sto_peat_scaled.append(hiia_sto_scaled[0][i]*standareas_hiia[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiia_co2eq=(sum(hiia_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiia)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiib_c_sto_peat_scaled=[]
for i in range(len(hiib_sto_scaled[0])):
    hiib_c_sto_peat_scaled.append(hiib_sto_scaled[0][i]*standareas_hiib[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiib_co2eq=(sum(hiib_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiib)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

print(scenarios[0], "peat co2eq", baua_co2eq)
print(scenarios[1], "peat co2eq", baub_co2eq)
print(scenarios[2], "peat co2eq", bioa_co2eq)
print(scenarios[3], "peat co2eq", biob_co2eq)
print(scenarios[4], "peat co2eq", hiia_co2eq)
print(scenarios[5], "peat co2eq", hiib_co2eq)

baua_c_sto_peat_scaled=[] #now mor
for i in range(len(baua_sto_scaled[1])):
    baua_c_sto_peat_scaled.append(baua_sto_scaled[1][i]*standareas_baua[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baua_co2eq=(sum(baua_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baua)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

baub_c_sto_peat_scaled=[]
for i in range(len(baub_sto_scaled[1])):
    baub_c_sto_peat_scaled.append(baub_sto_scaled[1][i]*standareas_baub[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

baub_co2eq=(sum(baub_c_sto_peat_scaled)*0.5/1000/(sum(standareas_baub)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

bioa_c_sto_peat_scaled=[]
for i in range(len(bioa_sto_scaled[1])):
    bioa_c_sto_peat_scaled.append(bioa_sto_scaled[1][i]*standareas_bioa[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

bioa_co2eq=(sum(bioa_c_sto_peat_scaled)*0.5/1000/(sum(standareas_bioa)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

biob_c_sto_peat_scaled=[]
for i in range(len(biob_sto_scaled[1])):
    biob_c_sto_peat_scaled.append(biob_sto_scaled[1][i]*standareas_biob[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

biob_co2eq=(sum(biob_c_sto_peat_scaled)*0.5/1000/(sum(standareas_biob)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiia_c_sto_peat_scaled=[]
for i in range(len(hiia_sto_scaled[1])):
    hiia_c_sto_peat_scaled.append(hiia_sto_scaled[1][i]*standareas_hiia[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiia_co2eq=(sum(hiia_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiia)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

hiib_c_sto_peat_scaled=[]
for i in range(len(hiib_sto_scaled[1])):
    hiib_c_sto_peat_scaled.append(hiib_sto_scaled[1][i]*standareas_hiib[i]*10000) #organic soil mass kg/m2, areas in m2, #kg

hiib_co2eq=(sum(hiib_c_sto_peat_scaled)*0.5/1000/(sum(standareas_hiib)))*(44/12) #t kg/ha to c (*05) ja sitten 44/12 to co2

print(scenarios[0],"mor co2eq", baua_co2eq)
print(scenarios[1],"mor co2eq", baub_co2eq)
print(scenarios[2],"mor co2eq", bioa_co2eq)
print(scenarios[3],"mor co2eq", biob_co2eq)
print(scenarios[4],"mor co2eq", hiia_co2eq)
print(scenarios[5],"mor co2eq", hiib_co2eq)


#baua=get_NP_exp2(scenarios[0]); 
#baub=get_NP_exp2(scenarios[1]); 
#bioa=get_NP(scenarios[2]);
##biob=get_NP_exp2(scenarios[3]); 
#hiia=get_NP_exp2(scenarios[4]);
#hiib=get_NP_exp2(scenarios[5]);
baua=get_NP(scenarios[0]);
baub=get_NP(scenarios[1]);
bioa=get_NP(scenarios[2]);
biob=get_NP(scenarios[3]);
hiia=get_NP(scenarios[4]);
hiib=get_NP(scenarios[5]);

print(scenarios[0], "nsum", sum(baua['nsum']))
print(scenarios[1], "nsum", sum(baub['nsum']))
print(scenarios[2], "nsum", sum(bioa['nsum']))
print(scenarios[3], "nsum", sum(biob['nsum']))
print(scenarios[4], "nsum", sum(hiia['nsum']))
print(scenarios[5], "nsum", sum(hiib['nsum']))

print(scenarios[0], "psum", sum(baua['psum']))
print(scenarios[1], "psum", sum(baub['psum']))
print(scenarios[2], "psum", sum(bioa['psum']))
print(scenarios[3], "psum", sum(biob['psum']))
print(scenarios[4], "psum", sum(hiia['psum']))
print(scenarios[5], "psum", sum(hiib['psum']))


print(scenarios[0],"n kg/ha/yr", sum(baua['nsum'])/sum(standareas_baua)/50)
print(scenarios[1],"n kg/ha/yr", sum(baub['nsum'])/sum(standareas_baub)/50)
print(scenarios[2],"n kg/ha/yr", sum(bioa['nsum'])/sum(standareas_bioa)/50)
print(scenarios[3],"n kg/ha/yr", sum(biob['nsum'])/sum(standareas_biob)/50)
print(scenarios[4],"n kg/ha/yr", sum(hiia['nsum'])/sum(standareas_hiia)/50)
print(scenarios[5],"n kg/ha/yr", sum(hiib['nsum'])/sum(standareas_hiib)/50)

print(scenarios[0],"p kg/ha/yr", sum(baua['psum'])/sum(standareas_baua)/50)
print(scenarios[1],"p kg/ha/yr", sum(baub['psum'])/sum(standareas_baub)/50)
print(scenarios[2],"p kg/ha/yr", sum(bioa['psum'])/sum(standareas_bioa)/50)
print(scenarios[3],"p kg/ha/yr", sum(biob['psum'])/sum(standareas_biob)/50)
print(scenarios[4],"p kg/ha/yr", sum(hiia['psum'])/sum(standareas_hiia)/50)
print(scenarios[5],"p kg/ha/yr", sum(hiib['psum'])/sum(standareas_hiib)/50)


#baua_d0_sto=get_hiili(scenarios[0],0); baua_d1_sto=get_hiili(scenarios[0],1); baua_d2_sto=get_hiili(scenarios[0],2);  baua_d3_sto=get_hiili(scenarios[0],3)
#baub_d0_sto=get_hiili(scenarios[1],0); baub_d1_sto=get_hiili(scenarios[1],1); baub_d2_sto=get_hiili(scenarios[1],2);  baub_d3_sto=get_hiili(scenarios[1],3)
#hiia_d0_sto=get_hiili(scenarios[4],0); hiia_d1_sto=get_hiili(scenarios[4],1); hiia_d2_sto=get_hiili(scenarios[4],2);  hiia_d3_sto=get_hiili(scenarios[4],3)
#hiib_d0_sto=get_hiili(scenarios[5],0); hiib_d1_sto=get_hiili(scenarios[5],1); hiib_d2_sto=get_hiili(scenarios[5],2);  hiib_d3_sto=get_hiili(scenarios[5],3)

#baua_c_sto_peat=[]
#for i in range(len(baua_sto[0])):
#    baua_c_sto_peat.append(baua_sto[0][i]*standareas_baua[i]*10000) #organic soil mass kg/m2, areas in m2

#sorvasc_stokgha=[]
#sorvasc_stokgha.append((sum(baua_c_sto_peat)*0.5)/1000000/(sum(standareas_baua))) #milj. kg/ha

#c_sto_peat=[]
#for i in range(len(baub_sto[0])):
#    c_sto_peat.append(baub_sto[0][i]*standareas_baub[i]*10000) #organic soil mass kg/m2, areas in m2

#sorvasc_stokgha.append((sum(c_sto_peat)*0.5)/1000000/(sum(standareas_baub))) #milj. kg/ha

#c_sto_peat=[]
#for i in range(len(hiia_sto[0])):
#    c_sto_peat.append(hiia_sto[0][i]*standareas_hiia[i]*10000) #organic soil mass kg/m2, areas in m2

#sorvasc_stokgha.append((sum(c_sto_peat)*0.5)/1000000/(sum(standareas_hiia))) #milj. kg/ha

#c_sto_peat=[]
#for i in range(len(hiib_sto[0])):
#    c_sto_peat.append(hiib_sto[0][i]*standareas_hiib[i]*10000) #organic soil mass kg/m2, areas in m2

#sorvasc_stokgha.append((sum(c_sto_peat)*0.5)/1000000/(sum(standareas_hiib))) #milj. kg/ha





#hiia_c_sto_peat=[]
#for i in range(len(hiia_d0_sto[0])):
#    hiia_c_sto_peat.append(hiia_d0_sto[0][i]*standareas_hiia[i]*10000) #sto in kg/m2 areas in m2

#(sum(hiia_c_sto_peat)*0.58)/1000000/(sum(standareas_hiia)) #milj. kg/ha

####################################################################################################
####################################################################################################

#get areas from motti for stand, ens multiply kg/ha7yr with ha to get results!


#np.mean(ncfs["n{0}".format(0)]['balance']['N']['to_water'][0][:,1:-1],axis=1)[1:]
        
#ff = r'path_to_file/file_name.nc'

#inipeat, inimor, endpeat, endmor = peat_ini_end(ff,0) # soil mass, kg/m2, peat = peat layers, mor = upper humus & litter layers

"""

#C -> CO2 muunnoskerroin:
#c_to_co2 = 44/12.0

