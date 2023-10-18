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
#cd /scratch/project_2002470/
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset  
import numpy as np
import pandas as pd

from SUSI_HIILIPOLKU.susi.figures import hydrology
from SUSI_HIILIPOLKU.susi.figures import stand
from SUSI_HIILIPOLKU.susi.figures import mass
#from SUSI_HIILIPOLKU.susi.figures import carbon
from SUSI_HIILIPOLKU.susi.figures import nutrient_balance
from SUSI_HIILIPOLKU.susi.figures import compare_1
from SUSI_HIILIPOLKU.susi.figures import compare_scens
from SUSI_HIILIPOLKU.susi.figures import create_profile_line
from SUSI_HIILIPOLKU.susi.figures import create_profile_boxplot
#from SUSI_HIILIPOLKU.outputs_hiilipolku_neuvontaan import carbon

def carbon(ff, scen,outfol,st,sc):
    import matplotlib.pylab as plt 
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    facecolor = '#f2f5eb'
    fs = 15
    fig = plt.figure(num='carbon', figsize=(15,18))   #width, height
    fig.suptitle('Carbon balance components', fontsize = fs+2)
    gs = gridspec.GridSpec(ncols=12, nrows=14, figure=fig, wspace=0.5, hspace=0.5)
    mass_to_c = 0.5
    
    #-------soil C balance kg/ha/yr---------------------------
    litter = (ncf['groundvegetation']['ds_litterfall'][scen,:, :]/10000.\
        + ncf['groundvegetation']['h_litterfall'][scen,:, :]/10000.\
        + ncf['groundvegetation']['s_litterfall'][scen,:, :]/10000.\
        + ncf['stand']['nonwoodylitter'][scen, :, :]/10000.\
        + ncf['stand']['woodylitter'][scen, :, :]/10000.)*mass_to_c
    
    soil = ncf['esom']['Mass']['out'][scen,:, :]/10000.*-1 * mass_to_c + litter
    
    out = ncf['esom']['Mass']['out'][scen,:, :]/10000.*-1 * mass_to_c
    
    wt = np.mean(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    cols = np.shape(wt)[0]
    sd = np.std(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    wtls = np.mean(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    sdls = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    wtmin = min(wtls) -0.2
    
    #-------water table as a reference-----------------------------
    ax = fig.add_subplot(gs[12:, :6])
    ax = create_profile_line(ax, wt, wtmin, sd, cols, 'WT m', 'annual', fs, facecolor, 'blue', hidex=True, hidey=False)
    
    #----------- water table in H------------------------------   
    elevation = np.array(ncf['strip']['elevation'])
    wtls = np.mean(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    sd = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    h = elevation + wtls
    
    ax = fig.add_subplot(gs[12:, 6:])
    ax = create_profile_line(ax, h, wtmin, sd, cols, 'WT m', 'annual', fs, facecolor, \
                             'blue', hidex=True, hidey=False, elevation = elevation)
    
    
    #-------------LMW to Ditch----------------------------
    lmwtoditch = ncf['balance']['C']['LMWdoc_to_water'][scen,:, :] *-1
    ax = fig.add_subplot(gs[10:12, :6])
    df = pd.DataFrame(data=lmwtoditch, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'brown', 'LMW to ditch', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-------------HMW to Ditch----------------------------
    hmwtoditch = ncf['balance']['C']['HMW_to_water'][scen,:, :]*-1
    ax = fig.add_subplot(gs[10:12, 6:])
    df = pd.DataFrame(data=hmwtoditch, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'brown', 'HMW to ditch', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------LMW to atmosphere--------------------
    lmwtoatm = ncf['balance']['C']['LMWdoc_to_atm'][scen,:, :]*-1
    ax = fig.add_subplot(gs[8:10, :6])
    df = pd.DataFrame(data=lmwtoatm, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'orange', 'LMW to atmosphere', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------HMW to atmosphere--------------------
    hmwtoatm = ncf['balance']['C']['HMW_to_atm'][scen,:, :]*-1
    ax = fig.add_subplot(gs[8:10, 6:])
    df = pd.DataFrame(data=hmwtoatm, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'orange', 'HMW to atmosphere', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------CO2C to atmosphere--------------------
    co2 = ncf['balance']['C']['co2c_release'][scen,:, :]*-1
    ax = fig.add_subplot(gs[6:8, :6])
    df = pd.DataFrame(data=co2, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'grey', '$CO_2C$ to atmosphere', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------CH4C to atmosphere--------------------
    co2 = ncf['balance']['C']['ch4c_release'][scen,:, :]*-1
    ax = fig.add_subplot(gs[6:8, 6:])
    df = pd.DataFrame(data=co2, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'grey', '$CH_4C$ to atmosphere', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------stand litter in--------------------
    standl = ncf['balance']['C']['stand_litter_in'][scen,:, :]
    ax = fig.add_subplot(gs[4:6, :6])
    df = pd.DataFrame(data=standl, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', 'Stand litter', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------ground vegetation litter in--------------------
    gvl = ncf['balance']['C']['gv_litter_in'][scen,:, :]
    ax = fig.add_subplot(gs[4:6, 6:])
    df = pd.DataFrame(data=gvl, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', 'Groundvegetation litter', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------soil balance c--------------------
    soilc = ncf['balance']['C']['soil_c_balance_c'][scen,:, :]
    ax = fig.add_subplot(gs[2:4, :6])
    df = pd.DataFrame(data=soilc, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'black', 'Soil C balance in C', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------soil balance co2 equivalents--------------------
    soilco2 = ncf['balance']['C']['soil_c_balance_co2eq'][scen,:, :]
    ax = fig.add_subplot(gs[2:4, 6:])
    df = pd.DataFrame(data=soilco2, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'black', 'Soil C balance in $CO_2$ eq', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------stand balance c--------------------
    standc = ncf['balance']['C']['stand_c_balance_c'][scen,:, :]
    ax = fig.add_subplot(gs[:2, :6])
    df = pd.DataFrame(data=standc, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', 'Stand C balance in C', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------stand balance co2 equivalents--------------------
    standco2 = ncf['balance']['C']['stand_c_balance_co2eq'][scen,:, :]
    ax = fig.add_subplot(gs[:2, 6:])
    df = pd.DataFrame(data=standco2, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', 'Stand C balance in $CO_2$ eq', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    fig.savefig(outfol+st+'_'+sc+'_carbon.png', dpi=300)
    ncf.close()
    return fig
    
#Leenalta kuva skripti
#stands=['24289396','24289397','24289702','24289703','24294089','24294099','24294100','24294101','24294105']
#stands=['24289396','24289397','24289702','24289703','24294089','24294099','24294100','24294101']
#allstands=['24294099','24289389','24289390','24289391','24289395','24289396','24289397','24289702','24289703','24294089','24294090','24294091','24294092','24294093','24294094','24294095','24294096','24294097','24294098','24294100','24294101','24294102','24294105','25732072','28755336','33614143']

#kolmella kuviolla ongelmia susi-tuloksissa 24294105 kuivatustilanne on 1 eli ei tuu susi-mallinnukseen (koska?) ja 24294099 BIOaskenaario puuttuu
#24294089 bioa kokonaan nan, samoin hiia kokonaan nan
#stands=['24294089','24294099']
stands=['24289396','24289397','24289702','24289703','24294100','24294101']
allstands=['24289389','24289390','24289391','24289395','24289396','24289397','24289702','24289703','24294090','24294091','24294092','24294093','24294094','24294095','24294096','24294097','24294098','24294100','24294101','24294102','25732072','28755336','33614143']
minstands=list(set(allstands) - set(stands))

scenarios=['Halvanjoki_BAU_A', 'Halvanjoki_BAU_B', 'Halvanjoki_BIO_A', 'Halvanjoki_BIO_B', 'Halvanjoki_HII_A', 'Halvanjoki_HII_B']
mottifol=r'/scratch/project_2002470/HIILIPOLKU_data/'
resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'

dscens_686=[]
for s in stands:
   for scen in scenarios:
      ldsce=get_ev_dscens(s,scen)
      dscens_686.append(ldsce)
#kaikki oli dscen=1, ei siis ojitusta


### kuvailuskiprti kuvio j aresepti kuvaus
def get_motti_data():
        
    kuviot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/Halvanjoki_BIO_A/motti/Halvanjoki_BIO_A_kuviot.csv', encoding='latin1', sep=';')
    puustot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/Halvanjoki_BIO_A/motti/Halvanjoki_BIO_A_puustot.csv', encoding='latin1', sep=';')
    poistumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/Halvanjoki_BIO_A/motti/Halvanjoki_BIO_A_poistumat.csv', encoding='latin1', sep=';')
    tapahtumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/Halvanjoki_BIO_A/motti/Halvanjoki_BIO_A_tapahtumat.csv', encoding='latin1', sep=';', index_col=False)
    tulot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/Halvanjoki_BIO_A/motti/Halvanjoki_BIO_A_kasvut.csv', encoding='latin1', sep=';')
    
    return kuviot, puustot, poistumat, tapahtumat, tulot


standis_hbioa=find_i(scenarios[2],['24290600'])

hbaua_kuviot, hbaua_puustot, hbaua_poistumat, hbaua_tapahtumat, hbaua_tulot= get_motti_data(scenarios[0])

hbaua_kuviot
st_tdata= hbaua_tapahtumat[hbaua_tapahtumat['kuvio']==int(stands[2])]


st_tdata['lann']>0  

standinfo=[]
standinfo.append("id "+stands[2]+", Ala: "+str(hbaua_kuviot[(hbaua_kuviot['kuvio']==int(stands[2]))]['ala'].values[0])+",  alaryhma "+str(hbaua_kuviot[(hbaua_kuviot['kuvio']==int(stands[2]))]['alaryh'].values[0])+", kuivatustilanne:"+str(hbaua_kuviot[(hbaua_kuviot['kuvio']==int(stands[2]))]['kuivatustilanne'].values[0])+
", kasvupaikka:"+str(hbaua_kuviot[(hbaua_kuviot['kuvio']==int(stands[2]))]['kasvup'].values[0])+", hakkuurajoite:"+str(hbaua_kuviot[(hbaua_kuviot['kuvio']==int(stands[2]))]['HakRajoite'].values[0])+", Puulaji:"+str(hbaua_kuviot[(hbaua_kuviot['kuvio']==int(stands[2]))]['puulaji'].values[0])+', paatehakkuu: vuonna'+str(st_tdata[st_tdata['_2005']>-1]['_2005'].values))


###


"""
tapahtumat = tapahtumat.rename(columns={'kuvio':'KUVIO',
                                  '_2003': '2003',
                                  '_2004': '2004',
                                  '_2005': '2005',
                                  '_2001': '2001'}) # 2001 lannoitus, #2002 kunnostusipjitus

"""

#scenario=scenarios[2]#
#sc=scenarios[2]
#files=glob(resfol+sc+'/*.nc')
#st=stands[0]
#events=[i[(len(resfol)+len(sc)+1+len(st)+1):-3] for i in files if i.startswith(resfol+sc+'/'+st)]
#events.sort()
#ev=events[0]
#for ev in events:
#    sto_values=peat_ini_end(resfol+scenario+"/"+f+"_"+ev+".nc", dscen[ev])
#    storages.append(sto_values)
check_resfiles(scenarios)

colors=['tab:blue', 'lightblue', 'tab:green','lightgreen','tab:orange','navajowhite']
labels = ['1: BAUA', '1.5:BAUB', '2:BIOA', '2.5:BIOB', '3:HIIA', '3.5:HIIB']
res=[]
indices_ls=[]
vol_ends=[]
vol_starts=[]
gr=[]
hakkuut=[]
cb_total_as_CO2=[]
cb_as_CO2=[]
nleachs1=[]
pleachs1=[]
nleachs=[]
pleachs=[]
#standco2=[]
#soilco2=[]
standco21=[]
soilco21=[]
hakkuut_ls=[]
kasvut_ls=[]
vols_inend=[]

#stands[7] in missing_stands_list[0]
########HUOM HAKKUUT EI VÄLTSII LASKETTU OIKEIN KUN HARVENNUKSET EI NÄY??? 24289703 kuviolla ainkaskin hiia skenaariossa harvennuksia, mutta ei tuu hakkuisiin
#st=stands[2]

for st in stands:
    for sc in scenarios:
        hakkuut_ls=[];kasvut_ls=[];indices_ls=[];vol_ends=[];vol_starts=[];nleachs1=0;pleachs1=0; standco21=0;soilco21=0;kuviot=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+sc[:-6]+'/'+sc+'/motti/'+sc+'_kuviot.csv', encoding='latin1', sep=';', index_col=False);kdata=kuviot[(kuviot['kuvio']==int(st))]; files=glob(resfol+sc+'/*.nc');events=[i[(len(resfol)+len(sc)+1+len(st)+1):-3] for i in files if i.startswith(resfol+sc+'/'+st)];events.sort();print(st,sc,events, 'alkutilavuus',kdata['alkutilavuus'].values[0]);hak_kertym=kdata['alkutilavuus'].values[0]; gr_kertym=0;print("hak_ker",hak_kertym,gr_kertym);ev_count=0;
        for ev in events:
            ev_count+=events.index(ev);print(ev,"starts", hak_kertym);ff = os.path.join(resfol,sc,st+'_'+ev+'.nc');ncf=Dataset(ff, mode='r');vol = ncf['stand']['volume'][1,:, :];vol00=np.mean(vol[0,:]);vol_starts.append(vol00);
            if np.ma.is_masked(ncf['stand']['volume'][1]):  
                indices=np.where(~ncf['stand']['volume'][1].mask);lastind=max(indices[0]);
            else:
                lastind=len(ncf['stand']['volume'][1])-1 
            indices_ls.append(lastind);
            if lastind==0:
                continue
            totvol = np.mean(vol[lastind,:]);vol_ends.append(totvol);
            if ev_count>0:
                hak_kertym+=vol_ends[ev_count-1]-vol_starts[ev_count]
                #hak_kertym+=-vol00+totvol;#pitaa korjata et laskee et laksee seuraavan eventon alku, ja eeltävän loppu, eli laskettava kun kaikki eventit käyty läpi
            else:
                hak_kertym=0
            gr_kertym+=+(totvol-vol00);NL0 = ncf['balance']['N']['to_water'][1,:, :];NL0_meana = np.mean(NL0[:,:]);PL0 = ncf['balance']['P']['to_water'][1,:, :];PL0_meana = np.mean(PL0[:,:]);nleachs1+=(NL0_meana*lastind);pleachs1+=(PL0_meana*lastind);standco2 = ncf['balance']['C']['stand_c_balance_co2eq'][1,:, :];standco2_meana = np.mean(standco2[:,:]);standco21+=(standco2_meana*lastind);soilco2 = ncf['balance']['C']['soil_c_balance_co2eq'][1,:, :];soilco2_meana = np.mean(soilco2[:,:]);soilco21+=(soilco2_meana*lastind); print("hak_kerty",hak_kertym)
        vols_inend.append(totvol);print(st,sc,"all events calc");hakkuut.append(hak_kertym);gr.append(gr_kertym);nleachs.append(nleachs1/sum(indices_ls));pleachs.append(pleachs1/sum(indices_ls));cb_total_as_CO2.append(standco21/sum(indices_ls));cb_as_CO2.append(soilco21/sum(indices_ls));print(st,sc,hakkuut, gr,nleachs,pleachs, cb_total_as_CO2,cb_as_CO2)



###
hakkuut.append(kdata['alkutilavuus'].values[0]-vol_starts[0])+(vol_ends[0]-vol_starts[1])+(vol_ends[1]-vol_starts[2])


ff='/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/Halvanjoki_BAU_A/24289396_n0.nc'
ff2='/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/Halvanjoki_BAU_A/24289396_n2.nc'
ff1='/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/Halvanjoki_BAU_A/24289396_n1.nc'
ncf=Dataset(ff, mode='r')                                        
ncf1=Dataset(ff1, mode='r')                                        
ncf2=Dataset(ff2, mode='r')                                        


kuviot=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+sc[:-6]+'/'+sc+'/motti/'+sc+'_kuviot.csv', encoding='latin1', sep=';', index_col=False)
kdata=kuviot[(kuviot['kuvio']==int(st))]


vol = ncf['stand']['volume'][1,:, :]
vol00=np.mean(vol[0,:])
indices=np.where(~ncf['stand']['volume'][1].mask)
lastind=max(indices[0])
#endpeat += ncf['esom']['Mass'][sto][dscen,lastind, :]/10000. #pitais ottaa vika jossa on arvoja 
totvol = np.mean(vol[lastind,:])

vol1 = ncf1['stand']['volume'][1,:, :]
vol10=np.mean(vol1[0,:])
indices1=np.where(~ncf1['stand']['volume'][1].mask)
lastind1=max(indices1[0])
totvol1 = np.mean(vol1[lastind1,:])

vol2 = ncf2['stand']['volume'][1,:, :]
vol20=np.mean(vol2[0,:])
indices2=np.where(~ncf2['stand']['volume'][1].mask)
lastind2=max(indices2[0])
totvol2 = np.mean(vol2[lastind2,:])

hakkuut=(kdata['alkutilavuus'].values[0]-vol00)+(totvol-vol10)+(totvol1-vol20)

v_end=totvol2

vg = ncf['stand']['volumegrowth'][scen,:, :]
vg00=np.mean(vg[0,:])
indices=np.where(~ncf['stand']['volumegrowth'][1].mask)
lastind=max(indices[0])
vg0end = np.mean(vg[:,:]) * lastind

totgr= totvol-(kdata['alkutilavuus'].values[0]-vol00)+(totvol1-vol10)+(totvol2-vol20)

#vg1 = ncf1['stand']['volumegrowth'][scen,:, :]
#vg2 = ncf2['stand']['volumegrowth'][scen,:, :]

NL0 = ncf['balance']['N']['to_water'][1,:, :]
NL0_meana = np.mean(NL0[:,:]) 
NL1 = ncf1['balance']['N']['to_water'][1,:, :]
NL1_meana = np.mean(NL1[:,:]) 
NL2 = ncf2['balance']['N']['to_water'][1,:, :]
NL2_meana = np.mean(NL2[:,:]) 

PL0 = ncf['balance']['P']['to_water'][1,:, :]
PL0_meana = np.mean(PL0[:,:]) 
PL1 = ncf1['balance']['P']['to_water'][1,:, :]
PL1_meana = np.mean(PL1[:,:]) 
PL2 = ncf2['balance']['P']['to_water'][1,:, :]
PL2_meana = np.mean(PL2[:,:]) 

##Nleach=np.mean([NL0_meana,NL1_meana,NL2_meana]) #vai sittenkin kuin alla?
Nleach=(NL0_meana*lastind+NL1_meana*lastind1+NL2_meana*lastind2)/(lastind+lastind1+lastind2)
Pleach=(PL0_meana*lastind+PL1_meana*lastind1+PL2_meana*lastind2)/(lastind+lastind1+lastind2)

standco2 = ncf['balance']['C']['stand_c_balance_co2eq'][1,:, :]
standco2_meana = np.mean(standco2[:,:]) 
standco21 = ncf1['balance']['C']['stand_c_balance_co2eq'][1,:, :]
standco21_meana = np.mean(standco21[:,:]) 
standco22 = ncf2['balance']['C']['stand_c_balance_co2eq'][1,:, :]
standco22_meana = np.mean(standco22[:,:]) 


cb_total_as_CO2=(standco2_meana*lastind+standco21_meana*lastind1+standco22_meana*lastind2)/(lastind+lastind1+lastind2)

#sitten sama muillekin scenaarioille ja kaikille kuvioille
##ja sitten sama mineraalimaiden kuvioille

soilco2 = ncf['balance']['C']['soil_c_balance_co2eq'][1,:, :]
soilco2_meana = np.mean(soilco2[:,:]) 
soilco21 = ncf1['balance']['C']['soil_c_balance_co2eq'][1,:, :]
soilco21_meana = np.mean(soilco21[:,:]) 
soilco22 = ncf2['balance']['C']['soil_c_balance_co2eq'][1,:, :]
soilco22_meana = np.mean(soilco22[:,:]) 

cb_as_CO2=(soilco2_meana*lastind+soilco21_meana*lastind1+soilco22_meana*lastind2)/(lastind+lastind1+lastind2)

#cb_trees_as_CO2=cb_total_as_CO2-cb_as_CO2

for st in stands:
    for sc in scenarios:
        files=glob(resfol+sc+'/*.nc')
        events=[i[(len(resfol)+len(sc)+1+len(st)+1):-3] for i in files if i.startswith(resfol+sc+'/'+st)]
        events.sort()
        for ev in events:
            ff = os.path.join(resfol,sc,st+'_'+ev+'.nc')



soilco2 = ncf['balance']['C']['soil_c_balance_co2eq'][1,:, :]
standco2 = ncf['balance']['C']['stand_c_balance_co2eq'][1,:, :]

vg = ncf['stand']['volumegrowth'][scen,:, :]

towater = ncf['balance'][substance]['to_water'][1,:, :]

wt = np.mean(ncf['strip']['dwtyr'][1,:, :], axis = 0)
cols = np.shape(wt)[0]
df = pd.DataFrame(data=standco2, columns=list(range(cols)))
 
###
fig, ax = plt.subplots(3, 2, figsize=(10,10),gridspec_kw={'width_ratios': [1, 1]}, dpi=300)
ax[0,0].bar(labels, cb_total_as_CO2[0:6], color=colors)
ax[0,0].set_ylabel('$(tn$ $CO_2$$-ekv.$ \n$ha^{-1}$ $a^{-1})$')
ax[0,0].set_xlabel(stands[0])
ax[0,0].set_title('Hiilitase, maaperä&puusto')

ax[0,1].bar(labels, cb_total_as_CO2[6:12], color=colors)
ax[0,1].set_ylabel('$(tn$ $CO_2$$-ekv.$ \n$ha^{-1}$ $a^{-1})$')
ax[0,1].set_xlabel(stands[1])
ax[0,1].set_title('Hiilitase, maaperä&puusto')

ax[1,0].bar(labels, cb_total_as_CO2[12:18], color=colors)
ax[1,0].set_ylabel('$(tn$ $CO_2$$-ekv.$ \n$ha^{-1}$ $a^{-1})$')
ax[1,0].set_xlabel(stands[2])
ax[1,0].set_title('Hiilitase, maaperä&puusto')

ax[1,1].bar(labels, cb_total_as_CO2[18:24], color=colors)
ax[1,1].set_ylabel('$(tn$ $CO_2$$-ekv.$ \n$ha^{-1}$ $a^{-1})$')
ax[1,1].set_xlabel(stands[3])
ax[1,1].set_title('Hiilitase, maaperä&puusto')

ax[2,0].bar(labels, cb_total_as_CO2[24:30], color=colors)
ax[2,0].set_ylabel('$(tn$ $CO_2$$-ekv.$ \n$ha^{-1}$ $a^{-1})$')
ax[2,0].set_xlabel(stands[4])
ax[2,0].set_title('Hiilitase, maaperä&puusto')

ax[2,1].bar(labels, cb_total_as_CO2[30:36], color=colors)
ax[2,1].set_ylabel('$(tn$ $CO_2$$-ekv.$ \n$ha^{-1}$ $a^{-1})$')
ax[2,1].set_xlabel(stands[5])
ax[2,1].set_title('Hiilitase, maaperä&puusto')
plt.tight_layout()
fig.savefig(r'/scratch/project_2002470/HIILIPOLKU_neuvonta/carbon_balances_peatsites.png', dpi=300)

###
###
fig, ax = plt.subplots(3, 2, figsize=(10,10),gridspec_kw={'width_ratios': [1, 1]}, dpi=300)
ax[0,0].bar(labels, nleachs[0:6], color=colors)
ax[0,0].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[0,0].set_xlabel(stands[0])
ax[0,0].set_title('N-huuhtouma')

ax[0,1].bar(labels, nleachs[6:12], color=colors)
ax[0,1].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[0,1].set_xlabel(stands[1])
ax[0,1].set_title('N-huuhtouma')

ax[1,0].bar(labels, nleachs[12:18], color=colors)
ax[1,0].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[1,0].set_xlabel(stands[2])
ax[1,0].set_title('N-huuhtouma')

ax[1,1].bar(labels, nleachs[18:24], color=colors)
ax[1,1].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[1,1].set_xlabel(stands[3])
ax[1,1].set_title('N-huuhtouma')

ax[2,0].bar(labels, nleachs[24:30], color=colors)
ax[2,0].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[2,0].set_xlabel(stands[4])
ax[2,0].set_title('N-huuhtouma')

ax[2,1].bar(labels, nleachs[30:36], color=colors)
ax[2,1].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[2,1].set_xlabel(stands[5])
ax[2,1].set_title('N-huuhtouma')

plt.tight_layout()
fig.savefig(r'/scratch/project_2002470/HIILIPOLKU_neuvonta/nleach_peatsites.png', dpi=300)

###
fig, ax = plt.subplots(3, 2, figsize=(10,10),gridspec_kw={'width_ratios': [1, 1]}, dpi=300)
ax[0,0].bar(labels, pleachs[0:6], color=colors)
ax[0,0].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[0,0].set_xlabel(stands[0])
ax[0,0].set_title('P-huuhtouma')

ax[0,1].bar(labels, pleachs[6:12], color=colors)
ax[0,1].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[0,1].set_xlabel(stands[1])
ax[0,1].set_title('P-huuhtouma')

ax[1,0].bar(labels, pleachs[12:18], color=colors)
ax[1,0].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[1,0].set_xlabel(stands[2])
ax[1,0].set_title('P-huuhtouma')

ax[1,1].bar(labels, pleachs[18:24], color=colors)
ax[1,1].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[1,1].set_xlabel(stands[3])
ax[1,1].set_title('P-huuhtouma')

ax[2,0].bar(labels, pleachs[24:30], color=colors)
ax[2,0].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[2,0].set_xlabel(stands[4])
ax[2,0].set_title('P-huuhtouma')

ax[2,1].bar(labels, pleachs[30:36], color=colors)
ax[2,1].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[2,1].set_xlabel(stands[5])
ax[2,1].set_title('P-huuhtouma')

plt.tight_layout()
fig.savefig(r'/scratch/project_2002470/HIILIPOLKU_neuvonta/pleach_peatsites.png', dpi=300)



###
###
###
fig, ax = plt.subplots(3, 2, figsize=(10,10),gridspec_kw={'width_ratios': [1, 1]}, dpi=300)
ax[0,0].bar(labels, gr[0:6], color=colors)
ax[0,0].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[0,0].set_xlabel(stands[0])
ax[0,0].set_title('Puuston kasvu')

ax[0,1].bar(labels, gr[6:12], color=colors)
ax[0,1].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[0,1].set_xlabel(stands[1])
ax[0,1].set_title('Puuston kasvu')

ax[1,0].bar(labels, gr[12:18], color=colors)
ax[1,0].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[1,0].set_xlabel(stands[2])
ax[1,0].set_title('Puuston kasvu')

ax[1,1].bar(labels, gr[18:24], color=colors)
ax[1,1].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[1,1].set_xlabel(stands[3])
ax[1,1].set_title('Puuston kasvu')

ax[2,0].bar(labels, gr[24:30], color=colors)
ax[2,0].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[2,0].set_xlabel(stands[4])
ax[2,0].set_title('Puuston kasvu')

ax[2,1].bar(labels, gr[30:36], color=colors)
ax[2,1].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[2,1].set_xlabel(stands[5])
ax[2,1].set_title('Puuston kasvu')

plt.tight_layout()
fig.savefig(r'/scratch/project_2002470/HIILIPOLKU_neuvonta/growth_peatsites.png', dpi=300)

###
fig, ax = plt.subplots(3, 2, figsize=(10,10),gridspec_kw={'width_ratios': [1, 1]}, dpi=300)
ax[0,0].bar(labels, hakkuut[0:6], color=colors)
ax[0,0].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[0,0].set_xlabel(stands[0])
ax[0,0].set_title('Hakkuut')

ax[0,1].bar(labels, hakkuut[6:12], color=colors)
ax[0,1].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[0,1].set_xlabel(stands[1])
ax[0,1].set_title('Hakkuut')

ax[1,0].bar(labels, hakkuut[12:18], color=colors)
ax[1,0].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[1,0].set_xlabel(stands[2])
ax[1,0].set_title('Hakkuut')

ax[1,1].bar(labels, hakkuut[18:24], color=colors)
ax[1,1].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[1,1].set_xlabel(stands[3])
ax[1,1].set_title('Hakkuut')

ax[2,0].bar(labels, hakkuut[24:30], color=colors)
ax[2,0].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[2,0].set_xlabel(stands[4])
ax[2,0].set_title('Hakkuut')

ax[2,1].bar(labels, hakkuut[30:36], color=colors)
ax[2,1].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[2,1].set_xlabel(stands[5])
ax[2,1].set_title('Hakkuut')

plt.tight_layout()
fig.savefig(r'/scratch/project_2002470/HIILIPOLKU_neuvonta/hakkuut_peatsites.png', dpi=300)

###

fig, ax = plt.subplots(3, 2, figsize=(10,10),gridspec_kw={'width_ratios': [1, 1]}, dpi=300)
ax[0,0].bar(labels, vols_inend[0:6], color=colors)
ax[0,0].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[0,0].set_xlabel(stands[0])
ax[0,0].set_title('Puuston tilavuus 50 v kuluttua')

ax[0,1].bar(labels, vols_inend[6:12], color=colors)
ax[0,1].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[0,1].set_xlabel(stands[1])
ax[0,1].set_title('Puuston tilavuus 50 v kuluttua')

ax[1,0].bar(labels, vols_inend[12:18], color=colors)
ax[1,0].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[1,0].set_xlabel(stands[2])
ax[1,0].set_title('Puuston tilavuus 50 v kuluttua')

ax[1,1].bar(labels, vols_inend[18:24], color=colors)
ax[1,1].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[1,1].set_xlabel(stands[3])
ax[1,1].set_title('Puuston tilavuus 50 v kuluttua')

ax[2,0].bar(labels, vols_inend[24:30], color=colors)
ax[2,0].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[2,0].set_xlabel(stands[4])
ax[2,0].set_title('Puuston tilavuus 50 v kuluttua')

ax[2,1].bar(labels, vols_inend[30:36], color=colors)
ax[2,1].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[2,1].set_xlabel(stands[5])
ax[2,1].set_title('Puuston tilavuus 50 v kuluttua')

plt.tight_layout()
fig.savefig(r'/scratch/project_2002470/HIILIPOLKU_neuvonta/puustontilavuus_lopussa_peatsites.png', dpi=300)


###
fig, ax = plt.subplots(4, 2, figsize=(5.5,7),gridspec_kw={'width_ratios': [1, 1]}, dpi=300)

ax[0,0].bar(labels, df['v_end'][0:7], color=colors)
ax[0,0].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[0,0].set_title('Puuston tilavuus 50 v kuluttua')

 

ax[0,1].bar(labels, df['totgr'][0:7], color=colors)
ax[0,1].set_ylabel('$(m^3$ $ha^{-1}$)')
ax[0,1].set_title('Puuston kasvu')

 

ax[2,0].bar(labels, df['Nleach'][0:7], color=colors)
ax[2,0].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[2,0].set_title('N-huuhtouma')

 

ax[1,0].bar(labels, df['Nleach'][0:7], color=colors) #v_cut
ax[1,0].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$') #$(m^3$ $ha^{-1}$)
ax[1,0].set_title('N-huuhtouma') #Hakkuut

 


ax[1,1].bar(labels, df['cbt_as_CO2'][0:7]/-1000., color=colors)
ax[1,1].set_ylabel('$(tn$ $CO_2$$-ekv.$ \n$ha^{-1}$ $a^{-1})$')
ax[1,1].set_title('Puuston hiilitase')

 

ax[3,0].bar(labels, df['Pleach'][0:7], color=colors)
ax[3,0].set_ylabel('$(kg$ $ha^{-1}$ $a^{-1})$')
ax[3,0].set_title('P-huuhtouma')

 

ax[2,1].bar(labels, df['cb_as_CO2'][0:7]/-1000., color=colors)
ax[2,1].set_ylabel('$(tn$ $CO_2$$-ekv.$ \n$ha^{-1}$ $a^{-1})$')
ax[2,1].set_title('Maaperän hiilitase')

 

ax[3,1].bar(labels, df['cb_total_as_CO2'][0:7]/-1000., color=colors)
ax[3,1].set_ylabel('$(tn$ $CO_2$$-ekv.$ \n$ha^{-1}$ $a^{-1})$')
ax[3,1].set_title('Hiilitase, maaperä&puusto')

 


plt.tight_layout()"""

#yksittaisen ncf tarkastelu
#scenarios=['Kuonanjoki_BAU_A', 'Kuonanjoki_BAU_B', 'Kuonanjoki_BIO_A', 'Kuonanjoki_BIO_B', 'Kuonanjoki_HII_A', 'Kuonanjoki_HII_B']
#scenarios=['Sorvasranta_BAU_A', 'Sorvasranta_BAU_B', 'Sorvasranta_BIO_A', 'Sorvasranta_BIO_B', 'Sorvasranta_HII_A', 'Sorvasranta_HII_B']
scenarios=['Halvanjoki_BAU_A', 'Halvanjoki_BAU_B', 'Halvanjoki_BIO_A', 'Halvanjoki_BIO_B', 'Halvanjoki_HII_A', 'Halvanjoki_HII_B']


#scenario=scenarios[0]#
#mottifol=r'/scratch/project_2002470/HIILIPOLKU_data/'
#resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'
#sc=scenarios[0]
#files=glob(resfol+sc+'/*.nc')
#st=stands[0]
#events=[i[(len(resfol)+len(sc)+1+len(st)+1):-3] for i in files if i.startswith(resfol+sc+'/'+st)]
#events.sort()
#ev=events[0]
#for ev in events:
#    sto_values=peat_ini_end(resfol+scenario+"/"+f+"_"+ev+".nc", dscen[ev])
#    storages.append(sto_values)
"""
for st in stands:
    for sc in scenarios:
        files=glob(resfol+sc+'/*.nc')
        events=[i[(len(resfol)+len(sc)+1+len(st)+1):-3] for i in files if i.startswith(resfol+sc+'/'+st)]
        events.sort()
        for ev in events:
            ff = os.path.join(resfol,sc,st+'_'+ev+'.nc')
            ncf=Dataset(ff, mode='r')
            dscen=1 #nää pitäis tsekata tosin et onko kaikissa tuo 1 vaan eli ei kunnostusojitusta
            wt = ncf['strip']['dwtyr'][dscen,:, :]
            esomp1= ncf['esom']['Mass']['P1'][dscen]
            baltowat= ncf['balance']['N']['to_water'][1]
            if esomp1.mask:# ei toimi jos mask löytyy
               print("nomask", esomp1.data)    
            indices=np.where(~ncf['balance']['N']['to_water'][1].mask)
            lastind=max(indices[0])
            if lastind>0:
                print(sc, st, "ok")
            else:
                print(sc, st, ev,  lastind,baltowt)


ffe = os.path.join(r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/outputs_1505/',sc,st+'_'+ev+'.nc')
ncfe=Dataset(ffe, mode='r')
wte = ncfe['strip']['dwtyr'][dscen,:, :]

1773544
1794219
1794384
1794415
1794416
24184122
24184740
24289196
24289706
24290070
24290595
24290600
24290601
24290773
24291015
24291078
24291080
24291407
24291422
24292574
24292576
24293035
24293036
24293309
24293355
24293360
24293389
24293876
24293889
24293910
24294390
24294418
24294477
24294971
25736742
33367261
35394796
35710103
35710136
35715397
35715967
35716692


"""
#car=carbon(ff, 1,r'/scratch/project_2002470/HIILIPOLKU_neuvonta/',st,sc)
#compare_scens(ff)
# hydrology(ff, 0)
# stand(ff, 0)
# mass(ff, 0)
#carbon(ff,1)
# nutrient_balance(ff, 'N', 0)


#scenario=scenarios[0]
missing_stands_list=[]
newrunslist=[]
missing_stands_list_motti=[]

def check_mottifiles(scenarios):
    for scenario in scenarios:
        mottistands=pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/'+scenario[:-6]+'/'+scenario+'/motti/'+scenario+'_kuviot.csv', sep=';', encoding='latin1')
        ls_stands=list(set(mottistands['kuvio'][mottistands['kuivatustilanne']>=3]))
        ls_stands2=[str(i) for i in ls_stands]
        ls_stands2= [*set(ls_stands2)]

        mottifiles=glob(mottifol+scenario[:-6]+'/'+scenario+'/*.xls')
        stands_mf=[i[(len(mottifol)+len(scenario[:-6])+1+len(scenario)+1):-7] for i in mottifiles]
        stands_mf = [*set(stands_mf)]
        stands_mf.sort()
        if len(ls_stands2)==len(stands_mf):
            print(scenario, " ok")
        else:
            missing_stands=list(set(ls_stands) - set(stands_mf))   #list(set(common_ms).intersection(dtw_list))
            missing_stands_list_motti.append(missing_stands)
            missing_stands_list_motti.sort()
            print(scenario, len(ls_stands)-len(stands_mf),"missing stands: ", missing_stands) 
    return missing_stands_list_motti


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


"""
#scenario=['Sorvasranta_BAU_A', 'Sorvasranta_BAU_B']
#scenarios=['Sorvasranta_BAU_A', 'Sorvasranta_BAU_B', 'Sorvasranta_BIO_A', 'Sorvasranta_BIO_B','Sorvasranta_HII_A', 'Sorvasranta_HII_B']

resfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'

"""
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

def find_i(scenario, problem_stands):
    #fprlines=[]
    flines=[]
    #with open('/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'+scenario+'_res_prob.txt', "a") as myprfile:
    #    for myline in myprfile:
    #        fprlines.append(myline.rstrip('\n'))  
    with open(r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/'+scenario+'/'+scenario+'.txt','rt') as myfile:
                for myline in myfile:
                    flines.append(myline.rstrip('\n'))
    #nmlist=["i="+str(i) for i in range(len(stands))]
    standis=[]
    for line in flines:
        for i in problem_stands:
            if line.find(i)!=-1:
                standis.append(line[6+len(i)+4:])
    print(standis)
    standis = [*set(standis)]
    standis = [int(i) for i in standis] 
    standis.sort()
    return standis


#scenarios=['Kuonanjoki_BAU_A', 'Kuonanjoki_BAU_B', 'Kuonanjoki_BIO_A', 'Kuonanjoki_BIO_B', 'Kuonanjoki_HII_A', 'Kuonanjoki_HII_B']

#k_baua_prs=checkres(scenarios[0])
#k_baub_prs=checkres(scenarios[1])
#k_bioa_prs=checkres(scenarios[2])
#k_biob_prs=checkres(scenarios[3])
#k_hiia_prs=checkres(scenarios[4])
#k_hiib_prs=checkres(scenarios[5])
"""
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
"""
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
"""
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

