# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 11:53:12 2021

@author: Leena Stenberg / LUKE
"""


import numpy as np
import pandas as pd
import os
import shutil
from glob import glob


"""
#cd /scratch/project_2002470/HIILIPOLKU_data/
#unzip -d Sorvasranta Sorvasranta_220323.zip
#unzip -d Kuonanjoki Kuonanjoki_220323.zip
#unzip -d Halvanjoki Halvanjoki_220323.zip

mkdir Sorvasranta/Sorvasranta_BAU_A;mkdir Sorvasranta/Sorvasranta_BAU_B; mkdir Sorvasranta/Sorvasranta_BIO_A;mkdir Sorvasranta/Sorvasranta_BIO_B;
mkdir Sorvasranta/Sorvasranta_HII_A;mkdir Sorvasranta/Sorvasranta_HII_B;
mkdir Sorvasranta/Sorvasranta_BAU_A/motti;mkdir Sorvasranta/Sorvasranta_BAU_B/motti; mkdir Sorvasranta/Sorvasranta_BIO_A/motti; mkdir Sorvasranta/Sorvasranta_BIO_B/motti;
mkdir Sorvasranta/Sorvasranta_HII_A/motti;mkdir Sorvasranta/Sorvasranta_HII_B/motti;
cp Sorvasranta/Sorvasranta_BAU_A.xlsx Sorvasranta/Sorvasranta_BAU_A/motti/;
cp Sorvasranta/Sorvasranta_BAU_B.xlsx Sorvasranta/Sorvasranta_BAU_B/motti/;
cp Sorvasranta/Sorvasranta_BIO_A.xlsx Sorvasranta/Sorvasranta_BIO_A/motti/;
cp Sorvasranta/Sorvasranta_BIO_B.xlsx Sorvasranta/Sorvasranta_BIO_B/motti/;
cp Sorvasranta/Sorvasranta_HII_A.xlsx Sorvasranta/Sorvasranta_HII_A/motti/;
cp Sorvasranta/Sorvasranta_HII_B.xlsx Sorvasranta/Sorvasranta_HII_B/motti/

mkdir Kuonanjoki/Kuonanjoki_BAU_A;mkdir Kuonanjoki/Kuonanjoki_BAU_B; mkdir Kuonanjoki/Kuonanjoki_BIO_A;mkdir Kuonanjoki/Kuonanjoki_BIO_B;
mkdir Kuonanjoki/Kuonanjoki_HII_A;mkdir Kuonanjoki/Kuonanjoki_HII_B;
mkdir Kuonanjoki/Kuonanjoki_BAU_A/motti;mkdir Kuonanjoki/Kuonanjoki_BAU_B/motti; mkdir Kuonanjoki/Kuonanjoki_BIO_A/motti;mkdir Kuonanjoki/Kuonanjoki_BIO_B/motti;
mkdir Kuonanjoki/Kuonanjoki_HII_A/motti;mkdir Kuonanjoki/Kuonanjoki_HII_B/motti;
cp Kuonanjoki/Kuonanjoki_BAU_A.xlsx Kuonanjoki/Kuonanjoki_BAU_A/motti/;
cp Kuonanjoki/Kuonanjoki_BAU_B.xlsx Kuonanjoki/Kuonanjoki_BAU_B/motti/;
cp Kuonanjoki/Kuonanjoki_BIO_A.xlsx Kuonanjoki/Kuonanjoki_BIO_A/motti/;
cp Kuonanjoki/Kuonanjoki_BIO_B.xlsx Kuonanjoki/Kuonanjoki_BIO_B/motti/;
cp Kuonanjoki/Kuonanjoki_HII_A.xlsx Kuonanjoki/Kuonanjoki_HII_A/motti/;
cp Kuonanjoki/Kuonanjoki_HII_B.xlsx Kuonanjoki/Kuonanjoki_HII_B/motti/

mkdir Halvanjoki/Halvanjoki_BAU_A;mkdir Halvanjoki/Halvanjoki_BAU_B; mkdir Halvanjoki/Halvanjoki_BIO_A;mkdir Halvanjoki/Halvanjoki_BIO_B;
mkdir Halvanjoki/Halvanjoki_HII_A;mkdir Halvanjoki/Halvanjoki_HII_B;
mkdir Halvanjoki/Halvanjoki_BAU_A/motti;mkdir Halvanjoki/Halvanjoki_BAU_B/motti; mkdir Halvanjoki/Halvanjoki_BIO_A/motti;mkdir Halvanjoki/Halvanjoki_BIO_B/motti;
mkdir Halvanjoki/Halvanjoki_HII_A/motti;mkdir Halvanjoki/Halvanjoki_HII_B/motti;

cp Halvanjoki/Halvanjoki_BAU_A.xlsx Halvanjoki/Halvanjoki_BAU_A/motti/;
cp Halvanjoki/Halvanjoki_BAU_B.xlsx Halvanjoki/Halvanjoki_BAU_B/motti/;
cp Halvanjoki/Halvanjoki_BIO_A.xlsx Halvanjoki/Halvanjoki_BIO_A/motti/;
cp Halvanjoki/Halvanjoki_BIO_B.xlsx Halvanjoki/Halvanjoki_BIO_B/motti/;
cp Halvanjoki/Halvanjoki_HII_A.xlsx Halvanjoki/Halvanjoki_HII_A/motti/;
cp Halvanjoki/Halvanjoki_HII_B.xlsx Halvanjoki/Halvanjoki_HII_B/motti/

#in python 
scenarios=['Sorvasranta_BAU_A', 'Sorvasranta_BAU_B', 'Sorvasranta_BIO_A', 'Sorvasranta_BIO_B', 'Sorvasranta_HII_A', 'Sorvasranta_HII_B']
sheets=['kuviot', 'puustot', 'poistumat', 'tapahtumat', 'kasvut']

for scen in scenarios:
    for sh in sheets:
        ifile = r'/scratch/project_2002470/HIILIPOLKU_data/Sorvasranta/'+scen+'/motti/'+scen+'.xlsx'
        df = pd.read_excel(ifile, sheet_name=sh)
        df.to_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Sorvasranta/'+scen+'/motti/'+scen+'_'+sh+'.csv', sep=';', encoding='latin1')

create_motti_files_silvi('hiilipolku', scenarios[0])
create_motti_files_silvi('hiilipolku', scenarios[1]);create_motti_files_silvi('hiilipolku', scenarios[2]);create_motti_files_silvi('hiilipolku', scenarios[3]);create_motti_files_silvi('hiilipolku', scenarios[4]);create_motti_files_silvi('hiilipolku', scenarios[5])


scenarios=['Kuonanjoki_BAU_A', 'Kuonanjoki_BAU_B', 'Kuonanjoki_BIO_A', 'Kuonanjoki_BIO_B', 'Kuonanjoki_HII_A', 'Kuonanjoki_HII_B']
for scen in scenarios:
    for sh in sheets:
        ifile = r'/scratch/project_2002470/HIILIPOLKU_data/Kuonanjoki/'+scen+'/motti/'+scen+'.xlsx'
        df = pd.read_excel(ifile, sheet_name=sh)
        df.to_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Kuonanjoki/'+scen+'/motti/'+scen+'_'+sh+'.csv', sep=';', encoding='latin1')

create_motti_files_silvi('hiilipolku', scenarios[0]); create_motti_files_silvi('hiilipolku', scenarios[1]);create_motti_files_silvi('hiilipolku', scenarios[2]);create_motti_files_silvi('hiilipolku', scenarios[3]);create_motti_files_silvi('hiilipolku', scenarios[4]);create_motti_files_silvi('hiilipolku', scenarios[5])
######################
###############
#############
############
#######



scenarios=['Halvanjoki_BAU_A', 'Halvanjoki_BAU_B', 'Halvanjoki_BIO_A', 'Halvanjoki_BIO_B', 'Halvanjoki_HII_A', 'Halvanjoki_HII_B']
for scen in scenarios:
    for sh in sheets:
        ifile = r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/'+scen+'/motti/'+scen+'.xlsx'
        df = pd.read_excel(ifile, sheet_name=sh)
        df.to_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/'+scen+'/motti/'+scen+'_'+sh+'.csv', sep=';', encoding='latin1')

create_motti_files_silvi('hiilipolku', scenarios[0]); create_motti_files_silvi('hiilipolku', scenarios[1]);create_motti_files_silvi('hiilipolku', scenarios[2]);create_motti_files_silvi('hiilipolku', scenarios[3]);create_motti_files_silvi('hiilipolku', scenarios[4]);create_motti_files_silvi('hiilipolku', scenarios[5])

#####################
for scen in scenarios:
    files=glob(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/*.xls')
    for file in files:
        os.remove(file)

for scen in scenarios:
    files=glob(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/*.txt')
    for file in files:
        os.remove(file)


for scen in scenarios:
    for sh in sheets:
        ifile = r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/Halvanjoki/'+scen+'/motti/'+'tila6tarkistetuinkasittelyohjein.xlsx'
        df = pd.read_excel(ifile, sheet_name=sh)
        if sh != 'kuviot':
            dfsc=df[df['sken']==scen[-5:].lower()]
            dfsc.to_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/Halvanjoki/'+scen+'/motti/'+scen+'_'+sh+'.csv', sep=';', encoding='latin1')
        else:
            df.to_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/Halvanjoki/'+scen+'/motti/'+scen+'_'+sh+'.csv', sep=';', encoding='latin1')


############
create_motti_files_silvi('hiilipolku', scenarios[0])
create_motti_files_silvi('hiilipolku', scenarios[1]);create_motti_files_silvi('hiilipolku', scenarios[2]);create_motti_files_silvi('hiilipolku', scenarios[3]);create_motti_files_silvi('hiilipolku', scenarios[4]);create_motti_files_silvi('hiilipolku', scenarios[5])



#ifile = r'/scratch/project_2002470/HIILIPOLKU_data/Sorvasranta/Sorvasranta_BAU_A/Sorvasranta_BAU_A.xlsx'
#df = pd.read_excel(ifile, sheet_name='kuviot')
#df.to_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Sorvasranta/Sorvasranta_BAU_A/Sorvasranta_BAU_A_kuviot.csv', sep=';')

##sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_baua.sh; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_baub.sh;sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_bioa.sh;sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_biob.sh; sleep 45m; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_hiia.sh; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_hiib.sh; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_baua.sh; sleep 1h; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_baub.sh; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_bioa.sh; sleep 1h; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_biob.sh; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_hiia.sh; sleep 1h; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_hiib.sh; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_baua.sh; sleep 1h; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_baub.sh; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_bioa.sh; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_biob.sh; sleep 1h; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_hiia.sh; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_hiib.sh
sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_baua.sh; sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_baub.sh;sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_bioa.sh;sbatch SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_biob.sh;


import os
import shutil
from glob import glob
import numpy as np
from numpy.ma import masked_array

scenarios=['Sorvasranta_BAU_A', 'Sorvasranta_BAU_B', 'Sorvasranta_BIO_A', 'Sorvasranta_BIO_B', 'Sorvasranta_HII_A', 'Sorvasranta_HII_B']
scen=scenarios[0]
files=glob(r'/scratch/project_2002470/HIILIPOLKU_data/Sorvasranta/'+scen+'/*.xls')

scen=scenarios[5]
files=glob(r'/scratch/project_2002470/HIILIPOLKU_data/'+scen[:-6]+'/'+scen+'/*.xls')
files[0:10]    
problem_files=[]
for fil in files:
    df = pd.read_excel(fil)
    if len(df)<=1:
        problem_files.append(fil)

problem_files

prob_stands=[i[71:-4] for i in problem_files]

"""











# -------------------------------------------------------------------------
def create_motti_files_silvi(project, scen, maku=None):    # "ubm" poistaa saman sisältöiset bm rivit poikkeuksena aiempaan      # Huom! ohje>1
    '''
    Modified: January 2023
        
    Create Motti input files for Susi simulations from "Skene-Motti" results.
    
    Works technically, but contains many small issues that should be fixed.
    
    project = 'life' (Hydrology Life -project) or 'suo' (SUO-project)
    project-parameter affects only main tree species (and column names)!
    
    
    Modified 8.3.2022:
        - korjattu pÃ¤Ã¤tehakkuuseen liittyviÃ¤ ongelmia; hdom pitÃ¤isi nyt saada
        jÃ¤rkeviÃ¤ arvoja, biomassojen ei pitÃ¤isi jÃ¤Ã¤dÃ¤ pÃ¤Ã¤tehakkuun jÃ¤lkeen nollaksi
        (ainakaan toisen pÃ¤Ã¤tehakkuun jÃ¤lkeen)
        
    Issues to fix:
        - Jos pÃ¤Ã¤tehakkuu tulee vuonna 2, pitÃ¤Ã¤ vuosille 0 ja 1 hakea muut komponentit
        totvol perusteella: hdom saa tÃ¤ssÃ¤ liian suuria arvoja. Tuotos jÃ¤Ã¤ nollaksi?!
        - See comments below...
        
    
    '''
    
    
    mottifolder = r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/'
        
    outfile = mottifolder + '/out.txt'
    
    project = 'hiilipolku'
    
    kuviot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', encoding='latin1', sep=';')
    puustot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_puustot.csv', encoding='latin1', sep=';')
    poistumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_poistumat.csv', encoding='latin1', sep=';')
    tapahtumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_tapahtumat.csv', encoding='latin1', sep=';', index_col=False)
    tulot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kasvut.csv', encoding='latin1', sep=';')
    
    
    if (project == 'life') | (project =='vmi') | (project =='hiilipolku'):
        kuviot = kuviot.rename(columns={'kuvio':'KUVIO'})
        puustot = puustot.rename(columns={'kuvio':'KUVIO',
                                          '_2mm_juuretkg':'>2mm_juuretkg',})
        poistumat = poistumat.rename(columns={'kuvio':'KUVIO'})
        tapahtumat = tapahtumat.rename(columns={'kuvio':'KUVIO',
                                      '_2003': '2003',
                                      '_2004': '2004',
                                      '_2005': '2005'})
        tulot = tulot.rename(columns={'kuvio':'KUVIO'})
    
    
    kuviolista = list(set(kuviot['KUVIO'][kuviot['kuivatustilanne']>=3])) #tassa valitaan 3 Ojitettu kangas ja muut?DRAINAGESTATE	, ojitettu kangas jos pÃ¤Ã¤ryhmÃ¤1 ja alaryhmÃ¤ 3 ja kuivatustilanne 5?? 1 Ojittamaton kangas, 2 Soistunut kangas,3 Ojitettu kangas,	6 luonnontilainen suo,	7 Ojikko, 8	Muuttuma,9	Turvekangas
    
    
    for i in range(0,len(kuviolista)):
    
        print(i, 'stand: ', kuviolista[i])
        
        # alta poistettu ohjeeseen ja skenaarioon viittaavat rajaukset,
        # oletetaan, ettÃ¤ hiilipolussa vain yksi ohje per kuvio
        
        sdata0 = puustot[(puustot['KUVIO']==kuviolista[i])]
        pdata = poistumat[(poistumat['KUVIO']==kuviolista[i]) & (poistumat['HARV']>1)] #mita harv koodit meinaa?
        pdata2 = poistumat[(poistumat['KUVIO']==kuviolista[i])]
        tdata = tapahtumat[(tapahtumat['KUVIO']==kuviolista[i])]
        tulotdata = tulot[(tulot['KUVIO']==kuviolista[i])]
        kdata = kuviot[kuviot['KUVIO']==kuviolista[i]]
        
        sdata0 = sdata0.reset_index(drop=True)
        pdata = pdata.reset_index(drop=True)
        pdata2 = pdata2.reset_index(drop=True)
        tdata = tdata.reset_index(drop=True)
        tulotdata = tulotdata.reset_index(drop=True)
        kdata = kdata.reset_index(drop=True)
        
        tdata = tdata[(tdata['2003']>-1) | (tdata['2004']>-1) | (tdata['2005']>-1)]
        
        ensih = tdata[tdata['2003']>-1]['2003'].values
        harv = tdata[tdata['2004']>-1]['2004'].values
        paateh = tdata[tdata['2005']>-1]['2005'].values
        #voisi kayda tsekkaamassa onko taimikonhoito tai varhaisperkaus niin sitten ei muuttaisi puulajia tai ottaisi sen muualta et sailyy kuusena tai mantyna? 
        #tdic={'_2001':'lannoitus','_2002':'kunnostusojitus','_2003':'ensiharvennus','_2004':'harvennus','_2005':'paatehakkuu','_2006':'varhaiperkaus','_2007':'taimikonhoito','_2008':'viljely','_2009':'luontainen uudistuminen'}
        # taimh = tdata[tdata['_2007']>-1]['_2007'].values
        # Pudotetaan viimeisen vuoden tapahtumat pois, koska sekoittavat
        # myÃ¶hemmin tulevan koodin
        ensih = ensih[ensih<100]
        harv = harv[harv<100]
        paateh = paateh[paateh<100]
        
        vuodet = np.unique(np.sort(np.concatenate((ensih, harv, paateh))))
        
        ap = 0
        
        pera = -1
        
        # Tsekataan, onko tapahtumissa perÃ¤kkÃ¤isiÃ¤ vuosia ja poistetaan niistÃ¤ jÃ¤lkimmÃ¤inen
        for v in range(0,len(vuodet)-1):
            if vuodet[v+1]-vuodet[v]==1:
                
                ensih = np.delete(ensih, np.where(ensih == vuodet[v+1]))
                harv = np.delete(harv, np.where(harv == vuodet[v+1]))
                paateh = np.delete(paateh, np.where(paateh == vuodet[v+1]))
                
                vuodet = np.delete(vuodet, [v+1], None)
                
                pera = v
                
                ap = ap + 1
                
            if v+1>=len(vuodet)-1:
                break
        
        if ap>1:
            with open(outfile, "a") as myfile:
                myfile.write('\n\nKaksi kertaa perakkaiset vuodet! stand=' + str(kuviolista[i]) + ', i=' + str(i))
    
        n_mottifiles = len(vuodet) + 1
    
    
        if len(sdata0)==0: # skipataan kuvio, jos dataa ei ole
            continue
    
        if len(sdata0)>21:
            print('sdata0: liikaa rivejÃ¤?!')
        
    
                                    
        # Luodaan pÃ¤Ã¤tehakkuun jÃ¤lkeiset rivit ja korjataan rivin kopioimisesta aiheutuneet
        # luvut
        
        # muokattu 6.8.2021: poistettu tuotos, tuotosKuolleet, tuotosTilavuus comps:sta
        comps = ['DomH','MeanDbh','medianH',\
                          'BA','Kannot_juuretkg',\
                              'tuotosTilavuus',\
                              'TukPL1', 'TukPL2','TukPL3',\
                                  'KuiPL1', 'KuiPL2','KuiPL3',\
                                      'HukPL1', 'HukPL2','HukPL3',\
                                          'VolPL1', 'VolPL2', 'VolPL3']  

        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg'] #tama lisatty tahan elokuu 23 jotta namakin tsekataan

        for h in range(0,len(paateh)):
            
            sdata_paateh= sdata0[sdata0['AIKA']==paateh[h]].reset_index(drop=True)
            
            go = False
            
            if len(sdata_paateh)>0:
                if sdata_paateh.loc[0, 'Age']>5:
                    go = True
        
            if (len(sdata_paateh)==0) | (go==True): # Luodaan rivit vain, jos ensimmÃ¤istÃ¤ riviÃ¤ ei ole mukana
                
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>paateh[h]), None)
                
                #sdata0 = sdata0.append([sdata0.iloc[paaterow]], ignore_index=True) # Rivi menee taulukon loppuun
                paaterow_df=pd.DataFrame([sdata0.iloc[paaterow]])
                sdata0 = pd.concat((sdata0,paaterow_df),axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = paateh[h]+1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==paateh[h]+1), None)
                
                # sdata0.loc[paaterow,'Age'] = 1 # TÃ¤mÃ¤ korvattu alla olevalla
                
                if paaterow>0:
                    if sdata0.loc[paaterow-1, 'Age']<5:
                        sdata0.loc[paaterow,'Age'] = sdata0.loc[paaterow-1, 'Age']+1
                    else:
                        sdata0.loc[paaterow,'Age'] = sdata0.loc[paaterow+1, 'Age']-1
                
                sdata0.loc[paaterow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][paaterow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[paaterow,comps[c]]>0:
                        sdata0.loc[paaterow,comps[c]] = sdata0['TotVol'][paaterow]/sdata0['TotVol'][paaterow+1]*sdata0.loc[paaterow+1, comps[c]]
                        print(kuviolista[i],comps[c], "rivi 303 ongelma")

            
                if sdata0['tuotosKuolleet'][paaterow] < sdata0['TotVol'][paaterow]:
                    sdata0.loc[paaterow,'N'] = (1 + sdata0['tuotosKuolleet'][paaterow]/sdata0['TotVol'][paaterow])*sdata0['N'][paaterow]        
                    
                    if paaterow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][paaterow] > sdata0['N'][paaterow-1]:
                            if paaterow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[paaterow,'N'] = (sdata0['N'][paaterow-1] + sdata0['N'][paaterow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[paaterow,'N'] = 0.99*sdata0['N'][paaterow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.
                            
                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[paaterow,biom_comp[c]]>0: 
                        sdata0.loc[paaterow,biom_comp[c]] = sdata0['TotVol'][paaterow]/sdata0['TotVol'][paaterow+1]*sdata0.loc[paaterow+1, biom_comp[c]] ####
                
            
        # Luodaan ensiharvennuksen jÃ¤lkeiset rivit (jos vuosi ei osu kohdallaan muuten)
        # EsimerkissÃ¤ nÃ¤ytti siltÃ¤, ettÃ¤ ensiharvennusvuoden riville tulostuu harventamaton tilavuus! Siksi tÃ¤ssÃ¤ kÃ¤sitelty riviÃ¤+1
        # Toisessa esimerkissÃ¤ tulostui harvennettu tilavuus! Joten vaihtelee.
        
        for h in range(0,len(ensih)):
    
    
            ensihrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==ensih[h]+1), None)
            ensihrow2 = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>ensih[h]), None)
            
            if ensihrow is None:
                
                #sdata0 = sdata0.append([sdata0.iloc[ensihrow2]], ignore_index=True) # Rivi menee taulukon loppuun
                ensihrow2_df=pd.DataFrame([sdata0.iloc[ensihrow2]])
                sdata0 = pd.concat((sdata0,ensihrow2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = ensih[h]+1
                
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                ensihrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==ensih[h]+1), None)
                
                # Jos ikÃ¤ lÃ¶ytyy riviÃ¤ alaempaa:
                if sdata0.loc[ensihrow+1,'Age']>0:
                    sdata0.loc[ensihrow,'Age'] = sdata0.loc[ensihrow+1,'Age']-(sdata0['AIKA'][ensihrow+1]-sdata0['AIKA'][ensihrow])
                # Muussa tapauksessa ikÃ¤ lÃ¶ytyy riviÃ¤ ylempÃ¤Ã¤:
                else:
                    if ensihrow>0:
                        sdata0.loc[ensihrow,'Age'] = sdata0.loc[ensihrow-1,'Age']+(sdata0['AIKA'][ensihrow]-sdata0['AIKA'][ensihrow-1])
                    else:
                        # Jos ei lÃ¶ydy kummaltakaan puolelta:
                        with open(outfile, "a") as myfile:
                            myfile.write('\nEnsiharvennuksen jÃ¤lkeistÃ¤ ikÃ¤Ã¤ ei lÃ¶ydy! kuvio=' + str(kuviolista[i]) + ' i:'+ str(i))
                        
                    
                sdata0.loc[ensihrow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][ensihrow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[ensihrow,comps[c]]>0:
                        sdata0.loc[ensihrow,comps[c]] = sdata0['TotVol'][ensihrow]/sdata0['TotVol'][ensihrow+1]*sdata0.loc[ensihrow+1, comps[c]]
                
                # PÃ¤ivitetÃ¤Ã¤n N
                if sdata0['tuotosKuolleet'][ensihrow] < sdata0['TotVol'][ensihrow]:
                    sdata0.loc[ensihrow,'N'] = (1 + sdata0['tuotosKuolleet'][ensihrow]/sdata0['TotVol'][ensihrow])*sdata0['N'][ensihrow]        
    
                    if ensihrow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][ensihrow] > sdata0['N'][ensihrow-1]:
                            if ensihrow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[ensihrow,'N'] = (sdata0['N'][ensihrow-1] + sdata0['N'][ensihrow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[ensihrow,'N'] = 0.99*sdata0['N'][ensihrow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.

                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[ensihrow,biom_comp[c]]>0: 
                        sdata0.loc[ensihrow,biom_comp[c]] = sdata0['TotVol'][ensihrow]/sdata0['TotVol'][ensihrow+1]*sdata0.loc[ensihrow+1, biom_comp[c]] ####
                            
        
        # Luodaan harvennuksen jÃ¤lkeiset rivit (jos vuosi ei osu kohdalleen muuten)
        # Harvennusvuoden rivillÃ¤ on harventamaton tilavuus. Siksi kÃ¤sitelty riviÃ¤+1.
        
        for h in range(0,len(harv)):
    
    
            hrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==harv[h]+1), None)
            hrow2 = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>harv[h]), None)
            
            if hrow is None:
                
                #sdata0 = sdata0.append([sdata0.iloc[hrow2]], ignore_index=True) # Rivi menee taulukon loppuun
                hrow2_df=pd.DataFrame([sdata0.iloc[hrow2]])
                sdata0 = pd.concat((sdata0,hrow2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = harv[h]+1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                hrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==harv[h]+1), None)
                
                # Jos ikÃ¤ lÃ¶ytyy riviÃ¤ alaempaa:
                if sdata0.loc[hrow+1,'Age']>0:
                    sdata0.loc[hrow,'Age'] = sdata0.loc[hrow+1,'Age']-(sdata0['AIKA'][hrow+1]-sdata0['AIKA'][hrow])
                # Muussa tapauksessa ikÃ¤ lÃ¶ytyy riviÃ¤ ylempÃ¤Ã¤:
                else:
                    if hrow>0:
                        sdata0.loc[hrow,'Age'] = sdata0.loc[hrow-1,'Age']+(sdata0['AIKA'][hrow]-sdata0['AIKA'][hrow-1])
                    else:
                    # Jos ei lÃ¶ydy kummaltakaan puolelta:
                        with open(outfile, "a") as myfile:
                            myfile.write('\nHarvennuksen jÃ¤lkeistÃ¤ ikÃ¤Ã¤ ei lÃ¶ydy! i=' + str(i) + ' kuvio:'+ str(kuviolista[i]))
                
                    
                sdata0.loc[hrow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][hrow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[hrow,comps[c]]>0:
                        sdata0.loc[hrow,comps[c]] = sdata0['TotVol'][hrow]/sdata0['TotVol'][hrow+1]*sdata0.loc[hrow+1, comps[c]]

                        
                # Korjataan N:
                if sdata0['tuotosKuolleet'][hrow] < sdata0['TotVol'][hrow]:
                    sdata0.loc[hrow,'N'] = (1 + sdata0['tuotosKuolleet'][hrow]/sdata0['TotVol'][hrow])*sdata0['N'][hrow]
                    
                    if hrow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][hrow] > sdata0['N'][hrow-1]:
                            if hrow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[hrow,'N'] = (sdata0['N'][hrow-1] + sdata0['N'][hrow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[hrow,'N'] = 0.99*sdata0['N'][hrow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.
                            
                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[hrow,biom_comp[c]]>0: 
                        sdata0.loc[hrow,biom_comp[c]] = sdata0['TotVol'][hrow]/sdata0['TotVol'][hrow+1]*sdata0.loc[hrow+1, biom_comp[c]] ####
    
        # Tsekataan, ettei ole ekan rivin TotVol pyÃ¶ristyksen kanssa ongelmaa:
        
        if (sdata0['TotVol'][0]>10.) & (-0.1 < sdata0['TotVol'][1]-sdata0['TotVol'][0]<0):
            sdata0.loc[0,'TotVol'] = sdata0['TotVol'][1]
                    
          
        # Jos yllÃ¤ olevien korjausten jÃ¤lkeekin aika on alussa > 1,
        # luodaan aika=0 -rivi.
    
        
        repair = 0
        repairB = 0
        
        if (sdata0['AIKA'][0]>1) & (sdata0['Age'][0]>=sdata0['AIKA'][0]): # Jos aika=0 tai aika=1, lisÃ¤tÃ¤Ã¤n alkuun rivi, joka kopioitu ensimmÃ¤iseltÃ¤ olemassa olevalta riviltÃ¤, kÃ¤ytÃ¤nnÃ¶ssÃ¤ usein varhaisperkauksen jÃ¤lkeinen tilanne
            # print('Aika alussa > 1, i=' + str(i))
            with open(outfile, "a") as myfile:
                myfile.write('\nAika alussa > 1, i=' + str(i) + ' kuvio: '+ str(kuviolista[i]))
             
            #sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
            nolla_df=pd.DataFrame([sdata0.iloc[0]])
            sdata0 = pd.concat((sdata0,nolla_df),axis=0, ignore_index=True) # Rivi menee taulukon loppuun
            sdata0['AIKA'][len(sdata0)-1] = 0 # tehdaan alkuun 0 vuosi
            
            sdata0 = sdata0.sort_values(by=['AIKA']) #sorattaan jotta 0-vuosi alkuun
            sdata0 = sdata0.reset_index(drop=True) #resetoidaan indeksi
            sdata0['Age'][0] = sdata0['Age'][0]-(sdata0['AIKA'][1] - sdata0['AIKA'][0]) #tehdaan n-vuodelle oikea ika
    
            
            # if (sdata0['Age'][0]<0) | (sdata0['Age'][1]<0):
            #     sdata0['Age'][0]=80 # Jos ikÃ¤Ã¤ ei saada puustot-tiedostosta, arvotaan 80 v (tÃ¤llÃ¤ ei pitÃ¤isi olla juurikaan merkitystÃ¤ mihinkÃ¤Ã¤n??)
    
                # # print('Second row age = 0 !!! i = ' +str(i))
                # with open(outfile, "a") as myfile:
                #     myfile.write('\nSecond row age = 0 !!! i = ' +str(i) + 'kuvio: '+str(kuviolista[i]))
            
            
            sdata0['TotVol'][0] = tulotdata['volume0'][0] #luodulle aika=0 riville tilavuus tulotdatan vuosittaisita ilavuusista
            
            # PÃ¤ivitetÃ¤Ã¤n muut komponentit
            for c in range(0,len(comps)):
                if sdata0.loc[0,comps[c]]>0: #jos arvo on yli 0
                    sdata0.loc[0,comps[c]] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, comps[c]] ####

            #/users/asalmiva/hiilipolku_leena_silvi_mottitiedostojenteko.py            
            # PÃ¤ivitetÃ¤Ã¤n N
            if sdata0['tuotosKuolleet'][0] < sdata0['TotVol'][0]:
                sdata0.loc[0,'N'] = (1 + sdata0['tuotosKuolleet'][0]/sdata0['TotVol'][0])*sdata0['N'][0]       ##miten tama on laskettu 
                     
            for c in range(0,len(biom_comp)):
                if sdata0.loc[0,biom_comp[c]]>0: #jos arvo on yli 0
                    sdata0.loc[0,biom_comp[c]] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, biom_comp[c]] ####runkohukka pitäisi olla ianakin ennen päätehhakkuta eripäin tuo suhde kuin muilla jten pitäisi skaalata eri tavoin... kasvaa 35 vuotaan saakka ehkä ja sitten laksee ainkain lehtipuilla? tai kuuso+lehti
    
            repair = 1
    
        # # Jos harvennus tulee jo ennen Motti-tulosteen toista riviÃ¤, luodaan AIKA=1 -rivi,
        # # jotta pystytÃ¤Ã¤n tehdÃ¤ Motti-tiedosto, josta saadaan interpolointifunktiot
        # # parille ensimmÃ¤iselle vuodelle
    
        # if len(vuodet)>0:
        #     if vuodet[0]>1: # LisÃ¤ys tehdÃ¤Ã¤n vain, jos hakkuu tulee vuonna 2 tai myÃ¶hemmin
        #         if sdata0['AIKA'][1]>=vuodet[0]: # Jos aika=0 tai aika=1, lisÃ¤tÃ¤Ã¤n alkuun rivi, joka kopioitu ensimmÃ¤iseltÃ¤ olemassa olevalta riviltÃ¤, kÃ¤ytÃ¤nnÃ¶ssÃ¤ varhaisperkauksen jÃ¤lkeinen tilanne
        
        #             with open(outfile, "a") as myfile:
        #                 myfile.write('\nAika[1]>vuodet[0]!!! ---- HUOM!!!, i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
                    
        #             sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
                    
        #             sdata0['AIKA'][len(sdata0)-1] = 1
                    
        #             sdata0 = sdata0.sort_values(by=['AIKA'])
        #             sdata0 = sdata0.reset_index(drop=True)
        #             sdata0['Age'][1] = sdata0['Age'][0]+1
        
                    
        #             # if (sdata0['Age'][1]<0) | (sdata0['Age'][1]<0):
        #             #     sdata0['Age'][1]=81 # Jos ikÃ¤Ã¤ ei saada puustot-tiedostosta, arvotaan 80 v (tÃ¤llÃ¤ ei pitÃ¤isi olla juurikaan merkitystÃ¤ mihinkÃ¤Ã¤n??)
        
        #                 # # print('Second row age = 0 !!! i = ' +str(i))
        #                 # with open(outfile, "a") as myfile:
        #                 #     myfile.write('\nSecond row age = 0 !!! i = ' +str(i))
                    
                    
        #             sdata0['TotVol'][1] = tulotdata['volume1'][0] 
            
        #             repairB = 1
    
    
        # Tsekataan, onko TotVol=0 -rivejÃ¤
        
        a = 0
        zeroRow = next((m for m, x in enumerate(sdata0['TotVol'][a:]) if x==0), None)
        
        if zeroRow is not None:
            
            zeroRow = a + zeroRow
                    
            while (zeroRow is not None) & (zeroRow<len(sdata0)-1): # Ei tarkisteta viimeistÃ¤ riviÃ¤ tÃ¤ssÃ¤, koska tuottaa ongelmia. Tsekataan alempana.
                
                new_vol = tulotdata['volume'+str(int(sdata0['AIKA'][zeroRow]))][0]
                
                if new_vol > 5: # PÃ¤Ã¤tehakkuun tapauksessa 0-rivi on Tulot-tiedostossa vasta kohdassa aika+1
                    new_vol = tulotdata['volume'+str(int(sdata0['AIKA'][zeroRow+1]))][0]
                
                sdata0.loc[zeroRow,'TotVol'] = new_vol
                
                # Haetaan seuraava 0-rivi:
                    
                a = zeroRow + 1    
                zeroRow = next((m for m, x in enumerate(sdata0['TotVol'][a:]) if x==0), None)
                
                if zeroRow is not None:
                    zeroRow = a + zeroRow
                else:
                    break
    
    
    
        # Jos 0-biomassoja, tarkistetaan myÃ¶s 0-iÃ¤t ja korjataan ikÃ¤ vastaamaan AIKA-jaksoja
    
        a = 0
        while a < len(sdata0):
            zeroRow = next((m for m, x in enumerate(sdata0['Age'][a:]) if x==0), None)
            if zeroRow is not None:
                zeroRow = zeroRow+a
                if sdata0['Neulasetkg'][zeroRow]==0: # Muokataan vain, jos neulasissa 0-rivejÃ¤
                    nonZeroRow = next((m for m, x in enumerate(sdata0['Age'][zeroRow:]) if x!=0), None)
                    if nonZeroRow is not None: # Jos 0-ikÃ¤ viimeisellÃ¤ rivillÃ¤, ei tehdÃ¤ tarkistuksia/korjauksia
                        nonZeroRow = nonZeroRow + zeroRow
                        if sdata0['AIKA'][nonZeroRow]-sdata0['AIKA'][zeroRow] != sdata0['Age'][nonZeroRow]-sdata0['Age'][zeroRow]:
    
                            sdata0.loc[zeroRow, 'Age'] = int(sdata0['Age'][nonZeroRow] - (sdata0['AIKA'][nonZeroRow]-sdata0['AIKA'][zeroRow]))
                            with open(outfile, "a") as myfile:
                                myfile.write('\n0-age corrected to: ' + str(sdata0.loc[zeroRow, 'Age']) + ', i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
    
                    else:
                        with open(outfile, "a") as myfile:
                            myfile.write('\n0-age detected at last row!! i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
    
                a = zeroRow + 1
            else:
                a = len(sdata0)+1     
                
                
    
    
        # Tsekataan TotVol ja DomH nollat ja muokataan
        # TÃ¤hÃ¤n tarvii ehkÃ¤ suuremman alarajan?
        # PitÃ¤Ã¤ lisÃ¤tÃ¤ myÃ¶s tyhjien tsekkaus!
        
        for k in range(0,len(sdata0)):
            
            if (sdata0.loc[k,'TotVol']<0.5) | (np.isnan(sdata0.loc[k,'TotVol']) == True): #voi olla 0.01 tilavuuksiakin
                
                if k < len(sdata0)-1:
                    sdata0.loc[k,'TotVol'] = max(0.5, 0.1*sdata0.loc[k+1,'TotVol']) # LisÃ¤tÃ¤Ã¤n totvol nollien tilalle muu pieni luku (0.001 oli liian pieni, kasvu ei lÃ¤htenyt kÃ¤yntiin, joten kasvatettu 0.2:een) paitsi et tulee 0.5
                else:
                    sdata0.loc[k,'TotVol'] = 0.5
                    
            if (sdata0.loc[k,'DomH']<0.3) | (np.isnan(sdata0.loc[k,'DomH']) == True):
                 
                 if k < len(sdata0)-1:
                     sdata0.loc[k,'DomH'] = max(0.3, 0.3*sdata0.loc[k+1,'DomH']) # LisÃ¤tÃ¤Ã¤n domh nollien tilalle muu pieni luku
                 else:
                     sdata0.loc[k,'DomH'] = 0.3       
       
        
        # Adding values to "missing" biomass components with similar increase-% per year as TotVol
        # Neulaset viimeisenÃ¤, jotta sitÃ¤ voidaan kÃ¤yttÃ¤Ã¤ arvioimaan, tarviiko biomassoja tÃ¤ydentÃ¤Ã¤.
        # Poistettu tÃ¤stÃ¤ listasta 'tuotosKuolleet', koska se on kumulatiivinen..
        
        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg', \
                     'TukPL1', 'TukPL2', 'TukPL3', 'TukPL4',\
                         'KuiPL1', 'KuiPL2', 'KuiPL3', 'KuiPL4',\
                         'HukPL1', 'HukPL2', 'HukPL3', 'HukPL4',\
                             'tuotosTilavuus',\
                                 'Runkopuukg', 'Hukkapuukg',\
                                     'MeanDbh','medianH', 'BA', 'tuotos']  
        #comps = ['DomH','MeanDbh','medianH',\
        #          'BA','Kannot_juuretkg',\
        #                  'tuotosTilavuus',\
        #                      'TukPL1', 'TukPL2','TukPL3',\
        #                          'KuiPL1', 'KuiPL2','KuiPL3',\
        #                              'HukPL1', 'HukPL2','HukPL3',\
        #                                  'VolPL1', 'VolPL2', 'VolPL3']  

        #biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg'] #tama lisatty tahan elokuu 23 jotta namakin tsekataan

        # if repairB==1: # Jos lisÃ¤tty 1-rivi  kopioimalla, korjataan sen komponentit ensin
        #     for j in range(0,len(biom_comp)):
        #         comp = biom_comp[j]
        #         sdata0.loc[1, comp] = sdata0['TotVol'][1]/sdata0['TotVol'][2]*sdata0.loc[2, comp]    
            
        # if repair==1: # Jos lisÃ¤tty 0-rivi kopioimalla, korjataan sen komponentit ensin
        #     for j in range(0,len(biom_comp)):
        #         comp = biom_comp[j]
        #         sdata0.loc[0, comp] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, comp]
                
                            
    
                
        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg']
    
        for j in range(0,len(biom_comp)):
        
            comp = biom_comp[j]
            k=0
    
            while k < len(sdata0)-1: 
                
                zeroRow = next((m for m, x in enumerate(sdata0[comp][k:]) if x==0), None)
    
                if zeroRow is not None: # If zeroRow exists
                    
                    zeroRow = zeroRow + k
                    
                    # if sdata0['Neulasetkg'][zeroRow]==0: # Muiden komponenttien numeroita muokatan vain, jos neulasissa on nollaa. TÃ¤mÃ¤ poistettu, koska nÃ¤issÃ¤ komponenteissa voi olla nollaa, vaikka neulasmassa ei olisikaan nollaa..
                    nonZeroRow = next((m for m, x in enumerate(sdata0[comp][zeroRow:]) if x!=0), None)
                    
                    if nonZeroRow is not None:
                        
                        nonZeroRow = nonZeroRow + zeroRow
                        sdata0.loc[zeroRow, comp] = sdata0['TotVol'][zeroRow]/sdata0['TotVol'][nonZeroRow]*sdata0.loc[nonZeroRow, comp] # tää vaikuttaa siihen myös et 0 riville tulee 0.5 vol ja seuraava on jotain
                    
                        k=k+1
                        
                    else:
                        
                        k=k+1                          
    
                else:
                    k=k+1
                            
    
                        
                        # else:
                            
                            # Jotta ei toistu jokaisella komponentilla erikseen, 
                            # nÃ¤mÃ¤ korjataan vasta lopuksi: 
                            # kopioidaan loppuun rivit, jotka vastaavat pÃ¤Ã¤hakkuun jÃ¤lkeistÃ¤ kehitystÃ¤
                            
                            # EtsitÃ¤Ã¤n toisen pÃ¤Ã¤hakkuun jÃ¤lkeinen rivi
            
    
             
                
             
    
    
    
        # TehdÃ¤Ã¤n sarake 'treesp'=pÃ¤Ã¤puulaji volumen perusteella
        
        sdata0['treesp'] = pd.Series([], dtype='int64') #        sdata0['treesp'] = pd.Series([])

        
        for t in range(0,len(sdata0)):
            
            
            volPine = sdata0['VolPL1'][t]
            volSpruce = sdata0['VolPL2'][t]
            volBirch = sdata0['VolPL3'][t] + sdata0['VolPL4'][t]
        
            spe = 3
        
            if volPine >= max(volSpruce, volBirch):
                spe = 1
            
            if volSpruce > max(volPine, volBirch):
                spe = 2    
    
            sdata0.loc[t, 'treesp']=spe
    
                            
                            
        # Tsekataan, lÃ¶ytyykÃ¶ viimeisen pÃ¤Ã¤tehakkuun jÃ¤lkeen pelkkiÃ¤ 0-rivejÃ¤. Jos nÃ¤in,
        # kopioidaan aiemmat pÃ¤Ã¤tehakkuun jÃ¤lkeiset rivit loppuun.
        # Korjataan sen jÃ¤lkeen myÃ¶s AIKA-sarake
    
        if len(paateh)>0:
            paate1 =  next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>paateh[0]), None) 
            nonZeroRow = next((m for m, x in enumerate(sdata0['Neulasetkg'][paate1:]) if x!=0), None)
            
            if len(paateh)>1: # Jos useampi kuin yksi pÃ¤Ã¤tehakkuu
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>=paateh[1]), None) # tÃ¤ssÃ¤ aiemmin x>paateh[1] !!
                nonZeroRow = next((m for m, x in enumerate(sdata0['Neulasetkg'][paaterow:]) if x!=0), None)
            
                if nonZeroRow is None:
                    
                    paate2aika = sdata0['AIKA'][paaterow]
                    
                    b = 0
                    
                    for p in range(paaterow,len(sdata0)):
                        print(p, paate2aika, "paatehakkuun jalkeiset nollarivit")
    
                        currow = paate1 + p - paaterow -1 # LisÃ¤tty tÃ¤hÃ¤n -1, koska rivillÃ¤ 617 lisÃ¤tty x>=paateh[1]
    
                        if currow<b: # LisÃ¤ys 17.10: Joissakin tapauksissa yllÃ¤ toimii -1, toisissa 0, joten muokataan currow tarvittaessa...
                            currow=b
                            
                        sdata0.loc[p]=sdata0.loc[currow]
                        
                        if b==0:
                            sdata0.loc[p, 'AIKA'] = int(paate2aika)
                            sdata0.loc[p, 'tuotosKuolleet'] = sdata0.loc[p-1, 'tuotosKuolleet']
                            b=b+1
                        else:
                            sdata0.loc[p, 'AIKA'] = int(paate2aika + sdata0.loc[currow, 'AIKA'] - sdata0.loc[currow-b, 'AIKA'])
                            sdata0.loc[p, 'tuotosKuolleet'] = sdata0.loc[p-1, 'tuotosKuolleet']
                            b=b+1
                
            else:
                # Jos vain yksi pÃ¤Ã¤tehakkuu, ei voida kopioida rivejÃ¤ aiemman pÃ¤Ã¤tehakkuun jÃ¤lkeen..
                print('MitÃ¤s tÃ¤hÃ¤n? kuvio:', str(kuviolista[i]),paateh)   
                if nonZeroRow is None:
                    with open(outfile, "a") as myfile:
                        myfile.write('\nPÃ¤Ã¤tehakkuun jÃ¤lkeisiÃ¤ 0-rivejÃ¤ ei saada korjattua! i=' + str(i)+ ' kuvio:'+ str(kuviolista[i])+ ' paateh: '+str(paateh)) #Ã¤Ã¤ 
    
        
        # Tsekataan, lÃ¶ytyykÃ¶ viimeiseltÃ¤ riviltÃ¤ (aika=100), nollariviÃ¤.
        # Jos lÃ¶ytyy, poistetaan.
        
    
        s100 = sdata0[(sdata0['AIKA']==100) & (sdata0['Age']==0)] # Vai vol=0?
        
        if len(s100)>0:
            sdata0 = sdata0[sdata0['AIKA']!=100]
            
        # -------------------------------------------------------------------
        # LisÃ¤ys 13.9.2021
        # Tarkistetaan, olisiko ensimmÃ¤inen Motti-tiedosto 1-rivin mittainen:
        
        vuodet_p = vuodet[vuodet>0]
        
        if (len(vuodet_p)>0):
            
            sdata = sdata0[sdata0['AIKA'] < vuodet_p[0]]
            
            if len(sdata)==1:
                #sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
                nolla2_df=pd.DataFrame([sdata0.iloc[0]])
                sdata0 = pd.concat((sdata0,nolla2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = vuodet_p[0]-1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                sdata0.loc[1,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][1]))][0]
                sdata0.loc[1,'Age'] = sdata0['Age'][0]+(sdata0['AIKA'][1]-sdata0['AIKA'][0])
    
                for c in range(0,len(comps)):
                    if sdata0.loc[0,comps[c]]>0:
                        # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
                        # sitÃ¤ voisi supistaa..??
                        sdata0.loc[1,comps[c]] = sdata0['TotVol'][1]/sdata0['TotVol'][0]*sdata0.loc[0, comps[c]]
    
                for c in range(0,len(biom_comp)):
                    if sdata0.loc[0,biom_comp[c]]>0:
                        # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
                        # sitÃ¤ voisi supistaa..??
                        sdata0.loc[1,biom_comp[c]] = sdata0['TotVol'][1]/sdata0['TotVol'][0]*sdata0.loc[0, biom_comp[c]]
    
    
        # -------------------------------------------------------------------
        # LisÃ¤ys 9.3.2022
        # Tarkistetaan, olisiko *toinen* Motti-tiedosto 1-rivin mittainen:
        
        #if (len(vuodet_p)>0):
        #    
        #    if len(vuodet_p)==1:
        #        end_year=101
        #    else:
        #        end_year = vuodet_p[1]
        #        
        #    sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[0]) & (sdata0['AIKA'] < end_year)] # Huom. tÃ¤ssÃ¤ ei saa resetoida indeksiÃ¤, jotta jatko toimii
        #    
        #    if len(sdata)==1:
        #        sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # LisÃ¤tÃ¤Ã¤n rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
        #        sdata0['AIKA'][len(sdata0)-1] = end_year - 1
        # 
        #         sdata0 = sdata0.sort_values(by=['AIKA'])
        #         sdata0 = sdata0.reset_index(drop=True)
        #         
        #        sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
        #        
        #        sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
        #        
        #        if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa riviÃ¤
        #            sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
        #
        #            for c in range(0,len(comps)):
        #                if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
        #                    # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
        #                    # sitÃ¤ voisi supistaa..??
        #                    sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
        # 
        #            for c in range(0,len(biom_comp)):
        #                if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
        #                    # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
        #                    # sitÃ¤ voisi supistaa..??
        #                    sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
        
        
        # -------------------------------------------------------------------
        # Update 6.3.2023: checking all 1-row-motti-files except the last one
        
        for v in range(1, len(vuodet_p)):
            
            if (len(vuodet_p)>0):
                
                if len(vuodet_p)==1:
                    end_year=101
                else:
                    end_year = vuodet_p[v]                
                
                sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[v-1]) & (sdata0['AIKA'] < end_year)] # Huom. tÃ¤ssÃ¤ ei saa resetoida indeksiÃ¤ jotta jatko toimii
                
                if len(sdata)==1:
                    #sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # LisÃ¤tÃ¤Ã¤n rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
                    sdataind_df=pd.DataFrame([sdata0.iloc[sdata.index[0]]])
                    sdata0 = pd.concat((sdata0,sdataind_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                    sdata0['AIKA'][len(sdata0)-1] = end_year - 1
                    
                    sdata0 = sdata0.sort_values(by=['AIKA'])
                    sdata0 = sdata0.reset_index(drop=True)
                    
                    sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
                    
                    sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
                    
                    if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa rivia
                        sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
                        
                        for c in range(0,len(comps)):
                            if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
                                # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                                # sita voisi supistaa..??
                                sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
                        
                        for c in range(0,len(biom_comp)):
                            if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
                                # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                                # sita voisi supistaa..??
                                sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
        
        # Update 6.3.2023: Checking if the last motti-file would be a 1-row-file
        
        if (len(vuodet_p)>1):
            
            sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[len(vuodet_p)-1])] # Huom. tassa ei saa resetoida indeksia, jotta jatko toimii
            end_year = 101
            
            if (len(sdata)==1) & (sdata0['AIKA'][len(sdata0)-1]<100):
                
                #sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # Lisataan rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
                sdataind2_df=pd.DataFrame([sdata0.iloc[sdata.index[0]]])
                sdata0 = pd.concat((sdata0, sdataind2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                    
                sdata0['AIKA'][len(sdata0)-1] = end_year - 1
                
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
                
                sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
                
                if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa rivia
                    sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
                    
                    for c in range(0,len(comps)):
                        if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
                            # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                            # sita voisi supistaa..??
                            sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
                    
                    for c in range(0,len(biom_comp)):
                        if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
                            # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                            # sita voisi supistaa..??
                            sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
                            
                
        
        # -------------------------------------------------------------------
        # LisÃ¤ys 9.3.2022
        # Tarkistetaan, tuleeko perÃ¤kkÃ¤isten harvennusten takia liian paljon
        # vÃ¤heneviÃ¤ tilavuuksia 
        
        comps_notk = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg', \
                             'tuotosTilavuus','MeanDbh','medianH', 'BA', 'tuotos', 'TotVol']  
        
        
        if pera >= 0:
            
            vuodet_p = vuodet[vuodet>=0] # LisÃ¤ys 17.10.2022: vuodet>=0 (jos tÃ¤ssÃ¤ olisi vuodet>0, vuodet_p saattaisi jÃ¤Ã¤dÃ¤ tyhjÃ¤ksi -> virhe seuraavalla rivillÃ¤)
            
            sdata = sdata0[sdata0['AIKA'] < vuodet_p[0]]
            
            if (len(sdata)<=2) & (len(sdata)>0): # LisÃ¤ys 17.10.2022: len(sdata)>0, koska jos sdata on tyhjÃ¤ -> error
                
                if (sdata['TotVol'][0] - sdata['TotVol'][1]) > 0:
                    
                    sdata0.loc[0]=sdata0.loc[1]
                    sdata0.loc[0, 'AIKA'] = sdata0['AIKA'][1] -1
                    sdata0.loc[0, 'Age'] = sdata0['Age'][1] -1
                    
                    
                    for c in range(0,len(comps_notk)):
                        
                        comp = comps_notk[c]
                        sdata0.loc[0, comp] = sdata0.loc[0, comp] * 0.97 # Arvio, jotta 
                                            
    
    
        
        
        # Muutetaan Age=-1 --> Age=0, nÃ¤itÃ¤ lÃ¶ytyi suoraan Motti-datasta
        
        sdata0['Age'] = np.clip(sdata0['Age'], 0, np.nan) 
        
        # -------------------------------------------------------------------
        # Lisays 28.3.2023:
        # Tarkistetaan, loytyyko iasta perakkaisia nollia
        
        a = 0
        while a < len(sdata0):
            zeroRow = next((m for m, x in enumerate(sdata0['Age'][a:]) if x==0), None)
            
            if zeroRow is not None:
                zeroRow = zeroRow+a
                try:
                    if sdata0['Age'][zeroRow+1]==0:
                        sdata0.loc[zeroRow, 'Age'] = 1
                        sdata0.loc[zeroRow+1, 'Age'] = 1 + int(sdata0['AIKA'][zeroRow+1]-sdata0['AIKA'][zeroRow])
                        with open(outfile, "a") as myfile:
                            myfile.write('\nDouble zero ages!! i=' + str(i))
                        
                    a = zeroRow + 1
                except:
                    a = zeroRow + 1
                        
            else:
                a = len(sdata0)+1     
        
        
        
        
        # Varsinaisten Motti-tiedostojen teko alkaa tÃ¤stÃ¤
        
        if len(sdata0)>0:
            
            if len(vuodet)>0:
                
                for j in range(0,n_mottifiles):
                
                    if j == 0:
                        sdata = sdata0[sdata0['AIKA'] < vuodet[j]]
                        if len(sdata)==0:
                            continue
                    
                    if (j>0) & (j<n_mottifiles-1):
                        sdata = sdata0[(sdata0['AIKA'] >= vuodet[j-1]) & (sdata0['AIKA'] < vuodet[j])]
                        sdata = sdata.reset_index(drop=True)
                        if len(sdata)==0:
                            continue
                    
                    if j == n_mottifiles-1:
                        sdata = sdata0[sdata0['AIKA'] >= vuodet[j-1]]
                        sdata = sdata.reset_index(drop=True)
                        if len(sdata)==0:
                            continue
                    
                    kasvatus = pd.Series(np.ones(len(sdata)))
                    vuosi = sdata['AIKA']
                    ika = sdata['Age']
                    N = sdata['N']
                    PPA = sdata['BA']
                    Hg = sdata['medianH']
                    Dg = sdata['MeanDbh']
                    Hdom = sdata['DomH']
                    Tilavuus = sdata['TotVol']
                    Tukki = sdata['TukPL1']+sdata['TukPL2']+sdata['TukPL3']+sdata['TukPL4']
                    Kuitu = sdata['KuiPL1']+sdata['KuiPL2']+sdata['KuiPL3']+sdata['KuiPL4']
                    Hukka = sdata['HukPL1']+sdata['HukPL2']+sdata['HukPL3']+sdata['HukPL4']
                    Tuotos = sdata['tuotosKuolleet'] + sdata['tuotosTilavuus']
                    Kuolleisuus = Tuotos*0. # Ei kÃ¤ytetÃ¤ Susissa
                    runko_aines = sdata['Runkopuukg']/1000. # YksikkÃ¶ tn
                    runko_hukka = sdata['Hukkapuukg']/1000.
                    elavat_oksat = sdata['Oksatkg']*0.8/1000. # Oksatkg = elÃ¤vÃ¤t ja kuolleet oksat yhteensÃ¤
                    kuolleet_oksat = sdata['Oksatkg']*0.2/1000.
                    lehdet = sdata['Neulasetkg']/1000.
                    kannot = sdata['kannot']/1000.
                    juuret = sdata['>2mm_juuretkg']/1000.
                    hienojuuret = sdata['h_juuretkg']/1000.
                    
                    # Tsekkaa, ettÃ¤ tÃ¤mÃ¤ on ok!!!
                    for t in range(0,len(sdata)):
                        if t==0:
                            Tuotos[0]=sdata['tuotosTilavuus'][0]
                        
                        else:
                            Tuotos[t] = sdata['tuotosTilavuus'][t] + (sdata['tuotosKuolleet'][t] - sdata['tuotosKuolleet'][t-1])
                    
                    # ---tÃ¤hÃ¤n asti Ã¤
                    
                    data = pd.DataFrame(np.transpose(np.array([kasvatus, vuosi, ika, N, PPA, Hg, Dg, Hdom, Tilavuus, Tukki, Kuitu, \
                                         Hukka, Tuotos, Kuolleisuus, runko_aines, runko_hukka, \
                                         elavat_oksat, kuolleet_oksat, lehdet, kannot, juuret, hienojuuret])), \
                                         columns=['Kasvatus', 'Vuosi', 'Ikä', 'N', 'PPA', 'Hg', 'Dg', 'Hdom', 'Tilavuus', 'Tukki', 'Kuitu',\
                                                  'Hukka', 'Tuotos', 'Kuolleisuus', 'runko(aines)', 'runko(hukka)', \
                                                  'elävät oksat', 'kuolleet oksat', 'lehdet', 'Kannot', 'Juuret >2mm', 'Hienojuuret'])
                    
                    standInfo = kuviot[kuviot['KUVIO']==kuviolista[i]]
                    standInfo = standInfo.reset_index(drop=True)
                    
                    """ TSEKKAA, ETTÃ„ PUULAJI TULEE OIKEIN! """
                    
                    
                    
                    if (project == 'life') | (project == 'hiilipolku'):
                        spe = sdata['treesp'][len(sdata)-1]
                    
                    else: 
                        if (project =='vmi') | (project=='suo'):
                            
                            spe = int(sdata['treesp'][len(sdata)-1]) # viimeisen rivin puulaji
                            # treesp=kdata['Selite'][0][28:33]
                            # spe = 3
                            
                            # if treesp=='MÃ¤nty':
                            #     spe = 1
                            # if treesp=='Kuusi':
                            #     spe = 2       
                                
                            
                        else:
                            
                            treesp = pdata2['Selite'][0][28:33]
                            
                            spe = 3
                            
                            if treesp=='Mänty': #onk tÃ¤ a va Ã¤#
                                spe = 1
                            if treesp=='Kuusi':
                                spe = 2 
                    
                    
                    # print(kuviolista[i], ': puulaji', spe)
                    with open(outfile, "a") as myfile:
                        myfile.write('\n' + str(kuviolista[i]) +  ': puulaji ' + str(spe)+ ' n'+str(j))
                    
                    
                    fout = mottifolder + '/' + str(kuviolista[i]) + '_n' + str(j) + '.xls'
                    
                    cols = pd.DataFrame([['','','','',spe,'','','','','','','','','']],columns=['Kasvatus',	'Vuosi','id Harvennus',	'Harvennus',	'id Puulaji',	'Puulaji',	'Tukki[m³/ha]',	'Pikkutukki[m³/ha]',	'Kuitu[m³/ha]',	'Energiapuu, runko(aines)[m³/ha]',	'Energiapuu, runko(hukka)[m³/ha]',	'Energiapuu, oksat(e)[m³/ha]',	'Energiapuu, oksat(k)[m³/ha]',	'Energiapuu, kannot ja juuret[m³/ha]']) #Â in [mÂ³/ha]
                    
                    writer = pd.ExcelWriter(fout, engine = 'xlsxwriter')
                    data.to_excel(writer, index=False, sheet_name = 'Puustotunnukset')
                    cols.to_excel(writer, index=False, sheet_name = r'Kertymät') #Ã¤
                    #writer.save()
                    writer.close()
            else:
                # print(str(kuviolista[i]) + ': poistumat puuttuvat! i=' + str(i))
                with open(outfile, "a") as myfile:
                    myfile.write('\n' + str(kuviolista[i]) + ': poistumat puuttuvat! i=' + str(i))
                
                sdata = sdata0.copy()
                
                kasvatus = pd.Series(np.ones(len(sdata)))
                vuosi = sdata['AIKA']
                ika = sdata['Age']
                N = sdata['N']
                PPA = sdata['BA']
                Hg = sdata['medianH']
                Dg = sdata['MeanDbh']
                Hdom = sdata['DomH']
                Tilavuus = sdata['TotVol']
                Tukki = sdata['TukPL1']+sdata['TukPL2']+sdata['TukPL3']+sdata['TukPL4']
                Kuitu = sdata['KuiPL1']+sdata['KuiPL2']+sdata['KuiPL3']+sdata['KuiPL4']
                Hukka = sdata['HukPL1']+sdata['HukPL2']+sdata['HukPL3']+sdata['HukPL4']
                Tuotos = sdata['tuotos']
                Kuolleisuus = Tuotos*0. # Ei kÃ¤ytetÃ¤ Susissa
                runko_aines = sdata['Runkopuukg']/1000. # YksikkÃ¶ tn
                runko_hukka = sdata['Hukkapuukg']/1000.
                elavat_oksat = sdata['Oksatkg']*0.8/1000. # Oksatkg = elÃ¤vÃ¤t ja kuolleet oksat yhteensÃ¤
                kuolleet_oksat = sdata['Oksatkg']*0.2/1000.
                lehdet = sdata['Neulasetkg']/1000.
                kannot = sdata['kannot']/1000.
                juuret = sdata['>2mm_juuretkg']/1000.
                hienojuuret = sdata['h_juuretkg']/1000.
                
                data = pd.DataFrame(np.transpose(np.array([kasvatus, vuosi, ika, N, PPA, Hg, Dg, Hdom, Tilavuus, Tukki, Kuitu, \
                                     Hukka, Tuotos, Kuolleisuus, runko_aines, runko_hukka, \
                                     elavat_oksat, kuolleet_oksat, lehdet, kannot, juuret, hienojuuret])), \
                                     columns=['Kasvatus', 'Vuosi', 'Ikä', 'N', 'PPA', 'Hg', 'Dg', 'Hdom', 'Tilavuus', 'Tukki', 'Kuitu',\
                                              'Hukka', 'Tuotos', 'Kuolleisuus', 'runko(aines)', 'runko(hukka)', \
                                              'elävät oksat', 'kuolleet oksat', 'lehdet', 'Kannot', 'Juuret >2mm', 'Hienojuuret'])
                
                standInfo = kuviot[kuviot['KUVIO']==kuviolista[i]]
                standInfo = standInfo.reset_index(drop=True)
                
                
                """ TSEKKAA, ETTÃ„ PUULAJI TULEE OIKEIN! """
                if (project == 'life') | (project == 'hiilipolku'):
                    spe = sdata['treesp'][len(sdata)-1]
                
                else:
                    
                    if project =='vmi':
                        
                        # spe = int(sdata['treesp'][len(sdata)-1]) # viimeisen rivin puulaji
                       
                        treesp=kdata['Selite'][0][28:33]
                        spe = 3
                        
                        if treesp=='Mänty': #Ã¤
                            spe = 1
                        if treesp=='Kuusi':
                            spe = 2    
                    else:
                        treesp = pdata2['Selite'][0][28:33]
                        
                        spe = 3
                        
                        if treesp=='Mänty':
                            spe = 1
                        if treesp=='Kuusi':
                            spe = 2 
                    
                # print(kuviolista[i], ': puulaji ', spe)
                with open(outfile, "a") as myfile:
                    myfile.write('\n' + str(kuviolista[i]) +  ':' + str(spe)+ ' n0')
                
                
                
                fout = mottifolder + '/' + str(kuviolista[i]) + '_n0.xls'
                
                cols = pd.DataFrame([['','','','',spe,'','','','','','','','','']],columns=['Kasvatus',	'Vuosi',	'id Harvennus',	'Harvennus',	'id Puulaji',	'Puulaji',	'Tukki[m³/ha]',	'Pikkutukki[m³/ha]',	'Kuitu[m³/ha]',	'Energiapuu, runko(aines)[m³/ha]',	'Energiapuu, runko(hukka)[m³/ha]',	'Energiapuu, oksat(e)[m³/ha]',	'Energiapuu, oksat(k)[m³/ha]',	'Energiapuu, kannot ja juuret[m³/ha]'])
                
                writer = pd.ExcelWriter(fout, engine = 'xlsxwriter')
                data.to_excel(writer, index=False, sheet_name = 'Puustotunnukset')
                cols.to_excel(writer, index=False, sheet_name = r'Kertymät')
                #writer.save()
                writer.close()
                
def create_motti_files_silvi_lp4(project, scen, maku=None):       #lp4 skaalaukseen otettu lp4 mukaan  # Huom! ohje>1
    '''
    Modified: January 2023
        
    Create Motti input files for Susi simulations from "Skene-Motti" results.
    
    Works technically, but contains many small issues that should be fixed.
    
    project = 'life' (Hydrology Life -project) or 'suo' (SUO-project)
    project-parameter affects only main tree species (and column names)!
    
    
    Modified 8.3.2022:
        - korjattu pÃ¤Ã¤tehakkuuseen liittyviÃ¤ ongelmia; hdom pitÃ¤isi nyt saada
        jÃ¤rkeviÃ¤ arvoja, biomassojen ei pitÃ¤isi jÃ¤Ã¤dÃ¤ pÃ¤Ã¤tehakkuun jÃ¤lkeen nollaksi
        (ainakaan toisen pÃ¤Ã¤tehakkuun jÃ¤lkeen)
        
    Issues to fix:
        - Jos pÃ¤Ã¤tehakkuu tulee vuonna 2, pitÃ¤Ã¤ vuosille 0 ja 1 hakea muut komponentit
        totvol perusteella: hdom saa tÃ¤ssÃ¤ liian suuria arvoja. Tuotos jÃ¤Ã¤ nollaksi?!
        - See comments below...
        
    
    '''
    
    inputmainfol=r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan_lp4/'
    mottifolder = inputmainfol+scen[:-6]+'/'+scen+'/'
    if not os.path.exists(inputmainfol):
        os.mkdir(inputmainfol)
    if not os.path.exists(inputmainfol+scen[:-6]+'/'):
        os.mkdir(inputmainfol+scen[:-6]+'/')
    if not os.path.exists(mottifolder):
        os.mkdir(mottifolder)
            
    outfile = mottifolder + '/out.txt'
    
    project = 'hiilipolku'
    
    kuviot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', encoding='latin1', sep=';')
    puustot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_puustot.csv', encoding='latin1', sep=';')
    poistumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_poistumat.csv', encoding='latin1', sep=';')
    tapahtumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_tapahtumat.csv', encoding='latin1', sep=';', index_col=False)
    tulot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kasvut.csv', encoding='latin1', sep=';')
    
    
    if (project == 'life') | (project =='vmi') | (project =='hiilipolku'):
        kuviot = kuviot.rename(columns={'kuvio':'KUVIO'})
        puustot = puustot.rename(columns={'kuvio':'KUVIO',
                                          '_2mm_juuretkg':'>2mm_juuretkg',})
        poistumat = poistumat.rename(columns={'kuvio':'KUVIO'})
        tapahtumat = tapahtumat.rename(columns={'kuvio':'KUVIO',
                                      '_2003': '2003',
                                      '_2004': '2004',
                                      '_2005': '2005'})
        tulot = tulot.rename(columns={'kuvio':'KUVIO'})
    
    
    kuviolista = list(set(kuviot['KUVIO'][kuviot['kuivatustilanne']>=3])) #tassa valitaan 3 Ojitettu kangas ja muut?DRAINAGESTATE	, ojitettu kangas jos pÃ¤Ã¤ryhmÃ¤1 ja alaryhmÃ¤ 3 ja kuivatustilanne 5?? 1 Ojittamaton kangas, 2 Soistunut kangas,3 Ojitettu kangas,	6 luonnontilainen suo,	7 Ojikko, 8	Muuttuma,9	Turvekangas
    
    
    for i in range(0,len(kuviolista)):
    
        print(i, 'stand: ', kuviolista[i])
        
        # alta poistettu ohjeeseen ja skenaarioon viittaavat rajaukset,
        # oletetaan, ettÃ¤ hiilipolussa vain yksi ohje per kuvio
        
        sdata0 = puustot[(puustot['KUVIO']==kuviolista[i])]
        pdata = poistumat[(poistumat['KUVIO']==kuviolista[i]) & (poistumat['HARV']>1)] #mita harv koodit meinaa?
        pdata2 = poistumat[(poistumat['KUVIO']==kuviolista[i])]
        tdata = tapahtumat[(tapahtumat['KUVIO']==kuviolista[i])]
        tulotdata = tulot[(tulot['KUVIO']==kuviolista[i])]
        kdata = kuviot[kuviot['KUVIO']==kuviolista[i]]
        
        sdata0 = sdata0.reset_index(drop=True)
        pdata = pdata.reset_index(drop=True)
        pdata2 = pdata2.reset_index(drop=True)
        tdata = tdata.reset_index(drop=True)
        tulotdata = tulotdata.reset_index(drop=True)
        kdata = kdata.reset_index(drop=True)
        
        tdata = tdata[(tdata['2003']>-1) | (tdata['2004']>-1) | (tdata['2005']>-1)]
        
        ensih = tdata[tdata['2003']>-1]['2003'].values
        harv = tdata[tdata['2004']>-1]['2004'].values
        paateh = tdata[tdata['2005']>-1]['2005'].values
        #voisi kayda tsekkaamassa onko taimikonhoito tai varhaisperkaus niin sitten ei muuttaisi puulajia tai ottaisi sen muualta et sailyy kuusena tai mantyna? 
        #tdic={'_2001':'lannoitus','_2002':'kunnostusojitus','_2003':'ensiharvennus','_2004':'harvennus','_2005':'paatehakkuu','_2006':'varhaiperkaus','_2007':'taimikonhoito','_2008':'viljely','_2009':'luontainen uudistuminen'}
        # taimh = tdata[tdata['_2007']>-1]['_2007'].values
        # Pudotetaan viimeisen vuoden tapahtumat pois, koska sekoittavat
        # myÃ¶hemmin tulevan koodin
        ensih = ensih[ensih<100]
        harv = harv[harv<100]
        paateh = paateh[paateh<100]
        
        vuodet = np.unique(np.sort(np.concatenate((ensih, harv, paateh))))
        
        ap = 0
        
        pera = -1
        
        # Tsekataan, onko tapahtumissa perÃ¤kkÃ¤isiÃ¤ vuosia ja poistetaan niistÃ¤ jÃ¤lkimmÃ¤inen
        for v in range(0,len(vuodet)-1):
            if vuodet[v+1]-vuodet[v]==1:
                
                ensih = np.delete(ensih, np.where(ensih == vuodet[v+1]))
                harv = np.delete(harv, np.where(harv == vuodet[v+1]))
                paateh = np.delete(paateh, np.where(paateh == vuodet[v+1]))
                
                vuodet = np.delete(vuodet, [v+1], None)
                
                pera = v
                
                ap = ap + 1
                
            if v+1>=len(vuodet)-1:
                break
        
        if ap>1:
            with open(outfile, "a") as myfile:
                myfile.write('\n\nKaksi kertaa perakkaiset vuodet! stand=' + str(kuviolista[i]) + ', i=' + str(i))
    
        n_mottifiles = len(vuodet) + 1
    
    
        if len(sdata0)==0: # skipataan kuvio, jos dataa ei ole
            continue
    
        if len(sdata0)>21:
            print('sdata0: liikaa rivejÃ¤?!')
        
    
                                    
        # Luodaan pÃ¤Ã¤tehakkuun jÃ¤lkeiset rivit ja korjataan rivin kopioimisesta aiheutuneet
        # luvut
        
        # muokattu 6.8.2021: poistettu tuotos, tuotosKuolleet, tuotosTilavuus comps:sta
        comps = ['DomH','MeanDbh','medianH',\
                          'BA','Kannot_juuretkg',\
                              'tuotosTilavuus',\
                              'TukPL1', 'TukPL2','TukPL3','TukPL4',\
                                  'KuiPL1', 'KuiPL2','KuiPL3','KuiPL4',\
                                      'HukPL1', 'HukPL2','HukPL3','HukPL4',\
                                          'VolPL1', 'VolPL2', 'VolPL3','VolPL4',]  #lisatty 091023 PL4 ositteet mukaan

        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg'] #tama lisatty tahan elokuu 23 jotta namakin tsekataan

        for h in range(0,len(paateh)):
            
            sdata_paateh= sdata0[sdata0['AIKA']==paateh[h]].reset_index(drop=True)
            
            go = False
            
            if len(sdata_paateh)>0:
                if sdata_paateh.loc[0, 'Age']>5:
                    go = True
        
            if (len(sdata_paateh)==0) | (go==True): # Luodaan rivit vain, jos ensimmÃ¤istÃ¤ riviÃ¤ ei ole mukana
                
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>paateh[h]), None)
                
                #sdata0 = sdata0.append([sdata0.iloc[paaterow]], ignore_index=True) # Rivi menee taulukon loppuun
                paaterow_df=pd.DataFrame([sdata0.iloc[paaterow]])
                sdata0 = pd.concat((sdata0,paaterow_df),axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = paateh[h]+1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==paateh[h]+1), None)
                
                # sdata0.loc[paaterow,'Age'] = 1 # TÃ¤mÃ¤ korvattu alla olevalla
                
                if paaterow>0:
                    if sdata0.loc[paaterow-1, 'Age']<5:
                        sdata0.loc[paaterow,'Age'] = sdata0.loc[paaterow-1, 'Age']+1
                    else:
                        sdata0.loc[paaterow,'Age'] = sdata0.loc[paaterow+1, 'Age']-1
                
                sdata0.loc[paaterow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][paaterow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[paaterow,comps[c]]>0:
                        sdata0.loc[paaterow,comps[c]] = sdata0['TotVol'][paaterow]/sdata0['TotVol'][paaterow+1]*sdata0.loc[paaterow+1, comps[c]]
                        print(kuviolista[i],comps[c], "rivi 303 ongelma")

            
                if sdata0['tuotosKuolleet'][paaterow] < sdata0['TotVol'][paaterow]:
                    sdata0.loc[paaterow,'N'] = (1 + sdata0['tuotosKuolleet'][paaterow]/sdata0['TotVol'][paaterow])*sdata0['N'][paaterow]        
                    
                    if paaterow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][paaterow] > sdata0['N'][paaterow-1]:
                            if paaterow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[paaterow,'N'] = (sdata0['N'][paaterow-1] + sdata0['N'][paaterow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[paaterow,'N'] = 0.99*sdata0['N'][paaterow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.
                            
                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[paaterow,biom_comp[c]]>0: 
                        sdata0.loc[paaterow,biom_comp[c]] = sdata0['TotVol'][paaterow]/sdata0['TotVol'][paaterow+1]*sdata0.loc[paaterow+1, biom_comp[c]] ####
                
            
        # Luodaan ensiharvennuksen jÃ¤lkeiset rivit (jos vuosi ei osu kohdallaan muuten)
        # EsimerkissÃ¤ nÃ¤ytti siltÃ¤, ettÃ¤ ensiharvennusvuoden riville tulostuu harventamaton tilavuus! Siksi tÃ¤ssÃ¤ kÃ¤sitelty riviÃ¤+1
        # Toisessa esimerkissÃ¤ tulostui harvennettu tilavuus! Joten vaihtelee.
        
        for h in range(0,len(ensih)):
    
    
            ensihrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==ensih[h]+1), None)
            ensihrow2 = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>ensih[h]), None)
            
            if ensihrow is None:
                
                #sdata0 = sdata0.append([sdata0.iloc[ensihrow2]], ignore_index=True) # Rivi menee taulukon loppuun
                ensihrow2_df=pd.DataFrame([sdata0.iloc[ensihrow2]])
                sdata0 = pd.concat((sdata0,ensihrow2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = ensih[h]+1
                
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                ensihrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==ensih[h]+1), None)
                
                # Jos ikÃ¤ lÃ¶ytyy riviÃ¤ alaempaa:
                if sdata0.loc[ensihrow+1,'Age']>0:
                    sdata0.loc[ensihrow,'Age'] = sdata0.loc[ensihrow+1,'Age']-(sdata0['AIKA'][ensihrow+1]-sdata0['AIKA'][ensihrow])
                # Muussa tapauksessa ikÃ¤ lÃ¶ytyy riviÃ¤ ylempÃ¤Ã¤:
                else:
                    if ensihrow>0:
                        sdata0.loc[ensihrow,'Age'] = sdata0.loc[ensihrow-1,'Age']+(sdata0['AIKA'][ensihrow]-sdata0['AIKA'][ensihrow-1])
                    else:
                        # Jos ei lÃ¶ydy kummaltakaan puolelta:
                        with open(outfile, "a") as myfile:
                            myfile.write('\nEnsiharvennuksen jÃ¤lkeistÃ¤ ikÃ¤Ã¤ ei lÃ¶ydy! kuvio=' + str(kuviolista[i]) + ' i:'+ str(i))
                        
                    
                sdata0.loc[ensihrow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][ensihrow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[ensihrow,comps[c]]>0:
                        sdata0.loc[ensihrow,comps[c]] = sdata0['TotVol'][ensihrow]/sdata0['TotVol'][ensihrow+1]*sdata0.loc[ensihrow+1, comps[c]]
                
                # PÃ¤ivitetÃ¤Ã¤n N
                if sdata0['tuotosKuolleet'][ensihrow] < sdata0['TotVol'][ensihrow]:
                    sdata0.loc[ensihrow,'N'] = (1 + sdata0['tuotosKuolleet'][ensihrow]/sdata0['TotVol'][ensihrow])*sdata0['N'][ensihrow]        
    
                    if ensihrow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][ensihrow] > sdata0['N'][ensihrow-1]:
                            if ensihrow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[ensihrow,'N'] = (sdata0['N'][ensihrow-1] + sdata0['N'][ensihrow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[ensihrow,'N'] = 0.99*sdata0['N'][ensihrow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.

                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[ensihrow,biom_comp[c]]>0: 
                        sdata0.loc[ensihrow,biom_comp[c]] = sdata0['TotVol'][ensihrow]/sdata0['TotVol'][ensihrow+1]*sdata0.loc[ensihrow+1, biom_comp[c]] ####
                            
        
        # Luodaan harvennuksen jÃ¤lkeiset rivit (jos vuosi ei osu kohdalleen muuten)
        # Harvennusvuoden rivillÃ¤ on harventamaton tilavuus. Siksi kÃ¤sitelty riviÃ¤+1.
        
        for h in range(0,len(harv)):
    
    
            hrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==harv[h]+1), None)
            hrow2 = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>harv[h]), None)
            
            if hrow is None:
                
                #sdata0 = sdata0.append([sdata0.iloc[hrow2]], ignore_index=True) # Rivi menee taulukon loppuun
                hrow2_df=pd.DataFrame([sdata0.iloc[hrow2]])
                sdata0 = pd.concat((sdata0,hrow2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = harv[h]+1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                hrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==harv[h]+1), None)
                
                # Jos ikÃ¤ lÃ¶ytyy riviÃ¤ alaempaa:
                if sdata0.loc[hrow+1,'Age']>0:
                    sdata0.loc[hrow,'Age'] = sdata0.loc[hrow+1,'Age']-(sdata0['AIKA'][hrow+1]-sdata0['AIKA'][hrow])
                # Muussa tapauksessa ikÃ¤ lÃ¶ytyy riviÃ¤ ylempÃ¤Ã¤:
                else:
                    if hrow>0:
                        sdata0.loc[hrow,'Age'] = sdata0.loc[hrow-1,'Age']+(sdata0['AIKA'][hrow]-sdata0['AIKA'][hrow-1])
                    else:
                    # Jos ei lÃ¶ydy kummaltakaan puolelta:
                        with open(outfile, "a") as myfile:
                            myfile.write('\nHarvennuksen jÃ¤lkeistÃ¤ ikÃ¤Ã¤ ei lÃ¶ydy! i=' + str(i) + ' kuvio:'+ str(kuviolista[i]))
                
                    
                sdata0.loc[hrow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][hrow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[hrow,comps[c]]>0:
                        sdata0.loc[hrow,comps[c]] = sdata0['TotVol'][hrow]/sdata0['TotVol'][hrow+1]*sdata0.loc[hrow+1, comps[c]]

                        
                # Korjataan N:
                if sdata0['tuotosKuolleet'][hrow] < sdata0['TotVol'][hrow]:
                    sdata0.loc[hrow,'N'] = (1 + sdata0['tuotosKuolleet'][hrow]/sdata0['TotVol'][hrow])*sdata0['N'][hrow]
                    
                    if hrow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][hrow] > sdata0['N'][hrow-1]:
                            if hrow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[hrow,'N'] = (sdata0['N'][hrow-1] + sdata0['N'][hrow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[hrow,'N'] = 0.99*sdata0['N'][hrow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.
                            
                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[hrow,biom_comp[c]]>0: 
                        sdata0.loc[hrow,biom_comp[c]] = sdata0['TotVol'][hrow]/sdata0['TotVol'][hrow+1]*sdata0.loc[hrow+1, biom_comp[c]] ####
    
        # Tsekataan, ettei ole ekan rivin TotVol pyÃ¶ristyksen kanssa ongelmaa:
        
        if (sdata0['TotVol'][0]>10.) & (-0.1 < sdata0['TotVol'][1]-sdata0['TotVol'][0]<0):
            sdata0.loc[0,'TotVol'] = sdata0['TotVol'][1]
                    
          
        # Jos yllÃ¤ olevien korjausten jÃ¤lkeekin aika on alussa > 1,
        # luodaan aika=0 -rivi.
    
        
        repair = 0
        repairB = 0
        
        if (sdata0['AIKA'][0]>1) & (sdata0['Age'][0]>=sdata0['AIKA'][0]): # Jos aika=0 tai aika=1, lisÃ¤tÃ¤Ã¤n alkuun rivi, joka kopioitu ensimmÃ¤iseltÃ¤ olemassa olevalta riviltÃ¤, kÃ¤ytÃ¤nnÃ¶ssÃ¤ usein varhaisperkauksen jÃ¤lkeinen tilanne
            # print('Aika alussa > 1, i=' + str(i))
            with open(outfile, "a") as myfile:
                myfile.write('\nAika alussa > 1, i=' + str(i) + ' kuvio: '+ str(kuviolista[i]))
             
            #sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
            nolla_df=pd.DataFrame([sdata0.iloc[0]])
            sdata0 = pd.concat((sdata0,nolla_df),axis=0, ignore_index=True) # Rivi menee taulukon loppuun
            sdata0['AIKA'][len(sdata0)-1] = 0 # tehdaan alkuun 0 vuosi
            
            sdata0 = sdata0.sort_values(by=['AIKA']) #sorattaan jotta 0-vuosi alkuun
            sdata0 = sdata0.reset_index(drop=True) #resetoidaan indeksi
            sdata0['Age'][0] = sdata0['Age'][0]-(sdata0['AIKA'][1] - sdata0['AIKA'][0]) #tehdaan n-vuodelle oikea ika
    
            
            # if (sdata0['Age'][0]<0) | (sdata0['Age'][1]<0):
            #     sdata0['Age'][0]=80 # Jos ikÃ¤Ã¤ ei saada puustot-tiedostosta, arvotaan 80 v (tÃ¤llÃ¤ ei pitÃ¤isi olla juurikaan merkitystÃ¤ mihinkÃ¤Ã¤n??)
    
                # # print('Second row age = 0 !!! i = ' +str(i))
                # with open(outfile, "a") as myfile:
                #     myfile.write('\nSecond row age = 0 !!! i = ' +str(i) + 'kuvio: '+str(kuviolista[i]))
            
            
            sdata0['TotVol'][0] = tulotdata['volume0'][0] #luodulle aika=0 riville tilavuus tulotdatan vuosittaisita ilavuusista
            
            # PÃ¤ivitetÃ¤Ã¤n muut komponentit
            for c in range(0,len(comps)):
                if sdata0.loc[0,comps[c]]>0: #jos arvo on yli 0
                    sdata0.loc[0,comps[c]] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, comps[c]] ####

            #/users/asalmiva/hiilipolku_leena_silvi_mottitiedostojenteko.py            
            # PÃ¤ivitetÃ¤Ã¤n N
            if sdata0['tuotosKuolleet'][0] < sdata0['TotVol'][0]:
                sdata0.loc[0,'N'] = (1 + sdata0['tuotosKuolleet'][0]/sdata0['TotVol'][0])*sdata0['N'][0]       ##miten tama on laskettu 
                     
            for c in range(0,len(biom_comp)):
                if sdata0.loc[0,biom_comp[c]]>0: #jos arvo on yli 0
                    sdata0.loc[0,biom_comp[c]] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, biom_comp[c]] ####runkohukka pitäisi olla ianakin ennen päätehhakkuta eripäin tuo suhde kuin muilla jten pitäisi skaalata eri tavoin... kasvaa 35 vuotaan saakka ehkä ja sitten laksee ainkain lehtipuilla? tai kuuso+lehti
    
            repair = 1
    
        # # Jos harvennus tulee jo ennen Motti-tulosteen toista riviÃ¤, luodaan AIKA=1 -rivi,
        # # jotta pystytÃ¤Ã¤n tehdÃ¤ Motti-tiedosto, josta saadaan interpolointifunktiot
        # # parille ensimmÃ¤iselle vuodelle
    
        # if len(vuodet)>0:
        #     if vuodet[0]>1: # LisÃ¤ys tehdÃ¤Ã¤n vain, jos hakkuu tulee vuonna 2 tai myÃ¶hemmin
        #         if sdata0['AIKA'][1]>=vuodet[0]: # Jos aika=0 tai aika=1, lisÃ¤tÃ¤Ã¤n alkuun rivi, joka kopioitu ensimmÃ¤iseltÃ¤ olemassa olevalta riviltÃ¤, kÃ¤ytÃ¤nnÃ¶ssÃ¤ varhaisperkauksen jÃ¤lkeinen tilanne
        
        #             with open(outfile, "a") as myfile:
        #                 myfile.write('\nAika[1]>vuodet[0]!!! ---- HUOM!!!, i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
                    
        #             sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
                    
        #             sdata0['AIKA'][len(sdata0)-1] = 1
                    
        #             sdata0 = sdata0.sort_values(by=['AIKA'])
        #             sdata0 = sdata0.reset_index(drop=True)
        #             sdata0['Age'][1] = sdata0['Age'][0]+1
        
                    
        #             # if (sdata0['Age'][1]<0) | (sdata0['Age'][1]<0):
        #             #     sdata0['Age'][1]=81 # Jos ikÃ¤Ã¤ ei saada puustot-tiedostosta, arvotaan 80 v (tÃ¤llÃ¤ ei pitÃ¤isi olla juurikaan merkitystÃ¤ mihinkÃ¤Ã¤n??)
        
        #                 # # print('Second row age = 0 !!! i = ' +str(i))
        #                 # with open(outfile, "a") as myfile:
        #                 #     myfile.write('\nSecond row age = 0 !!! i = ' +str(i))
                    
                    
        #             sdata0['TotVol'][1] = tulotdata['volume1'][0] 
            
        #             repairB = 1
    
    
        # Tsekataan, onko TotVol=0 -rivejÃ¤
        
        a = 0
        zeroRow = next((m for m, x in enumerate(sdata0['TotVol'][a:]) if x==0), None)
        
        if zeroRow is not None:
            
            zeroRow = a + zeroRow
                    
            while (zeroRow is not None) & (zeroRow<len(sdata0)-1): # Ei tarkisteta viimeistÃ¤ riviÃ¤ tÃ¤ssÃ¤, koska tuottaa ongelmia. Tsekataan alempana.
                
                new_vol = tulotdata['volume'+str(int(sdata0['AIKA'][zeroRow]))][0]
                
                if new_vol > 5: # PÃ¤Ã¤tehakkuun tapauksessa 0-rivi on Tulot-tiedostossa vasta kohdassa aika+1
                    new_vol = tulotdata['volume'+str(int(sdata0['AIKA'][zeroRow+1]))][0]
                
                sdata0.loc[zeroRow,'TotVol'] = new_vol
                
                # Haetaan seuraava 0-rivi:
                    
                a = zeroRow + 1    
                zeroRow = next((m for m, x in enumerate(sdata0['TotVol'][a:]) if x==0), None)
                
                if zeroRow is not None:
                    zeroRow = a + zeroRow
                else:
                    break
    
    
    
        # Jos 0-biomassoja, tarkistetaan myÃ¶s 0-iÃ¤t ja korjataan ikÃ¤ vastaamaan AIKA-jaksoja
    
        a = 0
        while a < len(sdata0):
            zeroRow = next((m for m, x in enumerate(sdata0['Age'][a:]) if x==0), None)
            if zeroRow is not None:
                zeroRow = zeroRow+a
                if sdata0['Neulasetkg'][zeroRow]==0: # Muokataan vain, jos neulasissa 0-rivejÃ¤
                    nonZeroRow = next((m for m, x in enumerate(sdata0['Age'][zeroRow:]) if x!=0), None)
                    if nonZeroRow is not None: # Jos 0-ikÃ¤ viimeisellÃ¤ rivillÃ¤, ei tehdÃ¤ tarkistuksia/korjauksia
                        nonZeroRow = nonZeroRow + zeroRow
                        if sdata0['AIKA'][nonZeroRow]-sdata0['AIKA'][zeroRow] != sdata0['Age'][nonZeroRow]-sdata0['Age'][zeroRow]:
    
                            sdata0.loc[zeroRow, 'Age'] = int(sdata0['Age'][nonZeroRow] - (sdata0['AIKA'][nonZeroRow]-sdata0['AIKA'][zeroRow]))
                            with open(outfile, "a") as myfile:
                                myfile.write('\n0-age corrected to: ' + str(sdata0.loc[zeroRow, 'Age']) + ', i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
    
                    else:
                        with open(outfile, "a") as myfile:
                            myfile.write('\n0-age detected at last row!! i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
    
                a = zeroRow + 1
            else:
                a = len(sdata0)+1     
                
                
    
    
        # Tsekataan TotVol ja DomH nollat ja muokataan
        # TÃ¤hÃ¤n tarvii ehkÃ¤ suuremman alarajan?
        # PitÃ¤Ã¤ lisÃ¤tÃ¤ myÃ¶s tyhjien tsekkaus!
        
        for k in range(0,len(sdata0)):
            
            if (sdata0.loc[k,'TotVol']<0.5) | (np.isnan(sdata0.loc[k,'TotVol']) == True): #voi olla 0.01 tilavuuksiakin
                
                if k < len(sdata0)-1:
                    sdata0.loc[k,'TotVol'] = max(0.5, 0.1*sdata0.loc[k+1,'TotVol']) # LisÃ¤tÃ¤Ã¤n totvol nollien tilalle muu pieni luku (0.001 oli liian pieni, kasvu ei lÃ¤htenyt kÃ¤yntiin, joten kasvatettu 0.2:een) paitsi et tulee 0.5
                else:
                    sdata0.loc[k,'TotVol'] = 0.5
                    
            if (sdata0.loc[k,'DomH']<0.3) | (np.isnan(sdata0.loc[k,'DomH']) == True):
                 
                 if k < len(sdata0)-1:
                     sdata0.loc[k,'DomH'] = max(0.3, 0.3*sdata0.loc[k+1,'DomH']) # LisÃ¤tÃ¤Ã¤n domh nollien tilalle muu pieni luku
                 else:
                     sdata0.loc[k,'DomH'] = 0.3       
       
        
        # Adding values to "missing" biomass components with similar increase-% per year as TotVol
        # Neulaset viimeisenÃ¤, jotta sitÃ¤ voidaan kÃ¤yttÃ¤Ã¤ arvioimaan, tarviiko biomassoja tÃ¤ydentÃ¤Ã¤.
        # Poistettu tÃ¤stÃ¤ listasta 'tuotosKuolleet', koska se on kumulatiivinen..
        
        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg', \
                     'TukPL1', 'TukPL2', 'TukPL3', 'TukPL4',\
                         'KuiPL1', 'KuiPL2', 'KuiPL3', 'KuiPL4',\
                         'HukPL1', 'HukPL2', 'HukPL3', 'HukPL4',\
                             'tuotosTilavuus',\
                                 'Runkopuukg', 'Hukkapuukg',\
                                     'MeanDbh','medianH', 'BA', 'tuotos']  
        #comps = ['DomH','MeanDbh','medianH',\
        #          'BA','Kannot_juuretkg',\
        #                  'tuotosTilavuus',\
        #                      'TukPL1', 'TukPL2','TukPL3',\
        #                          'KuiPL1', 'KuiPL2','KuiPL3',\
        #                              'HukPL1', 'HukPL2','HukPL3',\
        #                                  'VolPL1', 'VolPL2', 'VolPL3']  

        #biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg'] #tama lisatty tahan elokuu 23 jotta namakin tsekataan

        # if repairB==1: # Jos lisÃ¤tty 1-rivi  kopioimalla, korjataan sen komponentit ensin
        #     for j in range(0,len(biom_comp)):
        #         comp = biom_comp[j]
        #         sdata0.loc[1, comp] = sdata0['TotVol'][1]/sdata0['TotVol'][2]*sdata0.loc[2, comp]    
            
        # if repair==1: # Jos lisÃ¤tty 0-rivi kopioimalla, korjataan sen komponentit ensin
        #     for j in range(0,len(biom_comp)):
        #         comp = biom_comp[j]
        #         sdata0.loc[0, comp] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, comp]
                
                            
    
                
        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg']
    
        for j in range(0,len(biom_comp)):
        
            comp = biom_comp[j]
            k=0
    
            while k < len(sdata0)-1: 
                
                zeroRow = next((m for m, x in enumerate(sdata0[comp][k:]) if x==0), None)
    
                if zeroRow is not None: # If zeroRow exists
                    
                    zeroRow = zeroRow + k
                    
                    # if sdata0['Neulasetkg'][zeroRow]==0: # Muiden komponenttien numeroita muokatan vain, jos neulasissa on nollaa. TÃ¤mÃ¤ poistettu, koska nÃ¤issÃ¤ komponenteissa voi olla nollaa, vaikka neulasmassa ei olisikaan nollaa..
                    nonZeroRow = next((m for m, x in enumerate(sdata0[comp][zeroRow:]) if x!=0), None)
                    
                    if nonZeroRow is not None:
                        
                        nonZeroRow = nonZeroRow + zeroRow
                        sdata0.loc[zeroRow, comp] = sdata0['TotVol'][zeroRow]/sdata0['TotVol'][nonZeroRow]*sdata0.loc[nonZeroRow, comp] # tää vaikuttaa siihen myös et 0 riville tulee 0.5 vol ja seuraava on jotain
                    
                        k=k+1
                        
                    else:
                        
                        k=k+1                          
    
                else:
                    k=k+1
                            
    
                        
                        # else:
                            
                            # Jotta ei toistu jokaisella komponentilla erikseen, 
                            # nÃ¤mÃ¤ korjataan vasta lopuksi: 
                            # kopioidaan loppuun rivit, jotka vastaavat pÃ¤Ã¤hakkuun jÃ¤lkeistÃ¤ kehitystÃ¤
                            
                            # EtsitÃ¤Ã¤n toisen pÃ¤Ã¤hakkuun jÃ¤lkeinen rivi
            
    
             
                
             
    
    
    
        # TehdÃ¤Ã¤n sarake 'treesp'=pÃ¤Ã¤puulaji volumen perusteella
        
        sdata0['treesp'] = pd.Series([], dtype='int64') #        sdata0['treesp'] = pd.Series([])

        
        for t in range(0,len(sdata0)):
            
            
            volPine = sdata0['VolPL1'][t]
            volSpruce = sdata0['VolPL2'][t]
            volBirch = sdata0['VolPL3'][t] + sdata0['VolPL4'][t]
        
            spe = 3
        
            if volPine >= max(volSpruce, volBirch):
                spe = 1
            
            if volSpruce > max(volPine, volBirch):
                spe = 2    
    
            sdata0.loc[t, 'treesp']=spe
    
                            
                            
        # Tsekataan, lÃ¶ytyykÃ¶ viimeisen pÃ¤Ã¤tehakkuun jÃ¤lkeen pelkkiÃ¤ 0-rivejÃ¤. Jos nÃ¤in,
        # kopioidaan aiemmat pÃ¤Ã¤tehakkuun jÃ¤lkeiset rivit loppuun.
        # Korjataan sen jÃ¤lkeen myÃ¶s AIKA-sarake
    
        if len(paateh)>0:
            paate1 =  next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>paateh[0]), None) 
            nonZeroRow = next((m for m, x in enumerate(sdata0['Neulasetkg'][paate1:]) if x!=0), None)
            
            if len(paateh)>1: # Jos useampi kuin yksi pÃ¤Ã¤tehakkuu
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>=paateh[1]), None) # tÃ¤ssÃ¤ aiemmin x>paateh[1] !!
                nonZeroRow = next((m for m, x in enumerate(sdata0['Neulasetkg'][paaterow:]) if x!=0), None)
            
                if nonZeroRow is None:
                    
                    paate2aika = sdata0['AIKA'][paaterow]
                    
                    b = 0
                    
                    for p in range(paaterow,len(sdata0)):
                        print(p, paate2aika, "paatehakkuun jalkeiset nollarivit")
    
                        currow = paate1 + p - paaterow -1 # LisÃ¤tty tÃ¤hÃ¤n -1, koska rivillÃ¤ 617 lisÃ¤tty x>=paateh[1]
    
                        if currow<b: # LisÃ¤ys 17.10: Joissakin tapauksissa yllÃ¤ toimii -1, toisissa 0, joten muokataan currow tarvittaessa...
                            currow=b
                            
                        sdata0.loc[p]=sdata0.loc[currow]
                        
                        if b==0:
                            sdata0.loc[p, 'AIKA'] = int(paate2aika)
                            sdata0.loc[p, 'tuotosKuolleet'] = sdata0.loc[p-1, 'tuotosKuolleet']
                            b=b+1
                        else:
                            sdata0.loc[p, 'AIKA'] = int(paate2aika + sdata0.loc[currow, 'AIKA'] - sdata0.loc[currow-b, 'AIKA'])
                            sdata0.loc[p, 'tuotosKuolleet'] = sdata0.loc[p-1, 'tuotosKuolleet']
                            b=b+1
                
            else:
                # Jos vain yksi pÃ¤Ã¤tehakkuu, ei voida kopioida rivejÃ¤ aiemman pÃ¤Ã¤tehakkuun jÃ¤lkeen..
                print('MitÃ¤s tÃ¤hÃ¤n? kuvio:', str(kuviolista[i]),paateh)   
                if nonZeroRow is None:
                    with open(outfile, "a") as myfile:
                        myfile.write('\nPÃ¤Ã¤tehakkuun jÃ¤lkeisiÃ¤ 0-rivejÃ¤ ei saada korjattua! i=' + str(i)+ ' kuvio:'+ str(kuviolista[i])+ ' paateh: '+str(paateh)) #Ã¤Ã¤ 
    
        
        # Tsekataan, lÃ¶ytyykÃ¶ viimeiseltÃ¤ riviltÃ¤ (aika=100), nollariviÃ¤.
        # Jos lÃ¶ytyy, poistetaan.
        
    
        s100 = sdata0[(sdata0['AIKA']==100) & (sdata0['Age']==0)] # Vai vol=0?
        
        if len(s100)>0:
            sdata0 = sdata0[sdata0['AIKA']!=100]
            
        # -------------------------------------------------------------------
        # LisÃ¤ys 13.9.2021
        # Tarkistetaan, olisiko ensimmÃ¤inen Motti-tiedosto 1-rivin mittainen:
        
        vuodet_p = vuodet[vuodet>0]
        
        if (len(vuodet_p)>0):
            
            sdata = sdata0[sdata0['AIKA'] < vuodet_p[0]]
            
            if len(sdata)==1:
                #sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
                nolla2_df=pd.DataFrame([sdata0.iloc[0]])
                sdata0 = pd.concat((sdata0,nolla2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = vuodet_p[0]-1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                sdata0.loc[1,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][1]))][0]
                sdata0.loc[1,'Age'] = sdata0['Age'][0]+(sdata0['AIKA'][1]-sdata0['AIKA'][0])
    
                for c in range(0,len(comps)):
                    if sdata0.loc[0,comps[c]]>0:
                        # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
                        # sitÃ¤ voisi supistaa..??
                        sdata0.loc[1,comps[c]] = sdata0['TotVol'][1]/sdata0['TotVol'][0]*sdata0.loc[0, comps[c]]
    
                for c in range(0,len(biom_comp)):
                    if sdata0.loc[0,biom_comp[c]]>0:
                        # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
                        # sitÃ¤ voisi supistaa..??
                        sdata0.loc[1,biom_comp[c]] = sdata0['TotVol'][1]/sdata0['TotVol'][0]*sdata0.loc[0, biom_comp[c]]
    
    
        # -------------------------------------------------------------------
        # LisÃ¤ys 9.3.2022
        # Tarkistetaan, olisiko *toinen* Motti-tiedosto 1-rivin mittainen:
        
        #if (len(vuodet_p)>0):
        #    
        #    if len(vuodet_p)==1:
        #        end_year=101
        #    else:
        #        end_year = vuodet_p[1]
        #        
        #    sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[0]) & (sdata0['AIKA'] < end_year)] # Huom. tÃ¤ssÃ¤ ei saa resetoida indeksiÃ¤, jotta jatko toimii
        #    
        #    if len(sdata)==1:
        #        sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # LisÃ¤tÃ¤Ã¤n rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
        #        sdata0['AIKA'][len(sdata0)-1] = end_year - 1
        # 
        #         sdata0 = sdata0.sort_values(by=['AIKA'])
        #         sdata0 = sdata0.reset_index(drop=True)
        #         
        #        sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
        #        
        #        sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
        #        
        #        if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa riviÃ¤
        #            sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
        #
        #            for c in range(0,len(comps)):
        #                if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
        #                    # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
        #                    # sitÃ¤ voisi supistaa..??
        #                    sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
        # 
        #            for c in range(0,len(biom_comp)):
        #                if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
        #                    # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
        #                    # sitÃ¤ voisi supistaa..??
        #                    sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
        
        
        # -------------------------------------------------------------------
        # Update 6.3.2023: checking all 1-row-motti-files except the last one
        
        for v in range(1, len(vuodet_p)):
            
            if (len(vuodet_p)>0):
                
                if len(vuodet_p)==1:
                    end_year=101
                else:
                    end_year = vuodet_p[v]                
                
                sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[v-1]) & (sdata0['AIKA'] < end_year)] # Huom. tÃ¤ssÃ¤ ei saa resetoida indeksiÃ¤ jotta jatko toimii
                
                if len(sdata)==1:
                    #sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # LisÃ¤tÃ¤Ã¤n rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
                    sdataind_df=pd.DataFrame([sdata0.iloc[sdata.index[0]]])
                    sdata0 = pd.concat((sdata0,sdataind_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                    sdata0['AIKA'][len(sdata0)-1] = end_year - 1
                    
                    sdata0 = sdata0.sort_values(by=['AIKA'])
                    sdata0 = sdata0.reset_index(drop=True)
                    
                    sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
                    
                    sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
                    
                    if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa rivia
                        sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
                        
                        for c in range(0,len(comps)):
                            if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
                                # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                                # sita voisi supistaa..??
                                sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
                        
                        for c in range(0,len(biom_comp)):
                            if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
                                # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                                # sita voisi supistaa..??
                                sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
        
        # Update 6.3.2023: Checking if the last motti-file would be a 1-row-file
        
        if (len(vuodet_p)>1):
            
            sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[len(vuodet_p)-1])] # Huom. tassa ei saa resetoida indeksia, jotta jatko toimii
            end_year = 101
            
            if (len(sdata)==1) & (sdata0['AIKA'][len(sdata0)-1]<100):
                
                #sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # Lisataan rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
                sdataind2_df=pd.DataFrame([sdata0.iloc[sdata.index[0]]])
                sdata0 = pd.concat((sdata0, sdataind2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                    
                sdata0['AIKA'][len(sdata0)-1] = end_year - 1
                
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
                
                sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
                
                if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa rivia
                    sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
                    
                    for c in range(0,len(comps)):
                        if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
                            # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                            # sita voisi supistaa..??
                            sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
                    
                    for c in range(0,len(biom_comp)):
                        if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
                            # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                            # sita voisi supistaa..??
                            sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
                            
                
        
        # -------------------------------------------------------------------
        # LisÃ¤ys 9.3.2022
        # Tarkistetaan, tuleeko perÃ¤kkÃ¤isten harvennusten takia liian paljon
        # vÃ¤heneviÃ¤ tilavuuksia 
        
        comps_notk = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg', \
                             'tuotosTilavuus','MeanDbh','medianH', 'BA', 'tuotos', 'TotVol']  
        
        
        if pera >= 0:
            
            vuodet_p = vuodet[vuodet>=0] # LisÃ¤ys 17.10.2022: vuodet>=0 (jos tÃ¤ssÃ¤ olisi vuodet>0, vuodet_p saattaisi jÃ¤Ã¤dÃ¤ tyhjÃ¤ksi -> virhe seuraavalla rivillÃ¤)
            
            sdata = sdata0[sdata0['AIKA'] < vuodet_p[0]]
            
            if (len(sdata)<=2) & (len(sdata)>0): # LisÃ¤ys 17.10.2022: len(sdata)>0, koska jos sdata on tyhjÃ¤ -> error
                
                if (sdata['TotVol'][0] - sdata['TotVol'][1]) > 0:
                    
                    sdata0.loc[0]=sdata0.loc[1]
                    sdata0.loc[0, 'AIKA'] = sdata0['AIKA'][1] -1
                    sdata0.loc[0, 'Age'] = sdata0['Age'][1] -1
                    
                    
                    for c in range(0,len(comps_notk)):
                        
                        comp = comps_notk[c]
                        sdata0.loc[0, comp] = sdata0.loc[0, comp] * 0.97 # Arvio, jotta 
                                            
    
    
        
        
        # Muutetaan Age=-1 --> Age=0, nÃ¤itÃ¤ lÃ¶ytyi suoraan Motti-datasta
        
        sdata0['Age'] = np.clip(sdata0['Age'], 0, np.nan) 
        
        # -------------------------------------------------------------------
        # Lisays 28.3.2023:
        # Tarkistetaan, loytyyko iasta perakkaisia nollia
        
        a = 0
        while a < len(sdata0):
            zeroRow = next((m for m, x in enumerate(sdata0['Age'][a:]) if x==0), None)
            
            if zeroRow is not None:
                zeroRow = zeroRow+a
                try:
                    if sdata0['Age'][zeroRow+1]==0:
                        sdata0.loc[zeroRow, 'Age'] = 1
                        sdata0.loc[zeroRow+1, 'Age'] = 1 + int(sdata0['AIKA'][zeroRow+1]-sdata0['AIKA'][zeroRow])
                        with open(outfile, "a") as myfile:
                            myfile.write('\nDouble zero ages!! i=' + str(i))
                        
                    a = zeroRow + 1
                except:
                    a = zeroRow + 1
                        
            else:
                a = len(sdata0)+1     
        
        
        
        
        # Varsinaisten Motti-tiedostojen teko alkaa tÃ¤stÃ¤
        
        if len(sdata0)>0:
            
            if len(vuodet)>0:
                
                for j in range(0,n_mottifiles):
                
                    if j == 0:
                        sdata = sdata0[sdata0['AIKA'] < vuodet[j]]
                        if len(sdata)==0:
                            continue
                    
                    if (j>0) & (j<n_mottifiles-1):
                        sdata = sdata0[(sdata0['AIKA'] >= vuodet[j-1]) & (sdata0['AIKA'] < vuodet[j])]
                        sdata = sdata.reset_index(drop=True)
                        if len(sdata)==0:
                            continue
                    
                    if j == n_mottifiles-1:
                        sdata = sdata0[sdata0['AIKA'] >= vuodet[j-1]]
                        sdata = sdata.reset_index(drop=True)
                        if len(sdata)==0:
                            continue
                    
                    kasvatus = pd.Series(np.ones(len(sdata)))
                    vuosi = sdata['AIKA']
                    ika = sdata['Age']
                    N = sdata['N']
                    PPA = sdata['BA']
                    Hg = sdata['medianH']
                    Dg = sdata['MeanDbh']
                    Hdom = sdata['DomH']
                    Tilavuus = sdata['TotVol']
                    Tukki = sdata['TukPL1']+sdata['TukPL2']+sdata['TukPL3']+sdata['TukPL4']
                    Kuitu = sdata['KuiPL1']+sdata['KuiPL2']+sdata['KuiPL3']+sdata['KuiPL4']
                    Hukka = sdata['HukPL1']+sdata['HukPL2']+sdata['HukPL3']+sdata['HukPL4']
                    Tuotos = sdata['tuotosKuolleet'] + sdata['tuotosTilavuus']
                    Kuolleisuus = Tuotos*0. # Ei kÃ¤ytetÃ¤ Susissa
                    runko_aines = sdata['Runkopuukg']/1000. # YksikkÃ¶ tn
                    runko_hukka = sdata['Hukkapuukg']/1000.
                    elavat_oksat = sdata['Oksatkg']*0.8/1000. # Oksatkg = elÃ¤vÃ¤t ja kuolleet oksat yhteensÃ¤
                    kuolleet_oksat = sdata['Oksatkg']*0.2/1000.
                    lehdet = sdata['Neulasetkg']/1000.
                    kannot = sdata['kannot']/1000.
                    juuret = sdata['>2mm_juuretkg']/1000.
                    hienojuuret = sdata['h_juuretkg']/1000.
                    
                    # Tsekkaa, ettÃ¤ tÃ¤mÃ¤ on ok!!!
                    for t in range(0,len(sdata)):
                        if t==0:
                            Tuotos[0]=sdata['tuotosTilavuus'][0]
                        
                        else:
                            Tuotos[t] = sdata['tuotosTilavuus'][t] + (sdata['tuotosKuolleet'][t] - sdata['tuotosKuolleet'][t-1])
                    
                    # ---tÃ¤hÃ¤n asti
                    
                    data = pd.DataFrame(np.transpose(np.array([kasvatus, vuosi, ika, N, PPA, Hg, Dg, Hdom, Tilavuus, Tukki, Kuitu, \
                                         Hukka, Tuotos, Kuolleisuus, runko_aines, runko_hukka, \
                                         elavat_oksat, kuolleet_oksat, lehdet, kannot, juuret, hienojuuret])), \
                                         columns=['Kasvatus', 'Vuosi', 'Ikä', 'N', 'PPA', 'Hg', 'Dg', 'Hdom', 'Tilavuus', 'Tukki', 'Kuitu',\
                                                  'Hukka', 'Tuotos', 'Kuolleisuus', 'runko(aines)', 'runko(hukka)', \
                                                  'elävät oksat', 'kuolleet oksat', 'lehdet', 'Kannot', 'Juuret >2mm', 'Hienojuuret'])
                    
                    standInfo = kuviot[kuviot['KUVIO']==kuviolista[i]]
                    standInfo = standInfo.reset_index(drop=True)
                    
                    """ TSEKKAA, ETTÃ„ PUULAJI TULEE OIKEIN! """
                    
                    
                    
                    if (project == 'life') | (project == 'hiilipolku'):
                        spe = sdata['treesp'][len(sdata)-1]
                    
                    else: 
                        if (project =='vmi') | (project=='suo'):
                            
                            spe = int(sdata['treesp'][len(sdata)-1]) # viimeisen rivin puulaji
                            # treesp=kdata['Selite'][0][28:33]
                            # spe = 3
                            
                            # if treesp=='MÃ¤nty':
                            #     spe = 1
                            # if treesp=='Kuusi':
                            #     spe = 2       
                                
                            
                        else:
                            
                            treesp = pdata2['Selite'][0][28:33]
                            
                            spe = 3
                            
                            if treesp=='Mänty': #onk tÃ¤ a va Ã¤#
                                spe = 1
                            if treesp=='Kuusi':
                                spe = 2 
                    
                    
                    # print(kuviolista[i], ': puulaji', spe)
                    with open(outfile, "a") as myfile:
                        myfile.write('\n' + str(kuviolista[i]) +  ': puulaji ' + str(spe)+ ' n'+str(j))
                    
                    
                    fout = mottifolder + '/' + str(kuviolista[i]) + '_n' + str(j) + '.xls'
                    
                    cols = pd.DataFrame([['','','','',spe,'','','','','','','','','']],columns=['Kasvatus',	'Vuosi',	'id Harvennus',	'Harvennus',	'id Puulaji',	'Puulaji',	'Tukki[m³/ha]',	'Pikkutukki[m³/ha]',	'Kuitu[m³/ha]',	'Energiapuu, runko(aines)[m³/ha]',	'Energiapuu, runko(hukka)[m³/ha]',	'Energiapuu, oksat(e)[m³/ha]',	'Energiapuu, oksat(k)[m³/ha]',	'Energiapuu, kannot ja juuret[m³/ha]'])
                    
                    writer = pd.ExcelWriter(fout, engine = 'xlsxwriter')
                    data.to_excel(writer, index=False, sheet_name = 'Puustotunnukset')
                    cols.to_excel(writer, index=False, sheet_name = r'Kertymät')
                    #writer.save()
                    writer.close()
            else:
                # print(str(kuviolista[i]) + ': poistumat puuttuvat! i=' + str(i))
                with open(outfile, "a") as myfile:
                    myfile.write('\n' + str(kuviolista[i]) + ': poistumat puuttuvat! i=' + str(i))
                
                sdata = sdata0.copy()
                
                kasvatus = pd.Series(np.ones(len(sdata)))
                vuosi = sdata['AIKA']
                ika = sdata['Age']
                N = sdata['N']
                PPA = sdata['BA']
                Hg = sdata['medianH']
                Dg = sdata['MeanDbh']
                Hdom = sdata['DomH']
                Tilavuus = sdata['TotVol']
                Tukki = sdata['TukPL1']+sdata['TukPL2']+sdata['TukPL3']+sdata['TukPL4']
                Kuitu = sdata['KuiPL1']+sdata['KuiPL2']+sdata['KuiPL3']+sdata['KuiPL4']
                Hukka = sdata['HukPL1']+sdata['HukPL2']+sdata['HukPL3']+sdata['HukPL4']
                Tuotos = sdata['tuotos']
                Kuolleisuus = Tuotos*0. # Ei kÃ¤ytetÃ¤ Susissa
                runko_aines = sdata['Runkopuukg']/1000. # YksikkÃ¶ tn
                runko_hukka = sdata['Hukkapuukg']/1000.
                elavat_oksat = sdata['Oksatkg']*0.8/1000. # Oksatkg = elÃ¤vÃ¤t ja kuolleet oksat yhteensÃ¤
                kuolleet_oksat = sdata['Oksatkg']*0.2/1000.
                lehdet = sdata['Neulasetkg']/1000.
                kannot = sdata['kannot']/1000.
                juuret = sdata['>2mm_juuretkg']/1000.
                hienojuuret = sdata['h_juuretkg']/1000.
                
                data = pd.DataFrame(np.transpose(np.array([kasvatus, vuosi, ika, N, PPA, Hg, Dg, Hdom, Tilavuus, Tukki, Kuitu, \
                                     Hukka, Tuotos, Kuolleisuus, runko_aines, runko_hukka, \
                                     elavat_oksat, kuolleet_oksat, lehdet, kannot, juuret, hienojuuret])), \
                                     columns=['Kasvatus', 'Vuosi', 'Ikä', 'N', 'PPA', 'Hg', 'Dg', 'Hdom', 'Tilavuus', 'Tukki', 'Kuitu',\
                                              'Hukka', 'Tuotos', 'Kuolleisuus', 'runko(aines)', 'runko(hukka)', \
                                              'elävät oksat', 'kuolleet oksat', 'lehdet', 'Kannot', 'Juuret >2mm', 'Hienojuuret'])
                
                standInfo = kuviot[kuviot['KUVIO']==kuviolista[i]]
                standInfo = standInfo.reset_index(drop=True)
                
                
                """ TSEKKAA, ETTÃ„ PUULAJI TULEE OIKEIN! """
                if (project == 'life') | (project == 'hiilipolku'):
                    spe = sdata['treesp'][len(sdata)-1]
                
                else:
                    
                    if project =='vmi':
                        
                        # spe = int(sdata['treesp'][len(sdata)-1]) # viimeisen rivin puulaji
                       
                        treesp=kdata['Selite'][0][28:33]
                        spe = 3
                        
                        if treesp=='Mänty': #Ã¤
                            spe = 1
                        if treesp=='Kuusi':
                            spe = 2    
                    else:
                        treesp = pdata2['Selite'][0][28:33]
                        
                        spe = 3
                        
                        if treesp=='Mänty':
                            spe = 1
                        if treesp=='Kuusi':
                            spe = 2 
                    
                # print(kuviolista[i], ': puulaji ', spe)
                with open(outfile, "a") as myfile:
                    myfile.write('\n' + str(kuviolista[i]) +  ':' + str(spe)+ ' n0')
                
                
                
                fout = mottifolder + '/' + str(kuviolista[i]) + '_n0.xls'
                
                cols = pd.DataFrame([['','','','',spe,'','','','','','','','','']],columns=['Kasvatus',	'Vuosi',	'id Harvennus',	'Harvennus',	'id Puulaji',	'Puulaji',	'Tukki[m³/ha]',	'Pikkutukki[m³/ha]',	'Kuitu[m³/ha]',	'Energiapuu, runko(aines)[m³/ha]',	'Energiapuu, runko(hukka)[m³/ha]',	'Energiapuu, oksat(e)[m³/ha]',	'Energiapuu, oksat(k)[m³/ha]',	'Energiapuu, kannot ja juuret[m³/ha]'])
                
                writer = pd.ExcelWriter(fout, engine = 'xlsxwriter')
                data.to_excel(writer, index=False, sheet_name = 'Puustotunnukset')
                cols.to_excel(writer, index=False, sheet_name = r'Kertymät')
                #writer.save()
                writer.close()

def create_motti_files_silvi_lp4_vol01(project, scen, maku=None):       #alkutilavuus, annetaan olla 0.1tai 0.1*Totvol[1]lp4 skaalaukseen otettu lp4 mukaan  # Huom! ohje>1
    '''
    Modified: January 2023
        
    Create Motti input files for Susi simulations from "Skene-Motti" results.
    
    Works technically, but contains many small issues that should be fixed.
    
    project = 'life' (Hydrology Life -project) or 'suo' (SUO-project)
    project-parameter affects only main tree species (and column names)!
    
    
    Modified 8.3.2022:
        - korjattu pÃ¤Ã¤tehakkuuseen liittyviÃ¤ ongelmia; hdom pitÃ¤isi nyt saada
        jÃ¤rkeviÃ¤ arvoja, biomassojen ei pitÃ¤isi jÃ¤Ã¤dÃ¤ pÃ¤Ã¤tehakkuun jÃ¤lkeen nollaksi
        (ainakaan toisen pÃ¤Ã¤tehakkuun jÃ¤lkeen)
        
    Issues to fix:
        - Jos pÃ¤Ã¤tehakkuu tulee vuonna 2, pitÃ¤Ã¤ vuosille 0 ja 1 hakea muut komponentit
        totvol perusteella: hdom saa tÃ¤ssÃ¤ liian suuria arvoja. Tuotos jÃ¤Ã¤ nollaksi?!
        - See comments below...
        
    
    '''
    
    inputmainfol=r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan_vol01/'
    mottifolder = inputmainfol+scen[:-6]+'/'+scen+'/'
    if not os.path.exists(inputmainfol):
        os.mkdir(inputmainfol)
    if not os.path.exists(inputmainfol+scen[:-6]+'/'):
        os.mkdir(inputmainfol+scen[:-6]+'/')
    if not os.path.exists(mottifolder):
        os.mkdir(mottifolder)
            
    outfile = mottifolder + '/out.txt'
    
    project = 'hiilipolku'
    
    kuviot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', encoding='latin1', sep=';')
    puustot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_puustot.csv', encoding='latin1', sep=';')
    poistumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_poistumat.csv', encoding='latin1', sep=';')
    tapahtumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_tapahtumat.csv', encoding='latin1', sep=';', index_col=False)
    tulot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kasvut.csv', encoding='latin1', sep=';')
    
    
    if (project == 'life') | (project =='vmi') | (project =='hiilipolku'):
        kuviot = kuviot.rename(columns={'kuvio':'KUVIO'})
        puustot = puustot.rename(columns={'kuvio':'KUVIO',
                                          '_2mm_juuretkg':'>2mm_juuretkg',})
        poistumat = poistumat.rename(columns={'kuvio':'KUVIO'})
        tapahtumat = tapahtumat.rename(columns={'kuvio':'KUVIO',
                                      '_2003': '2003',
                                      '_2004': '2004',
                                      '_2005': '2005'})
        tulot = tulot.rename(columns={'kuvio':'KUVIO'})
    
    
    kuviolista = list(set(kuviot['KUVIO'][kuviot['kuivatustilanne']>=3])) #tassa valitaan 3 Ojitettu kangas ja muut?DRAINAGESTATE	, ojitettu kangas jos pÃ¤Ã¤ryhmÃ¤1 ja alaryhmÃ¤ 3 ja kuivatustilanne 5?? 1 Ojittamaton kangas, 2 Soistunut kangas,3 Ojitettu kangas,	6 luonnontilainen suo,	7 Ojikko, 8	Muuttuma,9	Turvekangas
    
    
    for i in range(0,len(kuviolista)):
    
        print(i, 'stand: ', kuviolista[i])
        
        # alta poistettu ohjeeseen ja skenaarioon viittaavat rajaukset,
        # oletetaan, ettÃ¤ hiilipolussa vain yksi ohje per kuvio
        
        sdata0 = puustot[(puustot['KUVIO']==kuviolista[i])]
        pdata = poistumat[(poistumat['KUVIO']==kuviolista[i]) & (poistumat['HARV']>1)] #mita harv koodit meinaa?
        pdata2 = poistumat[(poistumat['KUVIO']==kuviolista[i])]
        tdata = tapahtumat[(tapahtumat['KUVIO']==kuviolista[i])]
        tulotdata = tulot[(tulot['KUVIO']==kuviolista[i])]
        kdata = kuviot[kuviot['KUVIO']==kuviolista[i]]
        
        sdata0 = sdata0.reset_index(drop=True)
        pdata = pdata.reset_index(drop=True)
        pdata2 = pdata2.reset_index(drop=True)
        tdata = tdata.reset_index(drop=True)
        tulotdata = tulotdata.reset_index(drop=True)
        kdata = kdata.reset_index(drop=True)
        
        tdata = tdata[(tdata['2003']>-1) | (tdata['2004']>-1) | (tdata['2005']>-1)]
        
        ensih = tdata[tdata['2003']>-1]['2003'].values
        harv = tdata[tdata['2004']>-1]['2004'].values
        paateh = tdata[tdata['2005']>-1]['2005'].values
        #voisi kayda tsekkaamassa onko taimikonhoito tai varhaisperkaus niin sitten ei muuttaisi puulajia tai ottaisi sen muualta et sailyy kuusena tai mantyna? 
        #tdic={'_2001':'lannoitus','_2002':'kunnostusojitus','_2003':'ensiharvennus','_2004':'harvennus','_2005':'paatehakkuu','_2006':'varhaiperkaus','_2007':'taimikonhoito','_2008':'viljely','_2009':'luontainen uudistuminen'}
        # taimh = tdata[tdata['_2007']>-1]['_2007'].values
        # Pudotetaan viimeisen vuoden tapahtumat pois, koska sekoittavat
        # myÃ¶hemmin tulevan koodin
        ensih = ensih[ensih<100]
        harv = harv[harv<100]
        paateh = paateh[paateh<100]
        
        vuodet = np.unique(np.sort(np.concatenate((ensih, harv, paateh))))
        
        ap = 0
        
        pera = -1
        
        # Tsekataan, onko tapahtumissa perÃ¤kkÃ¤isiÃ¤ vuosia ja poistetaan niistÃ¤ jÃ¤lkimmÃ¤inen
        for v in range(0,len(vuodet)-1):
            if vuodet[v+1]-vuodet[v]==1:
                
                ensih = np.delete(ensih, np.where(ensih == vuodet[v+1]))
                harv = np.delete(harv, np.where(harv == vuodet[v+1]))
                paateh = np.delete(paateh, np.where(paateh == vuodet[v+1]))
                
                vuodet = np.delete(vuodet, [v+1], None)
                
                pera = v
                
                ap = ap + 1
                
            if v+1>=len(vuodet)-1:
                break
        
        if ap>1:
            with open(outfile, "a") as myfile:
                myfile.write('\n\nKaksi kertaa perakkaiset vuodet! stand=' + str(kuviolista[i]) + ', i=' + str(i))
    
        n_mottifiles = len(vuodet) + 1
    
    
        if len(sdata0)==0: # skipataan kuvio, jos dataa ei ole
            continue
    
        if len(sdata0)>21:
            print('sdata0: liikaa rivejÃ¤?!')
        
    
                                    
        # Luodaan pÃ¤Ã¤tehakkuun jÃ¤lkeiset rivit ja korjataan rivin kopioimisesta aiheutuneet
        # luvut
        
        # muokattu 6.8.2021: poistettu tuotos, tuotosKuolleet, tuotosTilavuus comps:sta
        comps = ['DomH','MeanDbh','medianH',\
                          'BA','Kannot_juuretkg',\
                              'tuotosTilavuus',\
                              'TukPL1', 'TukPL2','TukPL3','TukPL4',\
                                  'KuiPL1', 'KuiPL2','KuiPL3','KuiPL4',\
                                      'HukPL1', 'HukPL2','HukPL3','HukPL4',\
                                          'VolPL1', 'VolPL2', 'VolPL3','VolPL4',]  #lisatty 091023 PL4 ositteet mukaan

        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg'] #tama lisatty tahan elokuu 23 jotta namakin tsekataan

        for h in range(0,len(paateh)):
            
            sdata_paateh= sdata0[sdata0['AIKA']==paateh[h]].reset_index(drop=True)
            
            go = False
            
            if len(sdata_paateh)>0:
                if sdata_paateh.loc[0, 'Age']>5:
                    go = True
        
            if (len(sdata_paateh)==0) | (go==True): # Luodaan rivit vain, jos ensimmÃ¤istÃ¤ riviÃ¤ ei ole mukana
                
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>paateh[h]), None)
                
                #sdata0 = sdata0.append([sdata0.iloc[paaterow]], ignore_index=True) # Rivi menee taulukon loppuun
                paaterow_df=pd.DataFrame([sdata0.iloc[paaterow]])
                sdata0 = pd.concat((sdata0,paaterow_df),axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = paateh[h]+1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==paateh[h]+1), None)
                
                # sdata0.loc[paaterow,'Age'] = 1 # TÃ¤mÃ¤ korvattu alla olevalla
                
                if paaterow>0:
                    if sdata0.loc[paaterow-1, 'Age']<5:
                        sdata0.loc[paaterow,'Age'] = sdata0.loc[paaterow-1, 'Age']+1
                    else:
                        sdata0.loc[paaterow,'Age'] = sdata0.loc[paaterow+1, 'Age']-1
                
                sdata0.loc[paaterow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][paaterow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[paaterow,comps[c]]>0:
                        sdata0.loc[paaterow,comps[c]] = sdata0['TotVol'][paaterow]/sdata0['TotVol'][paaterow+1]*sdata0.loc[paaterow+1, comps[c]]
                        print(kuviolista[i],comps[c], "rivi 303 ongelma")

            
                if sdata0['tuotosKuolleet'][paaterow] < sdata0['TotVol'][paaterow]:
                    sdata0.loc[paaterow,'N'] = (1 + sdata0['tuotosKuolleet'][paaterow]/sdata0['TotVol'][paaterow])*sdata0['N'][paaterow]        
                    
                    if paaterow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][paaterow] > sdata0['N'][paaterow-1]:
                            if paaterow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[paaterow,'N'] = (sdata0['N'][paaterow-1] + sdata0['N'][paaterow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[paaterow,'N'] = 0.99*sdata0['N'][paaterow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.
                            
                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[paaterow,biom_comp[c]]>0: 
                        sdata0.loc[paaterow,biom_comp[c]] = sdata0['TotVol'][paaterow]/sdata0['TotVol'][paaterow+1]*sdata0.loc[paaterow+1, biom_comp[c]] ####
                
            
        # Luodaan ensiharvennuksen jÃ¤lkeiset rivit (jos vuosi ei osu kohdallaan muuten)
        # EsimerkissÃ¤ nÃ¤ytti siltÃ¤, ettÃ¤ ensiharvennusvuoden riville tulostuu harventamaton tilavuus! Siksi tÃ¤ssÃ¤ kÃ¤sitelty riviÃ¤+1
        # Toisessa esimerkissÃ¤ tulostui harvennettu tilavuus! Joten vaihtelee.
        
        for h in range(0,len(ensih)):
    
    
            ensihrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==ensih[h]+1), None)
            ensihrow2 = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>ensih[h]), None)
            
            if ensihrow is None:
                
                #sdata0 = sdata0.append([sdata0.iloc[ensihrow2]], ignore_index=True) # Rivi menee taulukon loppuun
                ensihrow2_df=pd.DataFrame([sdata0.iloc[ensihrow2]])
                sdata0 = pd.concat((sdata0,ensihrow2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = ensih[h]+1
                
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                ensihrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==ensih[h]+1), None)
                
                # Jos ikÃ¤ lÃ¶ytyy riviÃ¤ alaempaa:
                if sdata0.loc[ensihrow+1,'Age']>0:
                    sdata0.loc[ensihrow,'Age'] = sdata0.loc[ensihrow+1,'Age']-(sdata0['AIKA'][ensihrow+1]-sdata0['AIKA'][ensihrow])
                # Muussa tapauksessa ikÃ¤ lÃ¶ytyy riviÃ¤ ylempÃ¤Ã¤:
                else:
                    if ensihrow>0:
                        sdata0.loc[ensihrow,'Age'] = sdata0.loc[ensihrow-1,'Age']+(sdata0['AIKA'][ensihrow]-sdata0['AIKA'][ensihrow-1])
                    else:
                        # Jos ei lÃ¶ydy kummaltakaan puolelta:
                        with open(outfile, "a") as myfile:
                            myfile.write('\nEnsiharvennuksen jÃ¤lkeistÃ¤ ikÃ¤Ã¤ ei lÃ¶ydy! kuvio=' + str(kuviolista[i]) + ' i:'+ str(i))
                        
                    
                sdata0.loc[ensihrow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][ensihrow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[ensihrow,comps[c]]>0:
                        sdata0.loc[ensihrow,comps[c]] = sdata0['TotVol'][ensihrow]/sdata0['TotVol'][ensihrow+1]*sdata0.loc[ensihrow+1, comps[c]]
                
                # PÃ¤ivitetÃ¤Ã¤n N
                if sdata0['tuotosKuolleet'][ensihrow] < sdata0['TotVol'][ensihrow]:
                    sdata0.loc[ensihrow,'N'] = (1 + sdata0['tuotosKuolleet'][ensihrow]/sdata0['TotVol'][ensihrow])*sdata0['N'][ensihrow]        
    
                    if ensihrow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][ensihrow] > sdata0['N'][ensihrow-1]:
                            if ensihrow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[ensihrow,'N'] = (sdata0['N'][ensihrow-1] + sdata0['N'][ensihrow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[ensihrow,'N'] = 0.99*sdata0['N'][ensihrow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.

                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[ensihrow,biom_comp[c]]>0: 
                        sdata0.loc[ensihrow,biom_comp[c]] = sdata0['TotVol'][ensihrow]/sdata0['TotVol'][ensihrow+1]*sdata0.loc[ensihrow+1, biom_comp[c]] ####
                            
        
        # Luodaan harvennuksen jÃ¤lkeiset rivit (jos vuosi ei osu kohdalleen muuten)
        # Harvennusvuoden rivillÃ¤ on harventamaton tilavuus. Siksi kÃ¤sitelty riviÃ¤+1.
        
        for h in range(0,len(harv)):
    
    
            hrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==harv[h]+1), None)
            hrow2 = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>harv[h]), None)
            
            if hrow is None:
                
                #sdata0 = sdata0.append([sdata0.iloc[hrow2]], ignore_index=True) # Rivi menee taulukon loppuun
                hrow2_df=pd.DataFrame([sdata0.iloc[hrow2]])
                sdata0 = pd.concat((sdata0,hrow2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = harv[h]+1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                hrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==harv[h]+1), None)
                
                # Jos ikÃ¤ lÃ¶ytyy riviÃ¤ alaempaa:
                if sdata0.loc[hrow+1,'Age']>0:
                    sdata0.loc[hrow,'Age'] = sdata0.loc[hrow+1,'Age']-(sdata0['AIKA'][hrow+1]-sdata0['AIKA'][hrow])
                # Muussa tapauksessa ikÃ¤ lÃ¶ytyy riviÃ¤ ylempÃ¤Ã¤:
                else:
                    if hrow>0:
                        sdata0.loc[hrow,'Age'] = sdata0.loc[hrow-1,'Age']+(sdata0['AIKA'][hrow]-sdata0['AIKA'][hrow-1])
                    else:
                    # Jos ei lÃ¶ydy kummaltakaan puolelta:
                        with open(outfile, "a") as myfile:
                            myfile.write('\nHarvennuksen jÃ¤lkeistÃ¤ ikÃ¤Ã¤ ei lÃ¶ydy! i=' + str(i) + ' kuvio:'+ str(kuviolista[i]))
                
                    
                sdata0.loc[hrow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][hrow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[hrow,comps[c]]>0:
                        sdata0.loc[hrow,comps[c]] = sdata0['TotVol'][hrow]/sdata0['TotVol'][hrow+1]*sdata0.loc[hrow+1, comps[c]]

                        
                # Korjataan N:
                if sdata0['tuotosKuolleet'][hrow] < sdata0['TotVol'][hrow]:
                    sdata0.loc[hrow,'N'] = (1 + sdata0['tuotosKuolleet'][hrow]/sdata0['TotVol'][hrow])*sdata0['N'][hrow]
                    
                    if hrow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][hrow] > sdata0['N'][hrow-1]:
                            if hrow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[hrow,'N'] = (sdata0['N'][hrow-1] + sdata0['N'][hrow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[hrow,'N'] = 0.99*sdata0['N'][hrow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.
                            
                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[hrow,biom_comp[c]]>0: 
                        sdata0.loc[hrow,biom_comp[c]] = sdata0['TotVol'][hrow]/sdata0['TotVol'][hrow+1]*sdata0.loc[hrow+1, biom_comp[c]] ####
    
        # Tsekataan, ettei ole ekan rivin TotVol pyÃ¶ristyksen kanssa ongelmaa:
        
        if (sdata0['TotVol'][0]>10.) & (-0.1 < sdata0['TotVol'][1]-sdata0['TotVol'][0]<0):
            sdata0.loc[0,'TotVol'] = sdata0['TotVol'][1]
                    
          
        # Jos yllÃ¤ olevien korjausten jÃ¤lkeekin aika on alussa > 1,
        # luodaan aika=0 -rivi.
    
        
        repair = 0
        repairB = 0
        
        if (sdata0['AIKA'][0]>1) & (sdata0['Age'][0]>=sdata0['AIKA'][0]): # Jos aika=0 tai aika=1, lisÃ¤tÃ¤Ã¤n alkuun rivi, joka kopioitu ensimmÃ¤iseltÃ¤ olemassa olevalta riviltÃ¤, kÃ¤ytÃ¤nnÃ¶ssÃ¤ usein varhaisperkauksen jÃ¤lkeinen tilanne
            # print('Aika alussa > 1, i=' + str(i))
            with open(outfile, "a") as myfile:
                myfile.write('\nAika alussa > 1, i=' + str(i) + ' kuvio: '+ str(kuviolista[i]))
             
            #sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
            nolla_df=pd.DataFrame([sdata0.iloc[0]])
            sdata0 = pd.concat((sdata0,nolla_df),axis=0, ignore_index=True) # Rivi menee taulukon loppuun
            sdata0['AIKA'][len(sdata0)-1] = 0 # tehdaan alkuun 0 vuosi
            
            sdata0 = sdata0.sort_values(by=['AIKA']) #sorattaan jotta 0-vuosi alkuun
            sdata0 = sdata0.reset_index(drop=True) #resetoidaan indeksi
            sdata0['Age'][0] = sdata0['Age'][0]-(sdata0['AIKA'][1] - sdata0['AIKA'][0]) #tehdaan n-vuodelle oikea ika
    
            
            # if (sdata0['Age'][0]<0) | (sdata0['Age'][1]<0):
            #     sdata0['Age'][0]=80 # Jos ikÃ¤Ã¤ ei saada puustot-tiedostosta, arvotaan 80 v (tÃ¤llÃ¤ ei pitÃ¤isi olla juurikaan merkitystÃ¤ mihinkÃ¤Ã¤n??)
    
                # # print('Second row age = 0 !!! i = ' +str(i))
                # with open(outfile, "a") as myfile:
                #     myfile.write('\nSecond row age = 0 !!! i = ' +str(i) + 'kuvio: '+str(kuviolista[i]))
            
            
            sdata0['TotVol'][0] = tulotdata['volume0'][0] #luodulle aika=0 riville tilavuus tulotdatan vuosittaisita ilavuusista
            
            # PÃ¤ivitetÃ¤Ã¤n muut komponentit
            for c in range(0,len(comps)):
                if sdata0.loc[0,comps[c]]>0: #jos arvo on yli 0
                    sdata0.loc[0,comps[c]] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, comps[c]] ####

            #/users/asalmiva/hiilipolku_leena_silvi_mottitiedostojenteko.py            
            # PÃ¤ivitetÃ¤Ã¤n N
            if sdata0['tuotosKuolleet'][0] < sdata0['TotVol'][0]:
                sdata0.loc[0,'N'] = (1 + sdata0['tuotosKuolleet'][0]/sdata0['TotVol'][0])*sdata0['N'][0]       ##miten tama on laskettu 
                     
            for c in range(0,len(biom_comp)):
                if sdata0.loc[0,biom_comp[c]]>0: #jos arvo on yli 0
                    sdata0.loc[0,biom_comp[c]] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, biom_comp[c]] ####runkohukka pitäisi olla ianakin ennen päätehhakkuta eripäin tuo suhde kuin muilla jten pitäisi skaalata eri tavoin... kasvaa 35 vuotaan saakka ehkä ja sitten laksee ainkain lehtipuilla? tai kuuso+lehti
    
            repair = 1
    
        # # Jos harvennus tulee jo ennen Motti-tulosteen toista riviÃ¤, luodaan AIKA=1 -rivi,
        # # jotta pystytÃ¤Ã¤n tehdÃ¤ Motti-tiedosto, josta saadaan interpolointifunktiot
        # # parille ensimmÃ¤iselle vuodelle
    
        # if len(vuodet)>0:
        #     if vuodet[0]>1: # LisÃ¤ys tehdÃ¤Ã¤n vain, jos hakkuu tulee vuonna 2 tai myÃ¶hemmin
        #         if sdata0['AIKA'][1]>=vuodet[0]: # Jos aika=0 tai aika=1, lisÃ¤tÃ¤Ã¤n alkuun rivi, joka kopioitu ensimmÃ¤iseltÃ¤ olemassa olevalta riviltÃ¤, kÃ¤ytÃ¤nnÃ¶ssÃ¤ varhaisperkauksen jÃ¤lkeinen tilanne
        
        #             with open(outfile, "a") as myfile:
        #                 myfile.write('\nAika[1]>vuodet[0]!!! ---- HUOM!!!, i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
                    
        #             sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
                    
        #             sdata0['AIKA'][len(sdata0)-1] = 1
                    
        #             sdata0 = sdata0.sort_values(by=['AIKA'])
        #             sdata0 = sdata0.reset_index(drop=True)
        #             sdata0['Age'][1] = sdata0['Age'][0]+1
        
                    
        #             # if (sdata0['Age'][1]<0) | (sdata0['Age'][1]<0):
        #             #     sdata0['Age'][1]=81 # Jos ikÃ¤Ã¤ ei saada puustot-tiedostosta, arvotaan 80 v (tÃ¤llÃ¤ ei pitÃ¤isi olla juurikaan merkitystÃ¤ mihinkÃ¤Ã¤n??)
        
        #                 # # print('Second row age = 0 !!! i = ' +str(i))
        #                 # with open(outfile, "a") as myfile:
        #                 #     myfile.write('\nSecond row age = 0 !!! i = ' +str(i))
                    
                    
        #             sdata0['TotVol'][1] = tulotdata['volume1'][0] 
            
        #             repairB = 1
    
    
        # Tsekataan, onko TotVol=0 -rivejÃ¤
        
        a = 0
        zeroRow = next((m for m, x in enumerate(sdata0['TotVol'][a:]) if x==0), None)
        
        if zeroRow is not None:
            
            zeroRow = a + zeroRow
                    
            while (zeroRow is not None) & (zeroRow<len(sdata0)-1): # Ei tarkisteta viimeistÃ¤ riviÃ¤ tÃ¤ssÃ¤, koska tuottaa ongelmia. Tsekataan alempana.
                
                new_vol = tulotdata['volume'+str(int(sdata0['AIKA'][zeroRow]))][0]
                
                if new_vol > 5: # PÃ¤Ã¤tehakkuun tapauksessa 0-rivi on Tulot-tiedostossa vasta kohdassa aika+1
                    new_vol = tulotdata['volume'+str(int(sdata0['AIKA'][zeroRow+1]))][0]
                
                sdata0.loc[zeroRow,'TotVol'] = new_vol
                
                # Haetaan seuraava 0-rivi:
                    
                a = zeroRow + 1    
                zeroRow = next((m for m, x in enumerate(sdata0['TotVol'][a:]) if x==0), None)
                
                if zeroRow is not None:
                    zeroRow = a + zeroRow
                else:
                    break
    
    
    
        # Jos 0-biomassoja, tarkistetaan myÃ¶s 0-iÃ¤t ja korjataan ikÃ¤ vastaamaan AIKA-jaksoja
    
        a = 0
        while a < len(sdata0):
            zeroRow = next((m for m, x in enumerate(sdata0['Age'][a:]) if x==0), None)
            if zeroRow is not None:
                zeroRow = zeroRow+a
                if sdata0['Neulasetkg'][zeroRow]==0: # Muokataan vain, jos neulasissa 0-rivejÃ¤
                    nonZeroRow = next((m for m, x in enumerate(sdata0['Age'][zeroRow:]) if x!=0), None)
                    if nonZeroRow is not None: # Jos 0-ikÃ¤ viimeisellÃ¤ rivillÃ¤, ei tehdÃ¤ tarkistuksia/korjauksia
                        nonZeroRow = nonZeroRow + zeroRow
                        if sdata0['AIKA'][nonZeroRow]-sdata0['AIKA'][zeroRow] != sdata0['Age'][nonZeroRow]-sdata0['Age'][zeroRow]:
    
                            sdata0.loc[zeroRow, 'Age'] = int(sdata0['Age'][nonZeroRow] - (sdata0['AIKA'][nonZeroRow]-sdata0['AIKA'][zeroRow]))
                            with open(outfile, "a") as myfile:
                                myfile.write('\n0-age corrected to: ' + str(sdata0.loc[zeroRow, 'Age']) + ', i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
    
                    else:
                        with open(outfile, "a") as myfile:
                            myfile.write('\n0-age detected at last row!! i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
    
                a = zeroRow + 1
            else:
                a = len(sdata0)+1     
                
                
    
    
        # Tsekataan TotVol ja DomH nollat ja muokataan
        # TÃ¤hÃ¤n tarvii ehkÃ¤ suuremman alarajan?
        # PitÃ¤Ã¤ lisÃ¤tÃ¤ myÃ¶s tyhjien tsekkaus!
        
        for k in range(0,len(sdata0)):
            
            if (sdata0.loc[k,'TotVol']<0.5) | (np.isnan(sdata0.loc[k,'TotVol']) == True): #voi olla 0.01 tilavuuksiakin
                
                if k < len(sdata0)-1:
                    sdata0.loc[k,'TotVol'] = max(0.1, 0.1*sdata0.loc[k+1,'TotVol']) # LisÃ¤tÃ¤Ã¤n totvol nollien tilalle muu pieni luku (0.001 oli liian pieni, kasvu ei lÃ¤htenyt kÃ¤yntiin, joten kasvatettu 0.2:een) paitsi et tulee 0.5
                else:
                    sdata0.loc[k,'TotVol'] = 0.1
                    
            if (sdata0.loc[k,'DomH']<0.3) | (np.isnan(sdata0.loc[k,'DomH']) == True):
                 
                 if k < len(sdata0)-1:
                     sdata0.loc[k,'DomH'] = max(0.3, 0.3*sdata0.loc[k+1,'DomH']) # LisÃ¤tÃ¤Ã¤n domh nollien tilalle muu pieni luku
                 else:
                     sdata0.loc[k,'DomH'] = 0.3       
       
        
        # Adding values to "missing" biomass components with similar increase-% per year as TotVol
        # Neulaset viimeisenÃ¤, jotta sitÃ¤ voidaan kÃ¤yttÃ¤Ã¤ arvioimaan, tarviiko biomassoja tÃ¤ydentÃ¤Ã¤.
        # Poistettu tÃ¤stÃ¤ listasta 'tuotosKuolleet', koska se on kumulatiivinen..
        
        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg', \
                     'TukPL1', 'TukPL2', 'TukPL3', 'TukPL4',\
                         'KuiPL1', 'KuiPL2', 'KuiPL3', 'KuiPL4',\
                         'HukPL1', 'HukPL2', 'HukPL3', 'HukPL4',\
                             'tuotosTilavuus',\
                                 'Runkopuukg', 'Hukkapuukg',\
                                     'MeanDbh','medianH', 'BA', 'tuotos']  
        #comps = ['DomH','MeanDbh','medianH',\
        #          'BA','Kannot_juuretkg',\
        #                  'tuotosTilavuus',\
        #                      'TukPL1', 'TukPL2','TukPL3',\
        #                          'KuiPL1', 'KuiPL2','KuiPL3',\
        #                              'HukPL1', 'HukPL2','HukPL3',\
        #                                  'VolPL1', 'VolPL2', 'VolPL3']  

        #biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg'] #tama lisatty tahan elokuu 23 jotta namakin tsekataan

        # if repairB==1: # Jos lisÃ¤tty 1-rivi  kopioimalla, korjataan sen komponentit ensin
        #     for j in range(0,len(biom_comp)):
        #         comp = biom_comp[j]
        #         sdata0.loc[1, comp] = sdata0['TotVol'][1]/sdata0['TotVol'][2]*sdata0.loc[2, comp]    
            
        # if repair==1: # Jos lisÃ¤tty 0-rivi kopioimalla, korjataan sen komponentit ensin
        #     for j in range(0,len(biom_comp)):
        #         comp = biom_comp[j]
        #         sdata0.loc[0, comp] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, comp]
                
                            
    
                
        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg']
    
        for j in range(0,len(biom_comp)):
        
            comp = biom_comp[j]
            k=0
    
            while k < len(sdata0)-1: 
                
                zeroRow = next((m for m, x in enumerate(sdata0[comp][k:]) if x==0), None)
    
                if zeroRow is not None: # If zeroRow exists
                    
                    zeroRow = zeroRow + k
                    
                    # if sdata0['Neulasetkg'][zeroRow]==0: # Muiden komponenttien numeroita muokatan vain, jos neulasissa on nollaa. TÃ¤mÃ¤ poistettu, koska nÃ¤issÃ¤ komponenteissa voi olla nollaa, vaikka neulasmassa ei olisikaan nollaa..
                    nonZeroRow = next((m for m, x in enumerate(sdata0[comp][zeroRow:]) if x!=0), None)
                    
                    if nonZeroRow is not None:
                        
                        nonZeroRow = nonZeroRow + zeroRow
                        sdata0.loc[zeroRow, comp] = sdata0['TotVol'][zeroRow]/sdata0['TotVol'][nonZeroRow]*sdata0.loc[nonZeroRow, comp] # tää vaikuttaa siihen myös et 0 riville tulee 0.5 vol ja seuraava on jotain
                    
                        k=k+1
                        
                    else:
                        
                        k=k+1                          
    
                else:
                    k=k+1
                            
    
                        
                        # else:
                            
                            # Jotta ei toistu jokaisella komponentilla erikseen, 
                            # nÃ¤mÃ¤ korjataan vasta lopuksi: 
                            # kopioidaan loppuun rivit, jotka vastaavat pÃ¤Ã¤hakkuun jÃ¤lkeistÃ¤ kehitystÃ¤
                            
                            # EtsitÃ¤Ã¤n toisen pÃ¤Ã¤hakkuun jÃ¤lkeinen rivi
            
    
             
                
             
    
    
    
        # TehdÃ¤Ã¤n sarake 'treesp'=pÃ¤Ã¤puulaji volumen perusteella
        
        sdata0['treesp'] = pd.Series([], dtype='int64') #        sdata0['treesp'] = pd.Series([])

        
        for t in range(0,len(sdata0)):
            
            
            volPine = sdata0['VolPL1'][t]
            volSpruce = sdata0['VolPL2'][t]
            volBirch = sdata0['VolPL3'][t] + sdata0['VolPL4'][t]
        
            spe = 3
        
            if volPine >= max(volSpruce, volBirch):
                spe = 1
            
            if volSpruce > max(volPine, volBirch):
                spe = 2    
    
            sdata0.loc[t, 'treesp']=spe
    
                            
                            
        # Tsekataan, lÃ¶ytyykÃ¶ viimeisen pÃ¤Ã¤tehakkuun jÃ¤lkeen pelkkiÃ¤ 0-rivejÃ¤. Jos nÃ¤in,
        # kopioidaan aiemmat pÃ¤Ã¤tehakkuun jÃ¤lkeiset rivit loppuun.
        # Korjataan sen jÃ¤lkeen myÃ¶s AIKA-sarake
    
        if len(paateh)>0:
            paate1 =  next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>paateh[0]), None) 
            nonZeroRow = next((m for m, x in enumerate(sdata0['Neulasetkg'][paate1:]) if x!=0), None)
            
            if len(paateh)>1: # Jos useampi kuin yksi pÃ¤Ã¤tehakkuu
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>=paateh[1]), None) # tÃ¤ssÃ¤ aiemmin x>paateh[1] !!
                nonZeroRow = next((m for m, x in enumerate(sdata0['Neulasetkg'][paaterow:]) if x!=0), None)
            
                if nonZeroRow is None:
                    
                    paate2aika = sdata0['AIKA'][paaterow]
                    
                    b = 0
                    
                    for p in range(paaterow,len(sdata0)):
                        print(p, paate2aika, "paatehakkuun jalkeiset nollarivit")
    
                        currow = paate1 + p - paaterow -1 # LisÃ¤tty tÃ¤hÃ¤n -1, koska rivillÃ¤ 617 lisÃ¤tty x>=paateh[1]
    
                        if currow<b: # LisÃ¤ys 17.10: Joissakin tapauksissa yllÃ¤ toimii -1, toisissa 0, joten muokataan currow tarvittaessa...
                            currow=b
                            
                        sdata0.loc[p]=sdata0.loc[currow]
                        
                        if b==0:
                            sdata0.loc[p, 'AIKA'] = int(paate2aika)
                            sdata0.loc[p, 'tuotosKuolleet'] = sdata0.loc[p-1, 'tuotosKuolleet']
                            b=b+1
                        else:
                            sdata0.loc[p, 'AIKA'] = int(paate2aika + sdata0.loc[currow, 'AIKA'] - sdata0.loc[currow-b, 'AIKA'])
                            sdata0.loc[p, 'tuotosKuolleet'] = sdata0.loc[p-1, 'tuotosKuolleet']
                            b=b+1
                
            else:
                # Jos vain yksi pÃ¤Ã¤tehakkuu, ei voida kopioida rivejÃ¤ aiemman pÃ¤Ã¤tehakkuun jÃ¤lkeen..
                print('MitÃ¤s tÃ¤hÃ¤n? kuvio:', str(kuviolista[i]),paateh)   
                if nonZeroRow is None:
                    with open(outfile, "a") as myfile:
                        myfile.write('\nPÃ¤Ã¤tehakkuun jÃ¤lkeisiÃ¤ 0-rivejÃ¤ ei saada korjattua! i=' + str(i)+ ' kuvio:'+ str(kuviolista[i])+ ' paateh: '+str(paateh)) #Ã¤Ã¤ 
    
        
        # Tsekataan, lÃ¶ytyykÃ¶ viimeiseltÃ¤ riviltÃ¤ (aika=100), nollariviÃ¤.
        # Jos lÃ¶ytyy, poistetaan.
        
    
        s100 = sdata0[(sdata0['AIKA']==100) & (sdata0['Age']==0)] # Vai vol=0?
        
        if len(s100)>0:
            sdata0 = sdata0[sdata0['AIKA']!=100]
            
        # -------------------------------------------------------------------
        # LisÃ¤ys 13.9.2021
        # Tarkistetaan, olisiko ensimmÃ¤inen Motti-tiedosto 1-rivin mittainen:
        
        vuodet_p = vuodet[vuodet>0]
        
        if (len(vuodet_p)>0):
            
            sdata = sdata0[sdata0['AIKA'] < vuodet_p[0]]
            
            if len(sdata)==1:
                #sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
                nolla2_df=pd.DataFrame([sdata0.iloc[0]])
                sdata0 = pd.concat((sdata0,nolla2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = vuodet_p[0]-1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                sdata0.loc[1,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][1]))][0]
                sdata0.loc[1,'Age'] = sdata0['Age'][0]+(sdata0['AIKA'][1]-sdata0['AIKA'][0])
    
                for c in range(0,len(comps)):
                    if sdata0.loc[0,comps[c]]>0:
                        # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
                        # sitÃ¤ voisi supistaa..??
                        sdata0.loc[1,comps[c]] = sdata0['TotVol'][1]/sdata0['TotVol'][0]*sdata0.loc[0, comps[c]]
    
                for c in range(0,len(biom_comp)):
                    if sdata0.loc[0,biom_comp[c]]>0:
                        # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
                        # sitÃ¤ voisi supistaa..??
                        sdata0.loc[1,biom_comp[c]] = sdata0['TotVol'][1]/sdata0['TotVol'][0]*sdata0.loc[0, biom_comp[c]]
    
    
        # -------------------------------------------------------------------
        # LisÃ¤ys 9.3.2022
        # Tarkistetaan, olisiko *toinen* Motti-tiedosto 1-rivin mittainen:
        
        #if (len(vuodet_p)>0):
        #    
        #    if len(vuodet_p)==1:
        #        end_year=101
        #    else:
        #        end_year = vuodet_p[1]
        #        
        #    sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[0]) & (sdata0['AIKA'] < end_year)] # Huom. tÃ¤ssÃ¤ ei saa resetoida indeksiÃ¤, jotta jatko toimii
        #    
        #    if len(sdata)==1:
        #        sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # LisÃ¤tÃ¤Ã¤n rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
        #        sdata0['AIKA'][len(sdata0)-1] = end_year - 1
        # 
        #         sdata0 = sdata0.sort_values(by=['AIKA'])
        #         sdata0 = sdata0.reset_index(drop=True)
        #         
        #        sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
        #        
        #        sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
        #        
        #        if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa riviÃ¤
        #            sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
        #
        #            for c in range(0,len(comps)):
        #                if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
        #                    # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
        #                    # sitÃ¤ voisi supistaa..??
        #                    sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
        # 
        #            for c in range(0,len(biom_comp)):
        #                if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
        #                    # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
        #                    # sitÃ¤ voisi supistaa..??
        #                    sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
        
        
        # -------------------------------------------------------------------
        # Update 6.3.2023: checking all 1-row-motti-files except the last one
        
        for v in range(1, len(vuodet_p)):
            
            if (len(vuodet_p)>0):
                
                if len(vuodet_p)==1:
                    end_year=101
                else:
                    end_year = vuodet_p[v]                
                
                sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[v-1]) & (sdata0['AIKA'] < end_year)] # Huom. tÃ¤ssÃ¤ ei saa resetoida indeksiÃ¤ jotta jatko toimii
                
                if len(sdata)==1:
                    #sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # LisÃ¤tÃ¤Ã¤n rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
                    sdataind_df=pd.DataFrame([sdata0.iloc[sdata.index[0]]])
                    sdata0 = pd.concat((sdata0,sdataind_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                    sdata0['AIKA'][len(sdata0)-1] = end_year - 1
                    
                    sdata0 = sdata0.sort_values(by=['AIKA'])
                    sdata0 = sdata0.reset_index(drop=True)
                    
                    sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
                    
                    sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
                    
                    if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa rivia
                        sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
                        
                        for c in range(0,len(comps)):
                            if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
                                # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                                # sita voisi supistaa..??
                                sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
                        
                        for c in range(0,len(biom_comp)):
                            if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
                                # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                                # sita voisi supistaa..??
                                sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
        
        # Update 6.3.2023: Checking if the last motti-file would be a 1-row-file
        
        if (len(vuodet_p)>1):
            
            sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[len(vuodet_p)-1])] # Huom. tassa ei saa resetoida indeksia, jotta jatko toimii
            end_year = 101
            
            if (len(sdata)==1) & (sdata0['AIKA'][len(sdata0)-1]<100):
                
                #sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # Lisataan rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
                sdataind2_df=pd.DataFrame([sdata0.iloc[sdata.index[0]]])
                sdata0 = pd.concat((sdata0, sdataind2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                    
                sdata0['AIKA'][len(sdata0)-1] = end_year - 1
                
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
                
                sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
                
                if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa rivia
                    sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
                    
                    for c in range(0,len(comps)):
                        if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
                            # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                            # sita voisi supistaa..??
                            sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
                    
                    for c in range(0,len(biom_comp)):
                        if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
                            # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                            # sita voisi supistaa..??
                            sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
                            
                
        
        # -------------------------------------------------------------------
        # LisÃ¤ys 9.3.2022
        # Tarkistetaan, tuleeko perÃ¤kkÃ¤isten harvennusten takia liian paljon
        # vÃ¤heneviÃ¤ tilavuuksia 
        
        comps_notk = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg', \
                             'tuotosTilavuus','MeanDbh','medianH', 'BA', 'tuotos', 'TotVol']  
        
        
        if pera >= 0:
            
            vuodet_p = vuodet[vuodet>=0] # LisÃ¤ys 17.10.2022: vuodet>=0 (jos tÃ¤ssÃ¤ olisi vuodet>0, vuodet_p saattaisi jÃ¤Ã¤dÃ¤ tyhjÃ¤ksi -> virhe seuraavalla rivillÃ¤)
            
            sdata = sdata0[sdata0['AIKA'] < vuodet_p[0]]
            
            if (len(sdata)<=2) & (len(sdata)>0): # LisÃ¤ys 17.10.2022: len(sdata)>0, koska jos sdata on tyhjÃ¤ -> error
                
                if (sdata['TotVol'][0] - sdata['TotVol'][1]) > 0:
                    
                    sdata0.loc[0]=sdata0.loc[1]
                    sdata0.loc[0, 'AIKA'] = sdata0['AIKA'][1] -1
                    sdata0.loc[0, 'Age'] = sdata0['Age'][1] -1
                    
                    
                    for c in range(0,len(comps_notk)):
                        
                        comp = comps_notk[c]
                        sdata0.loc[0, comp] = sdata0.loc[0, comp] * 0.97 # Arvio, jotta 
                                            
    
    
        
        
        # Muutetaan Age=-1 --> Age=0, nÃ¤itÃ¤ lÃ¶ytyi suoraan Motti-datasta
        
        sdata0['Age'] = np.clip(sdata0['Age'], 0, np.nan) 
        
        # -------------------------------------------------------------------
        # Lisays 28.3.2023:
        # Tarkistetaan, loytyyko iasta perakkaisia nollia
        
        a = 0
        while a < len(sdata0):
            zeroRow = next((m for m, x in enumerate(sdata0['Age'][a:]) if x==0), None)
            
            if zeroRow is not None:
                zeroRow = zeroRow+a
                try:
                    if sdata0['Age'][zeroRow+1]==0:
                        sdata0.loc[zeroRow, 'Age'] = 1
                        sdata0.loc[zeroRow+1, 'Age'] = 1 + int(sdata0['AIKA'][zeroRow+1]-sdata0['AIKA'][zeroRow])
                        with open(outfile, "a") as myfile:
                            myfile.write('\nDouble zero ages!! i=' + str(i))
                        
                    a = zeroRow + 1
                except:
                    a = zeroRow + 1
                        
            else:
                a = len(sdata0)+1     
        
        
        
        
        # Varsinaisten Motti-tiedostojen teko alkaa tÃ¤stÃ¤
        
        if len(sdata0)>0:
            
            if len(vuodet)>0:
                
                for j in range(0,n_mottifiles):
                
                    if j == 0:
                        sdata = sdata0[sdata0['AIKA'] < vuodet[j]]
                        if len(sdata)==0:
                            continue
                    
                    if (j>0) & (j<n_mottifiles-1):
                        sdata = sdata0[(sdata0['AIKA'] >= vuodet[j-1]) & (sdata0['AIKA'] < vuodet[j])]
                        sdata = sdata.reset_index(drop=True)
                        if len(sdata)==0:
                            continue
                    
                    if j == n_mottifiles-1:
                        sdata = sdata0[sdata0['AIKA'] >= vuodet[j-1]]
                        sdata = sdata.reset_index(drop=True)
                        if len(sdata)==0:
                            continue
                    
                    kasvatus = pd.Series(np.ones(len(sdata)))
                    vuosi = sdata['AIKA']
                    ika = sdata['Age']
                    N = sdata['N']
                    PPA = sdata['BA']
                    Hg = sdata['medianH']
                    Dg = sdata['MeanDbh']
                    Hdom = sdata['DomH']
                    Tilavuus = sdata['TotVol']
                    Tukki = sdata['TukPL1']+sdata['TukPL2']+sdata['TukPL3']+sdata['TukPL4']
                    Kuitu = sdata['KuiPL1']+sdata['KuiPL2']+sdata['KuiPL3']+sdata['KuiPL4']
                    Hukka = sdata['HukPL1']+sdata['HukPL2']+sdata['HukPL3']+sdata['HukPL4']
                    Tuotos = sdata['tuotosKuolleet'] + sdata['tuotosTilavuus']
                    Kuolleisuus = Tuotos*0. # Ei kÃ¤ytetÃ¤ Susissa
                    runko_aines = sdata['Runkopuukg']/1000. # YksikkÃ¶ tn
                    runko_hukka = sdata['Hukkapuukg']/1000.
                    elavat_oksat = sdata['Oksatkg']*0.8/1000. # Oksatkg = elÃ¤vÃ¤t ja kuolleet oksat yhteensÃ¤
                    kuolleet_oksat = sdata['Oksatkg']*0.2/1000.
                    lehdet = sdata['Neulasetkg']/1000.
                    kannot = sdata['kannot']/1000.
                    juuret = sdata['>2mm_juuretkg']/1000.
                    hienojuuret = sdata['h_juuretkg']/1000.
                    
                    # Tsekkaa, ettÃ¤ tÃ¤mÃ¤ on ok!!!
                    for t in range(0,len(sdata)):
                        if t==0:
                            Tuotos[0]=sdata['tuotosTilavuus'][0]
                        
                        else:
                            Tuotos[t] = sdata['tuotosTilavuus'][t] + (sdata['tuotosKuolleet'][t] - sdata['tuotosKuolleet'][t-1])
                    
                    # ---tÃ¤hÃ¤n asti
                    
                    data = pd.DataFrame(np.transpose(np.array([kasvatus, vuosi, ika, N, PPA, Hg, Dg, Hdom, Tilavuus, Tukki, Kuitu, \
                                         Hukka, Tuotos, Kuolleisuus, runko_aines, runko_hukka, \
                                         elavat_oksat, kuolleet_oksat, lehdet, kannot, juuret, hienojuuret])), \
                                         columns=['Kasvatus', 'Vuosi', 'Ikä', 'N', 'PPA', 'Hg', 'Dg', 'Hdom', 'Tilavuus', 'Tukki', 'Kuitu',\
                                                  'Hukka', 'Tuotos', 'Kuolleisuus', 'runko(aines)', 'runko(hukka)', \
                                                  'elävät oksat', 'kuolleet oksat', 'lehdet', 'Kannot', 'Juuret >2mm', 'Hienojuuret'])
                    
                    standInfo = kuviot[kuviot['KUVIO']==kuviolista[i]]
                    standInfo = standInfo.reset_index(drop=True)
                    
                    """ TSEKKAA, ETTÃ„ PUULAJI TULEE OIKEIN! """
                    
                    
                    
                    if (project == 'life') | (project == 'hiilipolku'):
                        spe = sdata['treesp'][len(sdata)-1]
                    
                    else: 
                        if (project =='vmi') | (project=='suo'):
                            
                            spe = int(sdata['treesp'][len(sdata)-1]) # viimeisen rivin puulaji
                            # treesp=kdata['Selite'][0][28:33]
                            # spe = 3
                            
                            # if treesp=='MÃ¤nty':
                            #     spe = 1
                            # if treesp=='Kuusi':
                            #     spe = 2       
                                
                            
                        else:
                            
                            treesp = pdata2['Selite'][0][28:33]
                            
                            spe = 3
                            
                            if treesp=='Mänty': #onk tÃ¤ a va Ã¤#
                                spe = 1
                            if treesp=='Kuusi':
                                spe = 2 
                    
                    
                    # print(kuviolista[i], ': puulaji', spe)
                    with open(outfile, "a") as myfile:
                        myfile.write('\n' + str(kuviolista[i]) +  ': puulaji ' + str(spe)+ ' n'+str(j))
                    
                    
                    fout = mottifolder + '/' + str(kuviolista[i]) + '_n' + str(j) + '.xls'
                    
                    cols = pd.DataFrame([['','','','',spe,'','','','','','','','','']],columns=['Kasvatus',	'Vuosi',	'id Harvennus',	'Harvennus',	'id Puulaji',	'Puulaji',	'Tukki[m³/ha]',	'Pikkutukki[m³/ha]',	'Kuitu[m³/ha]',	'Energiapuu, runko(aines)[m³/ha]',	'Energiapuu, runko(hukka)[m³/ha]',	'Energiapuu, oksat(e)[m³/ha]',	'Energiapuu, oksat(k)[m³/ha]',	'Energiapuu, kannot ja juuret[m³/ha]'])
                    
                    writer = pd.ExcelWriter(fout, engine = 'xlsxwriter')
                    data.to_excel(writer, index=False, sheet_name = 'Puustotunnukset')
                    cols.to_excel(writer, index=False, sheet_name = r'Kertymät')
                    #writer.save()
                    writer.close()
            else:
                # print(str(kuviolista[i]) + ': poistumat puuttuvat! i=' + str(i))
                with open(outfile, "a") as myfile:
                    myfile.write('\n' + str(kuviolista[i]) + ': poistumat puuttuvat! i=' + str(i))
                
                sdata = sdata0.copy()
                
                kasvatus = pd.Series(np.ones(len(sdata)))
                vuosi = sdata['AIKA']
                ika = sdata['Age']
                N = sdata['N']
                PPA = sdata['BA']
                Hg = sdata['medianH']
                Dg = sdata['MeanDbh']
                Hdom = sdata['DomH']
                Tilavuus = sdata['TotVol']
                Tukki = sdata['TukPL1']+sdata['TukPL2']+sdata['TukPL3']+sdata['TukPL4']
                Kuitu = sdata['KuiPL1']+sdata['KuiPL2']+sdata['KuiPL3']+sdata['KuiPL4']
                Hukka = sdata['HukPL1']+sdata['HukPL2']+sdata['HukPL3']+sdata['HukPL4']
                Tuotos = sdata['tuotos']
                Kuolleisuus = Tuotos*0. # Ei kÃ¤ytetÃ¤ Susissa
                runko_aines = sdata['Runkopuukg']/1000. # YksikkÃ¶ tn
                runko_hukka = sdata['Hukkapuukg']/1000.
                elavat_oksat = sdata['Oksatkg']*0.8/1000. # Oksatkg = elÃ¤vÃ¤t ja kuolleet oksat yhteensÃ¤
                kuolleet_oksat = sdata['Oksatkg']*0.2/1000.
                lehdet = sdata['Neulasetkg']/1000.
                kannot = sdata['kannot']/1000.
                juuret = sdata['>2mm_juuretkg']/1000.
                hienojuuret = sdata['h_juuretkg']/1000.
                
                data = pd.DataFrame(np.transpose(np.array([kasvatus, vuosi, ika, N, PPA, Hg, Dg, Hdom, Tilavuus, Tukki, Kuitu, \
                                     Hukka, Tuotos, Kuolleisuus, runko_aines, runko_hukka, \
                                     elavat_oksat, kuolleet_oksat, lehdet, kannot, juuret, hienojuuret])), \
                                     columns=['Kasvatus', 'Vuosi', 'Ikä', 'N', 'PPA', 'Hg', 'Dg', 'Hdom', 'Tilavuus', 'Tukki', 'Kuitu',\
                                              'Hukka', 'Tuotos', 'Kuolleisuus', 'runko(aines)', 'runko(hukka)', \
                                              'elävät oksat', 'kuolleet oksat', 'lehdet', 'Kannot', 'Juuret >2mm', 'Hienojuuret'])
                
                standInfo = kuviot[kuviot['KUVIO']==kuviolista[i]]
                standInfo = standInfo.reset_index(drop=True)
                
                
                """ TSEKKAA, ETTÃ„ PUULAJI TULEE OIKEIN! """
                if (project == 'life') | (project == 'hiilipolku'):
                    spe = sdata['treesp'][len(sdata)-1]
                
                else:
                    
                    if project =='vmi':
                        
                        # spe = int(sdata['treesp'][len(sdata)-1]) # viimeisen rivin puulaji
                       
                        treesp=kdata['Selite'][0][28:33]
                        spe = 3
                        
                        if treesp=='Mänty': #Ã¤
                            spe = 1
                        if treesp=='Kuusi':
                            spe = 2    
                    else:
                        treesp = pdata2['Selite'][0][28:33]
                        
                        spe = 3
                        
                        if treesp=='Mänty':
                            spe = 1
                        if treesp=='Kuusi':
                            spe = 2 
                    
                # print(kuviolista[i], ': puulaji ', spe)
                with open(outfile, "a") as myfile:
                    myfile.write('\n' + str(kuviolista[i]) +  ':' + str(spe)+ ' n0')
                
                
                
                fout = mottifolder + '/' + str(kuviolista[i]) + '_n0.xls'
                
                cols = pd.DataFrame([['','','','',spe,'','','','','','','','','']],columns=['Kasvatus',	'Vuosi',	'id Harvennus',	'Harvennus',	'id Puulaji',	'Puulaji',	'Tukki[m³/ha]',	'Pikkutukki[m³/ha]',	'Kuitu[m³/ha]',	'Energiapuu, runko(aines)[m³/ha]',	'Energiapuu, runko(hukka)[m³/ha]',	'Energiapuu, oksat(e)[m³/ha]',	'Energiapuu, oksat(k)[m³/ha]',	'Energiapuu, kannot ja juuret[m³/ha]'])
                
                writer = pd.ExcelWriter(fout, engine = 'xlsxwriter')
                data.to_excel(writer, index=False, sheet_name = 'Puustotunnukset')
                cols.to_excel(writer, index=False, sheet_name = r'Kertymät')
                #writer.save()
                writer.close()


def create_motti_files_silvi_lp4_vol001(project, scen, maku=None):       #alkutilavuus, annetaan olla 0.01tai 0.01*Totvol[1]lp4 skaalaukseen otettu lp4 mukaan  # Huom! ohje>1
    '''
    Modified: January 2023
        
    Create Motti input files for Susi simulations from "Skene-Motti" results.
    
    Works technically, but contains many small issues that should be fixed.
    
    project = 'life' (Hydrology Life -project) or 'suo' (SUO-project)
    project-parameter affects only main tree species (and column names)!
    
    
    Modified 8.3.2022:
        - korjattu pÃ¤Ã¤tehakkuuseen liittyviÃ¤ ongelmia; hdom pitÃ¤isi nyt saada
        jÃ¤rkeviÃ¤ arvoja, biomassojen ei pitÃ¤isi jÃ¤Ã¤dÃ¤ pÃ¤Ã¤tehakkuun jÃ¤lkeen nollaksi
        (ainakaan toisen pÃ¤Ã¤tehakkuun jÃ¤lkeen)
        
    Issues to fix:
        - Jos pÃ¤Ã¤tehakkuu tulee vuonna 2, pitÃ¤Ã¤ vuosille 0 ja 1 hakea muut komponentit
        totvol perusteella: hdom saa tÃ¤ssÃ¤ liian suuria arvoja. Tuotos jÃ¤Ã¤ nollaksi?!
        - See comments below...
        
    
    '''
    
    inputmainfol=r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan_vol001/'
    mottifolder = inputmainfol+scen[:-6]+'/'+scen+'/'
    if not os.path.exists(inputmainfol):
        os.mkdir(inputmainfol)
    if not os.path.exists(inputmainfol+scen[:-6]+'/'):
        os.mkdir(inputmainfol+scen[:-6]+'/')
    if not os.path.exists(mottifolder):
        os.mkdir(mottifolder)
            
    outfile = mottifolder + '/out.txt'
    
    project = 'hiilipolku'
    
    kuviot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', encoding='latin1', sep=';')
    puustot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_puustot.csv', encoding='latin1', sep=';')
    poistumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_poistumat.csv', encoding='latin1', sep=';')
    tapahtumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_tapahtumat.csv', encoding='latin1', sep=';', index_col=False)
    tulot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kasvut.csv', encoding='latin1', sep=';')
    
    
    if (project == 'life') | (project =='vmi') | (project =='hiilipolku'):
        kuviot = kuviot.rename(columns={'kuvio':'KUVIO'})
        puustot = puustot.rename(columns={'kuvio':'KUVIO',
                                          '_2mm_juuretkg':'>2mm_juuretkg',})
        poistumat = poistumat.rename(columns={'kuvio':'KUVIO'})
        tapahtumat = tapahtumat.rename(columns={'kuvio':'KUVIO',
                                      '_2003': '2003',
                                      '_2004': '2004',
                                      '_2005': '2005'})
        tulot = tulot.rename(columns={'kuvio':'KUVIO'})
    
    
    kuviolista = list(set(kuviot['KUVIO'][kuviot['kuivatustilanne']>=3])) #tassa valitaan 3 Ojitettu kangas ja muut?DRAINAGESTATE	, ojitettu kangas jos pÃ¤Ã¤ryhmÃ¤1 ja alaryhmÃ¤ 3 ja kuivatustilanne 5?? 1 Ojittamaton kangas, 2 Soistunut kangas,3 Ojitettu kangas,	6 luonnontilainen suo,	7 Ojikko, 8	Muuttuma,9	Turvekangas
    
    
    for i in range(0,len(kuviolista)):
    
        print(i, 'stand: ', kuviolista[i])
        
        # alta poistettu ohjeeseen ja skenaarioon viittaavat rajaukset,
        # oletetaan, ettÃ¤ hiilipolussa vain yksi ohje per kuvio
        
        sdata0 = puustot[(puustot['KUVIO']==kuviolista[i])]
        pdata = poistumat[(poistumat['KUVIO']==kuviolista[i]) & (poistumat['HARV']>1)] #mita harv koodit meinaa?
        pdata2 = poistumat[(poistumat['KUVIO']==kuviolista[i])]
        tdata = tapahtumat[(tapahtumat['KUVIO']==kuviolista[i])]
        tulotdata = tulot[(tulot['KUVIO']==kuviolista[i])]
        kdata = kuviot[kuviot['KUVIO']==kuviolista[i]]
        
        sdata0 = sdata0.reset_index(drop=True)
        pdata = pdata.reset_index(drop=True)
        pdata2 = pdata2.reset_index(drop=True)
        tdata = tdata.reset_index(drop=True)
        tulotdata = tulotdata.reset_index(drop=True)
        kdata = kdata.reset_index(drop=True)
        
        tdata = tdata[(tdata['2003']>-1) | (tdata['2004']>-1) | (tdata['2005']>-1)]
        
        ensih = tdata[tdata['2003']>-1]['2003'].values
        harv = tdata[tdata['2004']>-1]['2004'].values
        paateh = tdata[tdata['2005']>-1]['2005'].values
        #voisi kayda tsekkaamassa onko taimikonhoito tai varhaisperkaus niin sitten ei muuttaisi puulajia tai ottaisi sen muualta et sailyy kuusena tai mantyna? 
        #tdic={'_2001':'lannoitus','_2002':'kunnostusojitus','_2003':'ensiharvennus','_2004':'harvennus','_2005':'paatehakkuu','_2006':'varhaiperkaus','_2007':'taimikonhoito','_2008':'viljely','_2009':'luontainen uudistuminen'}
        # taimh = tdata[tdata['_2007']>-1]['_2007'].values
        # Pudotetaan viimeisen vuoden tapahtumat pois, koska sekoittavat
        # myÃ¶hemmin tulevan koodin
        ensih = ensih[ensih<100]
        harv = harv[harv<100]
        paateh = paateh[paateh<100]
        
        vuodet = np.unique(np.sort(np.concatenate((ensih, harv, paateh))))
        
        ap = 0
        
        pera = -1
        
        # Tsekataan, onko tapahtumissa perÃ¤kkÃ¤isiÃ¤ vuosia ja poistetaan niistÃ¤ jÃ¤lkimmÃ¤inen
        for v in range(0,len(vuodet)-1):
            if vuodet[v+1]-vuodet[v]==1:
                
                ensih = np.delete(ensih, np.where(ensih == vuodet[v+1]))
                harv = np.delete(harv, np.where(harv == vuodet[v+1]))
                paateh = np.delete(paateh, np.where(paateh == vuodet[v+1]))
                
                vuodet = np.delete(vuodet, [v+1], None)
                
                pera = v
                
                ap = ap + 1
                
            if v+1>=len(vuodet)-1:
                break
        
        if ap>1:
            with open(outfile, "a") as myfile:
                myfile.write('\n\nKaksi kertaa perakkaiset vuodet! stand=' + str(kuviolista[i]) + ', i=' + str(i))
    
        n_mottifiles = len(vuodet) + 1
    
    
        if len(sdata0)==0: # skipataan kuvio, jos dataa ei ole
            continue
    
        if len(sdata0)>21:
            print('sdata0: liikaa rivejÃ¤?!')
        
    
                                    
        # Luodaan pÃ¤Ã¤tehakkuun jÃ¤lkeiset rivit ja korjataan rivin kopioimisesta aiheutuneet
        # luvut
        
        # muokattu 6.8.2021: poistettu tuotos, tuotosKuolleet, tuotosTilavuus comps:sta
        comps = ['DomH','MeanDbh','medianH',\
                          'BA','Kannot_juuretkg',\
                              'tuotosTilavuus',\
                              'TukPL1', 'TukPL2','TukPL3','TukPL4',\
                                  'KuiPL1', 'KuiPL2','KuiPL3','KuiPL4',\
                                      'HukPL1', 'HukPL2','HukPL3','HukPL4',\
                                          'VolPL1', 'VolPL2', 'VolPL3','VolPL4',]  #lisatty 091023 PL4 ositteet mukaan

        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg'] #tama lisatty tahan elokuu 23 jotta namakin tsekataan

        for h in range(0,len(paateh)):
            
            sdata_paateh= sdata0[sdata0['AIKA']==paateh[h]].reset_index(drop=True)
            
            go = False
            
            if len(sdata_paateh)>0:
                if sdata_paateh.loc[0, 'Age']>5:
                    go = True
        
            if (len(sdata_paateh)==0) | (go==True): # Luodaan rivit vain, jos ensimmÃ¤istÃ¤ riviÃ¤ ei ole mukana
                
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>paateh[h]), None)
                
                #sdata0 = sdata0.append([sdata0.iloc[paaterow]], ignore_index=True) # Rivi menee taulukon loppuun
                paaterow_df=pd.DataFrame([sdata0.iloc[paaterow]])
                sdata0 = pd.concat((sdata0,paaterow_df),axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = paateh[h]+1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==paateh[h]+1), None)
                
                # sdata0.loc[paaterow,'Age'] = 1 # TÃ¤mÃ¤ korvattu alla olevalla
                
                if paaterow>0:
                    if sdata0.loc[paaterow-1, 'Age']<5:
                        sdata0.loc[paaterow,'Age'] = sdata0.loc[paaterow-1, 'Age']+1
                    else:
                        sdata0.loc[paaterow,'Age'] = sdata0.loc[paaterow+1, 'Age']-1
                
                sdata0.loc[paaterow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][paaterow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[paaterow,comps[c]]>0:
                        sdata0.loc[paaterow,comps[c]] = sdata0['TotVol'][paaterow]/sdata0['TotVol'][paaterow+1]*sdata0.loc[paaterow+1, comps[c]]
                        print(kuviolista[i],comps[c], "rivi 303 ongelma")

            
                if sdata0['tuotosKuolleet'][paaterow] < sdata0['TotVol'][paaterow]:
                    sdata0.loc[paaterow,'N'] = (1 + sdata0['tuotosKuolleet'][paaterow]/sdata0['TotVol'][paaterow])*sdata0['N'][paaterow]        
                    
                    if paaterow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][paaterow] > sdata0['N'][paaterow-1]:
                            if paaterow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[paaterow,'N'] = (sdata0['N'][paaterow-1] + sdata0['N'][paaterow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[paaterow,'N'] = 0.99*sdata0['N'][paaterow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.
                            
                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[paaterow,biom_comp[c]]>0: 
                        sdata0.loc[paaterow,biom_comp[c]] = sdata0['TotVol'][paaterow]/sdata0['TotVol'][paaterow+1]*sdata0.loc[paaterow+1, biom_comp[c]] ####
                
            
        # Luodaan ensiharvennuksen jÃ¤lkeiset rivit (jos vuosi ei osu kohdallaan muuten)
        # EsimerkissÃ¤ nÃ¤ytti siltÃ¤, ettÃ¤ ensiharvennusvuoden riville tulostuu harventamaton tilavuus! Siksi tÃ¤ssÃ¤ kÃ¤sitelty riviÃ¤+1
        # Toisessa esimerkissÃ¤ tulostui harvennettu tilavuus! Joten vaihtelee.
        
        for h in range(0,len(ensih)):
    
    
            ensihrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==ensih[h]+1), None)
            ensihrow2 = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>ensih[h]), None)
            
            if ensihrow is None:
                
                #sdata0 = sdata0.append([sdata0.iloc[ensihrow2]], ignore_index=True) # Rivi menee taulukon loppuun
                ensihrow2_df=pd.DataFrame([sdata0.iloc[ensihrow2]])
                sdata0 = pd.concat((sdata0,ensihrow2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = ensih[h]+1
                
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                ensihrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==ensih[h]+1), None)
                
                # Jos ikÃ¤ lÃ¶ytyy riviÃ¤ alaempaa:
                if sdata0.loc[ensihrow+1,'Age']>0:
                    sdata0.loc[ensihrow,'Age'] = sdata0.loc[ensihrow+1,'Age']-(sdata0['AIKA'][ensihrow+1]-sdata0['AIKA'][ensihrow])
                # Muussa tapauksessa ikÃ¤ lÃ¶ytyy riviÃ¤ ylempÃ¤Ã¤:
                else:
                    if ensihrow>0:
                        sdata0.loc[ensihrow,'Age'] = sdata0.loc[ensihrow-1,'Age']+(sdata0['AIKA'][ensihrow]-sdata0['AIKA'][ensihrow-1])
                    else:
                        # Jos ei lÃ¶ydy kummaltakaan puolelta:
                        with open(outfile, "a") as myfile:
                            myfile.write('\nEnsiharvennuksen jÃ¤lkeistÃ¤ ikÃ¤Ã¤ ei lÃ¶ydy! kuvio=' + str(kuviolista[i]) + ' i:'+ str(i))
                        
                    
                sdata0.loc[ensihrow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][ensihrow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[ensihrow,comps[c]]>0:
                        sdata0.loc[ensihrow,comps[c]] = sdata0['TotVol'][ensihrow]/sdata0['TotVol'][ensihrow+1]*sdata0.loc[ensihrow+1, comps[c]]
                
                # PÃ¤ivitetÃ¤Ã¤n N
                if sdata0['tuotosKuolleet'][ensihrow] < sdata0['TotVol'][ensihrow]:
                    sdata0.loc[ensihrow,'N'] = (1 + sdata0['tuotosKuolleet'][ensihrow]/sdata0['TotVol'][ensihrow])*sdata0['N'][ensihrow]        
    
                    if ensihrow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][ensihrow] > sdata0['N'][ensihrow-1]:
                            if ensihrow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[ensihrow,'N'] = (sdata0['N'][ensihrow-1] + sdata0['N'][ensihrow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[ensihrow,'N'] = 0.99*sdata0['N'][ensihrow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.

                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[ensihrow,biom_comp[c]]>0: 
                        sdata0.loc[ensihrow,biom_comp[c]] = sdata0['TotVol'][ensihrow]/sdata0['TotVol'][ensihrow+1]*sdata0.loc[ensihrow+1, biom_comp[c]] ####
                            
        
        # Luodaan harvennuksen jÃ¤lkeiset rivit (jos vuosi ei osu kohdalleen muuten)
        # Harvennusvuoden rivillÃ¤ on harventamaton tilavuus. Siksi kÃ¤sitelty riviÃ¤+1.
        
        for h in range(0,len(harv)):
    
    
            hrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==harv[h]+1), None)
            hrow2 = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>harv[h]), None)
            
            if hrow is None:
                
                #sdata0 = sdata0.append([sdata0.iloc[hrow2]], ignore_index=True) # Rivi menee taulukon loppuun
                hrow2_df=pd.DataFrame([sdata0.iloc[hrow2]])
                sdata0 = pd.concat((sdata0,hrow2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = harv[h]+1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                hrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==harv[h]+1), None)
                
                # Jos ikÃ¤ lÃ¶ytyy riviÃ¤ alaempaa:
                if sdata0.loc[hrow+1,'Age']>0:
                    sdata0.loc[hrow,'Age'] = sdata0.loc[hrow+1,'Age']-(sdata0['AIKA'][hrow+1]-sdata0['AIKA'][hrow])
                # Muussa tapauksessa ikÃ¤ lÃ¶ytyy riviÃ¤ ylempÃ¤Ã¤:
                else:
                    if hrow>0:
                        sdata0.loc[hrow,'Age'] = sdata0.loc[hrow-1,'Age']+(sdata0['AIKA'][hrow]-sdata0['AIKA'][hrow-1])
                    else:
                    # Jos ei lÃ¶ydy kummaltakaan puolelta:
                        with open(outfile, "a") as myfile:
                            myfile.write('\nHarvennuksen jÃ¤lkeistÃ¤ ikÃ¤Ã¤ ei lÃ¶ydy! i=' + str(i) + ' kuvio:'+ str(kuviolista[i]))
                
                    
                sdata0.loc[hrow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][hrow]))][0]
                
                
                for c in range(0,len(comps)):
                    if sdata0.loc[hrow,comps[c]]>0:
                        sdata0.loc[hrow,comps[c]] = sdata0['TotVol'][hrow]/sdata0['TotVol'][hrow+1]*sdata0.loc[hrow+1, comps[c]]

                        
                # Korjataan N:
                if sdata0['tuotosKuolleet'][hrow] < sdata0['TotVol'][hrow]:
                    sdata0.loc[hrow,'N'] = (1 + sdata0['tuotosKuolleet'][hrow]/sdata0['TotVol'][hrow])*sdata0['N'][hrow]
                    
                    if hrow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][hrow] > sdata0['N'][hrow-1]:
                            if hrow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[hrow,'N'] = (sdata0['N'][hrow-1] + sdata0['N'][hrow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[hrow,'N'] = 0.99*sdata0['N'][hrow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.
                            
                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[hrow,biom_comp[c]]>0: 
                        sdata0.loc[hrow,biom_comp[c]] = sdata0['TotVol'][hrow]/sdata0['TotVol'][hrow+1]*sdata0.loc[hrow+1, biom_comp[c]] ####
    
        # Tsekataan, ettei ole ekan rivin TotVol pyÃ¶ristyksen kanssa ongelmaa:
        
        if (sdata0['TotVol'][0]>10.) & (-0.1 < sdata0['TotVol'][1]-sdata0['TotVol'][0]<0):
            sdata0.loc[0,'TotVol'] = sdata0['TotVol'][1]
                    
          
        # Jos yllÃ¤ olevien korjausten jÃ¤lkeekin aika on alussa > 1,
        # luodaan aika=0 -rivi.
    
        
        repair = 0
        repairB = 0
        
        if (sdata0['AIKA'][0]>1) & (sdata0['Age'][0]>=sdata0['AIKA'][0]): # Jos aika=0 tai aika=1, lisÃ¤tÃ¤Ã¤n alkuun rivi, joka kopioitu ensimmÃ¤iseltÃ¤ olemassa olevalta riviltÃ¤, kÃ¤ytÃ¤nnÃ¶ssÃ¤ usein varhaisperkauksen jÃ¤lkeinen tilanne
            # print('Aika alussa > 1, i=' + str(i))
            with open(outfile, "a") as myfile:
                myfile.write('\nAika alussa > 1, i=' + str(i) + ' kuvio: '+ str(kuviolista[i]))
             
            #sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
            nolla_df=pd.DataFrame([sdata0.iloc[0]])
            sdata0 = pd.concat((sdata0,nolla_df),axis=0, ignore_index=True) # Rivi menee taulukon loppuun
            sdata0['AIKA'][len(sdata0)-1] = 0 # tehdaan alkuun 0 vuosi
            
            sdata0 = sdata0.sort_values(by=['AIKA']) #sorattaan jotta 0-vuosi alkuun
            sdata0 = sdata0.reset_index(drop=True) #resetoidaan indeksi
            sdata0['Age'][0] = sdata0['Age'][0]-(sdata0['AIKA'][1] - sdata0['AIKA'][0]) #tehdaan n-vuodelle oikea ika
    
            
            # if (sdata0['Age'][0]<0) | (sdata0['Age'][1]<0):
            #     sdata0['Age'][0]=80 # Jos ikÃ¤Ã¤ ei saada puustot-tiedostosta, arvotaan 80 v (tÃ¤llÃ¤ ei pitÃ¤isi olla juurikaan merkitystÃ¤ mihinkÃ¤Ã¤n??)
    
                # # print('Second row age = 0 !!! i = ' +str(i))
                # with open(outfile, "a") as myfile:
                #     myfile.write('\nSecond row age = 0 !!! i = ' +str(i) + 'kuvio: '+str(kuviolista[i]))
            
            
            sdata0['TotVol'][0] = tulotdata['volume0'][0] #luodulle aika=0 riville tilavuus tulotdatan vuosittaisita ilavuusista
            
            # PÃ¤ivitetÃ¤Ã¤n muut komponentit
            for c in range(0,len(comps)):
                if sdata0.loc[0,comps[c]]>0: #jos arvo on yli 0
                    sdata0.loc[0,comps[c]] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, comps[c]] ####

            #/users/asalmiva/hiilipolku_leena_silvi_mottitiedostojenteko.py            
            # PÃ¤ivitetÃ¤Ã¤n N
            if sdata0['tuotosKuolleet'][0] < sdata0['TotVol'][0]:
                sdata0.loc[0,'N'] = (1 + sdata0['tuotosKuolleet'][0]/sdata0['TotVol'][0])*sdata0['N'][0]       ##miten tama on laskettu 
                     
            for c in range(0,len(biom_comp)):
                if sdata0.loc[0,biom_comp[c]]>0: #jos arvo on yli 0
                    sdata0.loc[0,biom_comp[c]] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, biom_comp[c]] ####runkohukka pitäisi olla ianakin ennen päätehhakkuta eripäin tuo suhde kuin muilla jten pitäisi skaalata eri tavoin... kasvaa 35 vuotaan saakka ehkä ja sitten laksee ainkain lehtipuilla? tai kuuso+lehti
    
            repair = 1
    
        # # Jos harvennus tulee jo ennen Motti-tulosteen toista riviÃ¤, luodaan AIKA=1 -rivi,
        # # jotta pystytÃ¤Ã¤n tehdÃ¤ Motti-tiedosto, josta saadaan interpolointifunktiot
        # # parille ensimmÃ¤iselle vuodelle
    
        # if len(vuodet)>0:
        #     if vuodet[0]>1: # LisÃ¤ys tehdÃ¤Ã¤n vain, jos hakkuu tulee vuonna 2 tai myÃ¶hemmin
        #         if sdata0['AIKA'][1]>=vuodet[0]: # Jos aika=0 tai aika=1, lisÃ¤tÃ¤Ã¤n alkuun rivi, joka kopioitu ensimmÃ¤iseltÃ¤ olemassa olevalta riviltÃ¤, kÃ¤ytÃ¤nnÃ¶ssÃ¤ varhaisperkauksen jÃ¤lkeinen tilanne
        
        #             with open(outfile, "a") as myfile:
        #                 myfile.write('\nAika[1]>vuodet[0]!!! ---- HUOM!!!, i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
                    
        #             sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
                    
        #             sdata0['AIKA'][len(sdata0)-1] = 1
                    
        #             sdata0 = sdata0.sort_values(by=['AIKA'])
        #             sdata0 = sdata0.reset_index(drop=True)
        #             sdata0['Age'][1] = sdata0['Age'][0]+1
        
                    
        #             # if (sdata0['Age'][1]<0) | (sdata0['Age'][1]<0):
        #             #     sdata0['Age'][1]=81 # Jos ikÃ¤Ã¤ ei saada puustot-tiedostosta, arvotaan 80 v (tÃ¤llÃ¤ ei pitÃ¤isi olla juurikaan merkitystÃ¤ mihinkÃ¤Ã¤n??)
        
        #                 # # print('Second row age = 0 !!! i = ' +str(i))
        #                 # with open(outfile, "a") as myfile:
        #                 #     myfile.write('\nSecond row age = 0 !!! i = ' +str(i))
                    
                    
        #             sdata0['TotVol'][1] = tulotdata['volume1'][0] 
            
        #             repairB = 1
    
    
        # Tsekataan, onko TotVol=0 -rivejÃ¤
        
        a = 0
        zeroRow = next((m for m, x in enumerate(sdata0['TotVol'][a:]) if x==0), None)
        
        if zeroRow is not None:
            
            zeroRow = a + zeroRow
                    
            while (zeroRow is not None) & (zeroRow<len(sdata0)-1): # Ei tarkisteta viimeistÃ¤ riviÃ¤ tÃ¤ssÃ¤, koska tuottaa ongelmia. Tsekataan alempana.
                
                new_vol = tulotdata['volume'+str(int(sdata0['AIKA'][zeroRow]))][0]
                
                if new_vol > 5: # PÃ¤Ã¤tehakkuun tapauksessa 0-rivi on Tulot-tiedostossa vasta kohdassa aika+1
                    new_vol = tulotdata['volume'+str(int(sdata0['AIKA'][zeroRow+1]))][0]
                
                sdata0.loc[zeroRow,'TotVol'] = new_vol
                
                # Haetaan seuraava 0-rivi:
                    
                a = zeroRow + 1    
                zeroRow = next((m for m, x in enumerate(sdata0['TotVol'][a:]) if x==0), None)
                
                if zeroRow is not None:
                    zeroRow = a + zeroRow
                else:
                    break
    
    
    
        # Jos 0-biomassoja, tarkistetaan myÃ¶s 0-iÃ¤t ja korjataan ikÃ¤ vastaamaan AIKA-jaksoja
    
        a = 0
        while a < len(sdata0):
            zeroRow = next((m for m, x in enumerate(sdata0['Age'][a:]) if x==0), None)
            if zeroRow is not None:
                zeroRow = zeroRow+a
                if sdata0['Neulasetkg'][zeroRow]==0: # Muokataan vain, jos neulasissa 0-rivejÃ¤
                    nonZeroRow = next((m for m, x in enumerate(sdata0['Age'][zeroRow:]) if x!=0), None)
                    if nonZeroRow is not None: # Jos 0-ikÃ¤ viimeisellÃ¤ rivillÃ¤, ei tehdÃ¤ tarkistuksia/korjauksia
                        nonZeroRow = nonZeroRow + zeroRow
                        if sdata0['AIKA'][nonZeroRow]-sdata0['AIKA'][zeroRow] != sdata0['Age'][nonZeroRow]-sdata0['Age'][zeroRow]:
    
                            sdata0.loc[zeroRow, 'Age'] = int(sdata0['Age'][nonZeroRow] - (sdata0['AIKA'][nonZeroRow]-sdata0['AIKA'][zeroRow]))
                            with open(outfile, "a") as myfile:
                                myfile.write('\n0-age corrected to: ' + str(sdata0.loc[zeroRow, 'Age']) + ', i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
    
                    else:
                        with open(outfile, "a") as myfile:
                            myfile.write('\n0-age detected at last row!! i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
    
                a = zeroRow + 1
            else:
                a = len(sdata0)+1     
                
                
    
    
        # Tsekataan TotVol ja DomH nollat ja muokataan
        # TÃ¤hÃ¤n tarvii ehkÃ¤ suuremman alarajan?
        # PitÃ¤Ã¤ lisÃ¤tÃ¤ myÃ¶s tyhjien tsekkaus!
        
        for k in range(0,len(sdata0)):
            
            if (sdata0.loc[k,'TotVol']<0.5) | (np.isnan(sdata0.loc[k,'TotVol']) == True): #voi olla 0.01 tilavuuksiakin
                
                if k < len(sdata0)-1:
                    sdata0.loc[k,'TotVol'] = max(0.01, 0.01*sdata0.loc[k+1,'TotVol']) # LisÃ¤tÃ¤Ã¤n totvol nollien tilalle muu pieni luku (0.001 oli liian pieni, kasvu ei lÃ¤htenyt kÃ¤yntiin, joten kasvatettu 0.2:een) paitsi et tulee 0.5
                else:
                    sdata0.loc[k,'TotVol'] = 0.01
                    
            if (sdata0.loc[k,'DomH']<0.3) | (np.isnan(sdata0.loc[k,'DomH']) == True):
                 
                 if k < len(sdata0)-1:
                     sdata0.loc[k,'DomH'] = max(0.3, 0.3*sdata0.loc[k+1,'DomH']) # LisÃ¤tÃ¤Ã¤n domh nollien tilalle muu pieni luku
                 else:
                     sdata0.loc[k,'DomH'] = 0.3       
       
        
        # Adding values to "missing" biomass components with similar increase-% per year as TotVol
        # Neulaset viimeisenÃ¤, jotta sitÃ¤ voidaan kÃ¤yttÃ¤Ã¤ arvioimaan, tarviiko biomassoja tÃ¤ydentÃ¤Ã¤.
        # Poistettu tÃ¤stÃ¤ listasta 'tuotosKuolleet', koska se on kumulatiivinen..
        
        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg', \
                     'TukPL1', 'TukPL2', 'TukPL3', 'TukPL4',\
                         'KuiPL1', 'KuiPL2', 'KuiPL3', 'KuiPL4',\
                         'HukPL1', 'HukPL2', 'HukPL3', 'HukPL4',\
                             'tuotosTilavuus',\
                                 'Runkopuukg', 'Hukkapuukg',\
                                     'MeanDbh','medianH', 'BA', 'tuotos']  
        #comps = ['DomH','MeanDbh','medianH',\
        #          'BA','Kannot_juuretkg',\
        #                  'tuotosTilavuus',\
        #                      'TukPL1', 'TukPL2','TukPL3',\
        #                          'KuiPL1', 'KuiPL2','KuiPL3',\
        #                              'HukPL1', 'HukPL2','HukPL3',\
        #                                  'VolPL1', 'VolPL2', 'VolPL3']  

        #biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg'] #tama lisatty tahan elokuu 23 jotta namakin tsekataan

        # if repairB==1: # Jos lisÃ¤tty 1-rivi  kopioimalla, korjataan sen komponentit ensin
        #     for j in range(0,len(biom_comp)):
        #         comp = biom_comp[j]
        #         sdata0.loc[1, comp] = sdata0['TotVol'][1]/sdata0['TotVol'][2]*sdata0.loc[2, comp]    
            
        # if repair==1: # Jos lisÃ¤tty 0-rivi kopioimalla, korjataan sen komponentit ensin
        #     for j in range(0,len(biom_comp)):
        #         comp = biom_comp[j]
        #         sdata0.loc[0, comp] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, comp]
                
                            
    
                
        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg']
    
        for j in range(0,len(biom_comp)):
        
            comp = biom_comp[j]
            k=0
    
            while k < len(sdata0)-1: 
                
                zeroRow = next((m for m, x in enumerate(sdata0[comp][k:]) if x==0), None)
    
                if zeroRow is not None: # If zeroRow exists
                    
                    zeroRow = zeroRow + k
                    
                    # if sdata0['Neulasetkg'][zeroRow]==0: # Muiden komponenttien numeroita muokatan vain, jos neulasissa on nollaa. TÃ¤mÃ¤ poistettu, koska nÃ¤issÃ¤ komponenteissa voi olla nollaa, vaikka neulasmassa ei olisikaan nollaa..
                    nonZeroRow = next((m for m, x in enumerate(sdata0[comp][zeroRow:]) if x!=0), None)
                    
                    if nonZeroRow is not None:
                        
                        nonZeroRow = nonZeroRow + zeroRow
                        sdata0.loc[zeroRow, comp] = sdata0['TotVol'][zeroRow]/sdata0['TotVol'][nonZeroRow]*sdata0.loc[nonZeroRow, comp] # tää vaikuttaa siihen myös et 0 riville tulee 0.5 vol ja seuraava on jotain
                    
                        k=k+1
                        
                    else:
                        
                        k=k+1                          
    
                else:
                    k=k+1
                            
    
                        
                        # else:
                            
                            # Jotta ei toistu jokaisella komponentilla erikseen, 
                            # nÃ¤mÃ¤ korjataan vasta lopuksi: 
                            # kopioidaan loppuun rivit, jotka vastaavat pÃ¤Ã¤hakkuun jÃ¤lkeistÃ¤ kehitystÃ¤
                            
                            # EtsitÃ¤Ã¤n toisen pÃ¤Ã¤hakkuun jÃ¤lkeinen rivi
            
    
             
                
             
    
    
    
        # TehdÃ¤Ã¤n sarake 'treesp'=pÃ¤Ã¤puulaji volumen perusteella
        
        sdata0['treesp'] = pd.Series([], dtype='int64') #        sdata0['treesp'] = pd.Series([])

        
        for t in range(0,len(sdata0)):
            
            
            volPine = sdata0['VolPL1'][t]
            volSpruce = sdata0['VolPL2'][t]
            volBirch = sdata0['VolPL3'][t] + sdata0['VolPL4'][t]
        
            spe = 3
        
            if volPine >= max(volSpruce, volBirch):
                spe = 1
            
            if volSpruce > max(volPine, volBirch):
                spe = 2    
    
            sdata0.loc[t, 'treesp']=spe
    
                            
                            
        # Tsekataan, lÃ¶ytyykÃ¶ viimeisen pÃ¤Ã¤tehakkuun jÃ¤lkeen pelkkiÃ¤ 0-rivejÃ¤. Jos nÃ¤in,
        # kopioidaan aiemmat pÃ¤Ã¤tehakkuun jÃ¤lkeiset rivit loppuun.
        # Korjataan sen jÃ¤lkeen myÃ¶s AIKA-sarake
    
        if len(paateh)>0:
            paate1 =  next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>paateh[0]), None) 
            nonZeroRow = next((m for m, x in enumerate(sdata0['Neulasetkg'][paate1:]) if x!=0), None)
            
            if len(paateh)>1: # Jos useampi kuin yksi pÃ¤Ã¤tehakkuu
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>=paateh[1]), None) # tÃ¤ssÃ¤ aiemmin x>paateh[1] !!
                nonZeroRow = next((m for m, x in enumerate(sdata0['Neulasetkg'][paaterow:]) if x!=0), None)
            
                if nonZeroRow is None:
                    
                    paate2aika = sdata0['AIKA'][paaterow]
                    
                    b = 0
                    
                    for p in range(paaterow,len(sdata0)):
                        print(p, paate2aika, "paatehakkuun jalkeiset nollarivit")
    
                        currow = paate1 + p - paaterow -1 # LisÃ¤tty tÃ¤hÃ¤n -1, koska rivillÃ¤ 617 lisÃ¤tty x>=paateh[1]
    
                        if currow<b: # LisÃ¤ys 17.10: Joissakin tapauksissa yllÃ¤ toimii -1, toisissa 0, joten muokataan currow tarvittaessa...
                            currow=b
                            
                        sdata0.loc[p]=sdata0.loc[currow]
                        
                        if b==0:
                            sdata0.loc[p, 'AIKA'] = int(paate2aika)
                            sdata0.loc[p, 'tuotosKuolleet'] = sdata0.loc[p-1, 'tuotosKuolleet']
                            b=b+1
                        else:
                            sdata0.loc[p, 'AIKA'] = int(paate2aika + sdata0.loc[currow, 'AIKA'] - sdata0.loc[currow-b, 'AIKA'])
                            sdata0.loc[p, 'tuotosKuolleet'] = sdata0.loc[p-1, 'tuotosKuolleet']
                            b=b+1
                
            else:
                # Jos vain yksi pÃ¤Ã¤tehakkuu, ei voida kopioida rivejÃ¤ aiemman pÃ¤Ã¤tehakkuun jÃ¤lkeen..
                print('MitÃ¤s tÃ¤hÃ¤n? kuvio:', str(kuviolista[i]),paateh)   
                if nonZeroRow is None:
                    with open(outfile, "a") as myfile:
                        myfile.write('\nPÃ¤Ã¤tehakkuun jÃ¤lkeisiÃ¤ 0-rivejÃ¤ ei saada korjattua! i=' + str(i)+ ' kuvio:'+ str(kuviolista[i])+ ' paateh: '+str(paateh)) #Ã¤Ã¤ 
    
        
        # Tsekataan, lÃ¶ytyykÃ¶ viimeiseltÃ¤ riviltÃ¤ (aika=100), nollariviÃ¤.
        # Jos lÃ¶ytyy, poistetaan.
        
    
        s100 = sdata0[(sdata0['AIKA']==100) & (sdata0['Age']==0)] # Vai vol=0?
        
        if len(s100)>0:
            sdata0 = sdata0[sdata0['AIKA']!=100]
            
        # -------------------------------------------------------------------
        # LisÃ¤ys 13.9.2021
        # Tarkistetaan, olisiko ensimmÃ¤inen Motti-tiedosto 1-rivin mittainen:
        
        vuodet_p = vuodet[vuodet>0]
        
        if (len(vuodet_p)>0):
            
            sdata = sdata0[sdata0['AIKA'] < vuodet_p[0]]
            
            if len(sdata)==1:
                #sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
                nolla2_df=pd.DataFrame([sdata0.iloc[0]])
                sdata0 = pd.concat((sdata0,nolla2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = vuodet_p[0]-1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                sdata0.loc[1,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][1]))][0]
                sdata0.loc[1,'Age'] = sdata0['Age'][0]+(sdata0['AIKA'][1]-sdata0['AIKA'][0])
    
                for c in range(0,len(comps)):
                    if sdata0.loc[0,comps[c]]>0:
                        # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
                        # sitÃ¤ voisi supistaa..??
                        sdata0.loc[1,comps[c]] = sdata0['TotVol'][1]/sdata0['TotVol'][0]*sdata0.loc[0, comps[c]]
    
                for c in range(0,len(biom_comp)):
                    if sdata0.loc[0,biom_comp[c]]>0:
                        # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
                        # sitÃ¤ voisi supistaa..??
                        sdata0.loc[1,biom_comp[c]] = sdata0['TotVol'][1]/sdata0['TotVol'][0]*sdata0.loc[0, biom_comp[c]]
    
    
        # -------------------------------------------------------------------
        # LisÃ¤ys 9.3.2022
        # Tarkistetaan, olisiko *toinen* Motti-tiedosto 1-rivin mittainen:
        
        #if (len(vuodet_p)>0):
        #    
        #    if len(vuodet_p)==1:
        #        end_year=101
        #    else:
        #        end_year = vuodet_p[1]
        #        
        #    sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[0]) & (sdata0['AIKA'] < end_year)] # Huom. tÃ¤ssÃ¤ ei saa resetoida indeksiÃ¤, jotta jatko toimii
        #    
        #    if len(sdata)==1:
        #        sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # LisÃ¤tÃ¤Ã¤n rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
        #        sdata0['AIKA'][len(sdata0)-1] = end_year - 1
        # 
        #         sdata0 = sdata0.sort_values(by=['AIKA'])
        #         sdata0 = sdata0.reset_index(drop=True)
        #         
        #        sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
        #        
        #        sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
        #        
        #        if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa riviÃ¤
        #            sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
        #
        #            for c in range(0,len(comps)):
        #                if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
        #                    # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
        #                    # sitÃ¤ voisi supistaa..??
        #                    sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
        # 
        #            for c in range(0,len(biom_comp)):
        #                if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
        #                    # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
        #                    # sitÃ¤ voisi supistaa..??
        #                    sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
        
        
        # -------------------------------------------------------------------
        # Update 6.3.2023: checking all 1-row-motti-files except the last one
        
        for v in range(1, len(vuodet_p)):
            
            if (len(vuodet_p)>0):
                
                if len(vuodet_p)==1:
                    end_year=101
                else:
                    end_year = vuodet_p[v]                
                
                sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[v-1]) & (sdata0['AIKA'] < end_year)] # Huom. tÃ¤ssÃ¤ ei saa resetoida indeksiÃ¤ jotta jatko toimii
                
                if len(sdata)==1:
                    #sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # LisÃ¤tÃ¤Ã¤n rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
                    sdataind_df=pd.DataFrame([sdata0.iloc[sdata.index[0]]])
                    sdata0 = pd.concat((sdata0,sdataind_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                    sdata0['AIKA'][len(sdata0)-1] = end_year - 1
                    
                    sdata0 = sdata0.sort_values(by=['AIKA'])
                    sdata0 = sdata0.reset_index(drop=True)
                    
                    sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
                    
                    sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
                    
                    if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa rivia
                        sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
                        
                        for c in range(0,len(comps)):
                            if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
                                # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                                # sita voisi supistaa..??
                                sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
                        
                        for c in range(0,len(biom_comp)):
                            if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
                                # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                                # sita voisi supistaa..??
                                sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
        
        # Update 6.3.2023: Checking if the last motti-file would be a 1-row-file
        
        if (len(vuodet_p)>1):
            
            sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[len(vuodet_p)-1])] # Huom. tassa ei saa resetoida indeksia, jotta jatko toimii
            end_year = 101
            
            if (len(sdata)==1) & (sdata0['AIKA'][len(sdata0)-1]<100):
                
                #sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # Lisataan rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
                sdataind2_df=pd.DataFrame([sdata0.iloc[sdata.index[0]]])
                sdata0 = pd.concat((sdata0, sdataind2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                    
                sdata0['AIKA'][len(sdata0)-1] = end_year - 1
                
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
                
                sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
                
                if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa rivia
                    sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
                    
                    for c in range(0,len(comps)):
                        if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
                            # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                            # sita voisi supistaa..??
                            sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
                    
                    for c in range(0,len(biom_comp)):
                        if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
                            # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                            # sita voisi supistaa..??
                            sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
                            
                
        
        # -------------------------------------------------------------------
        # LisÃ¤ys 9.3.2022
        # Tarkistetaan, tuleeko perÃ¤kkÃ¤isten harvennusten takia liian paljon
        # vÃ¤heneviÃ¤ tilavuuksia 
        
        comps_notk = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg', \
                             'tuotosTilavuus','MeanDbh','medianH', 'BA', 'tuotos', 'TotVol']  
        
        
        if pera >= 0:
            
            vuodet_p = vuodet[vuodet>=0] # LisÃ¤ys 17.10.2022: vuodet>=0 (jos tÃ¤ssÃ¤ olisi vuodet>0, vuodet_p saattaisi jÃ¤Ã¤dÃ¤ tyhjÃ¤ksi -> virhe seuraavalla rivillÃ¤)
            
            sdata = sdata0[sdata0['AIKA'] < vuodet_p[0]]
            
            if (len(sdata)<=2) & (len(sdata)>0): # LisÃ¤ys 17.10.2022: len(sdata)>0, koska jos sdata on tyhjÃ¤ -> error
                
                if (sdata['TotVol'][0] - sdata['TotVol'][1]) > 0:
                    
                    sdata0.loc[0]=sdata0.loc[1]
                    sdata0.loc[0, 'AIKA'] = sdata0['AIKA'][1] -1
                    sdata0.loc[0, 'Age'] = sdata0['Age'][1] -1
                    
                    
                    for c in range(0,len(comps_notk)):
                        
                        comp = comps_notk[c]
                        sdata0.loc[0, comp] = sdata0.loc[0, comp] * 0.97 # Arvio, jotta 
                                            
    
    
        
        
        # Muutetaan Age=-1 --> Age=0, nÃ¤itÃ¤ lÃ¶ytyi suoraan Motti-datasta
        
        sdata0['Age'] = np.clip(sdata0['Age'], 0, np.nan) 
        
        # -------------------------------------------------------------------
        # Lisays 28.3.2023:
        # Tarkistetaan, loytyyko iasta perakkaisia nollia
        
        a = 0
        while a < len(sdata0):
            zeroRow = next((m for m, x in enumerate(sdata0['Age'][a:]) if x==0), None)
            
            if zeroRow is not None:
                zeroRow = zeroRow+a
                try:
                    if sdata0['Age'][zeroRow+1]==0:
                        sdata0.loc[zeroRow, 'Age'] = 1
                        sdata0.loc[zeroRow+1, 'Age'] = 1 + int(sdata0['AIKA'][zeroRow+1]-sdata0['AIKA'][zeroRow])
                        with open(outfile, "a") as myfile:
                            myfile.write('\nDouble zero ages!! i=' + str(i))
                        
                    a = zeroRow + 1
                except:
                    a = zeroRow + 1
                        
            else:
                a = len(sdata0)+1     
        
        
        
        
        # Varsinaisten Motti-tiedostojen teko alkaa tÃ¤stÃ¤
        
        if len(sdata0)>0:
            
            if len(vuodet)>0:
                
                for j in range(0,n_mottifiles):
                
                    if j == 0:
                        sdata = sdata0[sdata0['AIKA'] < vuodet[j]]
                        if len(sdata)==0:
                            continue
                    
                    if (j>0) & (j<n_mottifiles-1):
                        sdata = sdata0[(sdata0['AIKA'] >= vuodet[j-1]) & (sdata0['AIKA'] < vuodet[j])]
                        sdata = sdata.reset_index(drop=True)
                        if len(sdata)==0:
                            continue
                    
                    if j == n_mottifiles-1:
                        sdata = sdata0[sdata0['AIKA'] >= vuodet[j-1]]
                        sdata = sdata.reset_index(drop=True)
                        if len(sdata)==0:
                            continue
                    
                    kasvatus = pd.Series(np.ones(len(sdata)))
                    vuosi = sdata['AIKA']
                    ika = sdata['Age']
                    N = sdata['N']
                    PPA = sdata['BA']
                    Hg = sdata['medianH']
                    Dg = sdata['MeanDbh']
                    Hdom = sdata['DomH']
                    Tilavuus = sdata['TotVol']
                    Tukki = sdata['TukPL1']+sdata['TukPL2']+sdata['TukPL3']+sdata['TukPL4']
                    Kuitu = sdata['KuiPL1']+sdata['KuiPL2']+sdata['KuiPL3']+sdata['KuiPL4']
                    Hukka = sdata['HukPL1']+sdata['HukPL2']+sdata['HukPL3']+sdata['HukPL4']
                    Tuotos = sdata['tuotosKuolleet'] + sdata['tuotosTilavuus']
                    Kuolleisuus = Tuotos*0. # Ei kÃ¤ytetÃ¤ Susissa
                    runko_aines = sdata['Runkopuukg']/1000. # YksikkÃ¶ tn
                    runko_hukka = sdata['Hukkapuukg']/1000.
                    elavat_oksat = sdata['Oksatkg']*0.8/1000. # Oksatkg = elÃ¤vÃ¤t ja kuolleet oksat yhteensÃ¤
                    kuolleet_oksat = sdata['Oksatkg']*0.2/1000.
                    lehdet = sdata['Neulasetkg']/1000.
                    kannot = sdata['kannot']/1000.
                    juuret = sdata['>2mm_juuretkg']/1000.
                    hienojuuret = sdata['h_juuretkg']/1000.
                    
                    # Tsekkaa, ettÃ¤ tÃ¤mÃ¤ on ok!!!
                    for t in range(0,len(sdata)):
                        if t==0:
                            Tuotos[0]=sdata['tuotosTilavuus'][0]
                        
                        else:
                            Tuotos[t] = sdata['tuotosTilavuus'][t] + (sdata['tuotosKuolleet'][t] - sdata['tuotosKuolleet'][t-1])
                    
                    # ---tÃ¤hÃ¤n asti
                    
                    data = pd.DataFrame(np.transpose(np.array([kasvatus, vuosi, ika, N, PPA, Hg, Dg, Hdom, Tilavuus, Tukki, Kuitu, \
                                         Hukka, Tuotos, Kuolleisuus, runko_aines, runko_hukka, \
                                         elavat_oksat, kuolleet_oksat, lehdet, kannot, juuret, hienojuuret])), \
                                         columns=['Kasvatus', 'Vuosi', 'Ikä', 'N', 'PPA', 'Hg', 'Dg', 'Hdom', 'Tilavuus', 'Tukki', 'Kuitu',\
                                                  'Hukka', 'Tuotos', 'Kuolleisuus', 'runko(aines)', 'runko(hukka)', \
                                                  'elävät oksat', 'kuolleet oksat', 'lehdet', 'Kannot', 'Juuret >2mm', 'Hienojuuret'])
                    
                    standInfo = kuviot[kuviot['KUVIO']==kuviolista[i]]
                    standInfo = standInfo.reset_index(drop=True)
                    
                    """ TSEKKAA, ETTÃ„ PUULAJI TULEE OIKEIN! """
                    
                    
                    
                    if (project == 'life') | (project == 'hiilipolku'):
                        spe = sdata['treesp'][len(sdata)-1]
                    
                    else: 
                        if (project =='vmi') | (project=='suo'):
                            
                            spe = int(sdata['treesp'][len(sdata)-1]) # viimeisen rivin puulaji
                            # treesp=kdata['Selite'][0][28:33]
                            # spe = 3
                            
                            # if treesp=='MÃ¤nty':
                            #     spe = 1
                            # if treesp=='Kuusi':
                            #     spe = 2       
                                
                            
                        else:
                            
                            treesp = pdata2['Selite'][0][28:33]
                            
                            spe = 3
                            
                            if treesp=='Mänty': #onk tÃ¤ a va Ã¤#
                                spe = 1
                            if treesp=='Kuusi':
                                spe = 2 
                    
                    
                    # print(kuviolista[i], ': puulaji', spe)
                    with open(outfile, "a") as myfile:
                        myfile.write('\n' + str(kuviolista[i]) +  ': puulaji ' + str(spe)+ ' n'+str(j))
                    
                    
                    fout = mottifolder + '/' + str(kuviolista[i]) + '_n' + str(j) + '.xls'
                    
                    cols = pd.DataFrame([['','','','',spe,'','','','','','','','','']],columns=['Kasvatus',	'Vuosi',	'id Harvennus',	'Harvennus',	'id Puulaji',	'Puulaji',	'Tukki[m³/ha]',	'Pikkutukki[m³/ha]',	'Kuitu[m³/ha]',	'Energiapuu, runko(aines)[m³/ha]',	'Energiapuu, runko(hukka)[m³/ha]',	'Energiapuu, oksat(e)[m³/ha]',	'Energiapuu, oksat(k)[m³/ha]',	'Energiapuu, kannot ja juuret[m³/ha]'])
                    
                    writer = pd.ExcelWriter(fout, engine = 'xlsxwriter')
                    data.to_excel(writer, index=False, sheet_name = 'Puustotunnukset')
                    cols.to_excel(writer, index=False, sheet_name = r'Kertymät')
                    #writer.save()
                    writer.close()
            else:
                # print(str(kuviolista[i]) + ': poistumat puuttuvat! i=' + str(i))
                with open(outfile, "a") as myfile:
                    myfile.write('\n' + str(kuviolista[i]) + ': poistumat puuttuvat! i=' + str(i))
                
                sdata = sdata0.copy()
                
                kasvatus = pd.Series(np.ones(len(sdata)))
                vuosi = sdata['AIKA']
                ika = sdata['Age']
                N = sdata['N']
                PPA = sdata['BA']
                Hg = sdata['medianH']
                Dg = sdata['MeanDbh']
                Hdom = sdata['DomH']
                Tilavuus = sdata['TotVol']
                Tukki = sdata['TukPL1']+sdata['TukPL2']+sdata['TukPL3']+sdata['TukPL4']
                Kuitu = sdata['KuiPL1']+sdata['KuiPL2']+sdata['KuiPL3']+sdata['KuiPL4']
                Hukka = sdata['HukPL1']+sdata['HukPL2']+sdata['HukPL3']+sdata['HukPL4']
                Tuotos = sdata['tuotos']
                Kuolleisuus = Tuotos*0. # Ei kÃ¤ytetÃ¤ Susissa
                runko_aines = sdata['Runkopuukg']/1000. # YksikkÃ¶ tn
                runko_hukka = sdata['Hukkapuukg']/1000.
                elavat_oksat = sdata['Oksatkg']*0.8/1000. # Oksatkg = elÃ¤vÃ¤t ja kuolleet oksat yhteensÃ¤
                kuolleet_oksat = sdata['Oksatkg']*0.2/1000.
                lehdet = sdata['Neulasetkg']/1000.
                kannot = sdata['kannot']/1000.
                juuret = sdata['>2mm_juuretkg']/1000.
                hienojuuret = sdata['h_juuretkg']/1000.
                
                data = pd.DataFrame(np.transpose(np.array([kasvatus, vuosi, ika, N, PPA, Hg, Dg, Hdom, Tilavuus, Tukki, Kuitu, \
                                     Hukka, Tuotos, Kuolleisuus, runko_aines, runko_hukka, \
                                     elavat_oksat, kuolleet_oksat, lehdet, kannot, juuret, hienojuuret])), \
                                     columns=['Kasvatus', 'Vuosi', 'Ikä', 'N', 'PPA', 'Hg', 'Dg', 'Hdom', 'Tilavuus', 'Tukki', 'Kuitu',\
                                              'Hukka', 'Tuotos', 'Kuolleisuus', 'runko(aines)', 'runko(hukka)', \
                                              'elävät oksat', 'kuolleet oksat', 'lehdet', 'Kannot', 'Juuret >2mm', 'Hienojuuret'])
                
                standInfo = kuviot[kuviot['KUVIO']==kuviolista[i]]
                standInfo = standInfo.reset_index(drop=True)
                
                
                """ TSEKKAA, ETTÃ„ PUULAJI TULEE OIKEIN! """
                if (project == 'life') | (project == 'hiilipolku'):
                    spe = sdata['treesp'][len(sdata)-1]
                
                else:
                    
                    if project =='vmi':
                        
                        # spe = int(sdata['treesp'][len(sdata)-1]) # viimeisen rivin puulaji
                       
                        treesp=kdata['Selite'][0][28:33]
                        spe = 3
                        
                        if treesp=='Mänty': #Ã¤
                            spe = 1
                        if treesp=='Kuusi':
                            spe = 2    
                    else:
                        treesp = pdata2['Selite'][0][28:33]
                        
                        spe = 3
                        
                        if treesp=='Mänty':
                            spe = 1
                        if treesp=='Kuusi':
                            spe = 2 
                    
                # print(kuviolista[i], ': puulaji ', spe)
                with open(outfile, "a") as myfile:
                    myfile.write('\n' + str(kuviolista[i]) +  ':' + str(spe)+ ' n0')
                
                
                
                fout = mottifolder + '/' + str(kuviolista[i]) + '_n0.xls'
                
                cols = pd.DataFrame([['','','','',spe,'','','','','','','','','']],columns=['Kasvatus',	'Vuosi',	'id Harvennus',	'Harvennus',	'id Puulaji',	'Puulaji',	'Tukki[m³/ha]',	'Pikkutukki[m³/ha]',	'Kuitu[m³/ha]',	'Energiapuu, runko(aines)[m³/ha]',	'Energiapuu, runko(hukka)[m³/ha]',	'Energiapuu, oksat(e)[m³/ha]',	'Energiapuu, oksat(k)[m³/ha]',	'Energiapuu, kannot ja juuret[m³/ha]'])
                
                writer = pd.ExcelWriter(fout, engine = 'xlsxwriter')
                data.to_excel(writer, index=False, sheet_name = 'Puustotunnukset')
                cols.to_excel(writer, index=False, sheet_name = r'Kertymät')
                #writer.save()
                writer.close()


def create_motti_files_silvi_lp4_volnage(project, scen, maku=None):       #alkutilavuus, annetaan olla 0.01tai 0.01*Totvol[1]lp4 skaalaukseen otettu lp4 mukaan  # Huom! ohje>1
    '''
    Modified: January 2023
        
    Create Motti input files for Susi simulations from "Skene-Motti" results.
    
    Works technically, but contains many small issues that should be fixed.
    
    project = 'life' (Hydrology Life -project) or 'suo' (SUO-project)
    project-parameter affects only main tree species (and column names)!
    
    
    Modified 8.3.2022:
        - korjattu pÃ¤Ã¤tehakkuuseen liittyviÃ¤ ongelmia; hdom pitÃ¤isi nyt saada
        jÃ¤rkeviÃ¤ arvoja, biomassojen ei pitÃ¤isi jÃ¤Ã¤dÃ¤ pÃ¤Ã¤tehakkuun jÃ¤lkeen nollaksi
        (ainakaan toisen pÃ¤Ã¤tehakkuun jÃ¤lkeen)
        
    Issues to fix:
        - Jos pÃ¤Ã¤tehakkuu tulee vuonna 2, pitÃ¤Ã¤ vuosille 0 ja 1 hakea muut komponentit
        totvol perusteella: hdom saa tÃ¤ssÃ¤ liian suuria arvoja. Tuotos jÃ¤Ã¤ nollaksi?!
        - See comments below...
        
    
    '''
    
    inputmainfol=r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan_volnage/'
    mottifolder = inputmainfol+scen[:-6]+'/'+scen+'/'
    if not os.path.exists(inputmainfol):
        os.mkdir(inputmainfol)
    if not os.path.exists(inputmainfol+scen[:-6]+'/'):
        os.mkdir(inputmainfol+scen[:-6]+'/')
    if not os.path.exists(mottifolder):
        os.mkdir(mottifolder)
            
    outfile = mottifolder + '/out.txt'
    
    project = 'hiilipolku'
    
    kuviot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kuviot.csv', encoding='latin1', sep=';')
    puustot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_puustot.csv', encoding='latin1', sep=';')
    poistumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_poistumat.csv', encoding='latin1', sep=';')
    tapahtumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_tapahtumat.csv', encoding='latin1', sep=';', index_col=False)
    tulot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/'+scen[:-6]+'/'+scen+'/motti/'+scen+'_kasvut.csv', encoding='latin1', sep=';')
    
    
    if (project == 'life') | (project =='vmi') | (project =='hiilipolku'):
        kuviot = kuviot.rename(columns={'kuvio':'KUVIO'})
        puustot = puustot.rename(columns={'kuvio':'KUVIO',
                                          '_2mm_juuretkg':'>2mm_juuretkg',})
        poistumat = poistumat.rename(columns={'kuvio':'KUVIO'})
        tapahtumat = tapahtumat.rename(columns={'kuvio':'KUVIO',
                                      '_2003': '2003',
                                      '_2004': '2004',
                                      '_2005': '2005'})
        tulot = tulot.rename(columns={'kuvio':'KUVIO'})
    
    
    kuviolista = list(set(kuviot['KUVIO'][kuviot['kuivatustilanne']>=3])) #tassa valitaan 3 Ojitettu kangas ja muut?DRAINAGESTATE	, ojitettu kangas jos pÃ¤Ã¤ryhmÃ¤1 ja alaryhmÃ¤ 3 ja kuivatustilanne 5?? 1 Ojittamaton kangas, 2 Soistunut kangas,3 Ojitettu kangas,	6 luonnontilainen suo,	7 Ojikko, 8	Muuttuma,9	Turvekangas
    
    
    for i in range(0,len(kuviolista)):
    
        print(i, 'stand: ', kuviolista[i])
        
        # alta poistettu ohjeeseen ja skenaarioon viittaavat rajaukset,
        # oletetaan, ettÃ¤ hiilipolussa vain yksi ohje per kuvio
        
        sdata0 = puustot[(puustot['KUVIO']==kuviolista[i])]
        pdata = poistumat[(poistumat['KUVIO']==kuviolista[i]) & (poistumat['HARV']>1)] #mita harv koodit meinaa?
        pdata2 = poistumat[(poistumat['KUVIO']==kuviolista[i])]
        tdata = tapahtumat[(tapahtumat['KUVIO']==kuviolista[i])]
        tulotdata = tulot[(tulot['KUVIO']==kuviolista[i])]
        kdata = kuviot[kuviot['KUVIO']==kuviolista[i]]
        
        sdata0 = sdata0.reset_index(drop=True)
        pdata = pdata.reset_index(drop=True)
        pdata2 = pdata2.reset_index(drop=True)
        tdata = tdata.reset_index(drop=True)
        tulotdata = tulotdata.reset_index(drop=True)
        kdata = kdata.reset_index(drop=True)
        
        tdata = tdata[(tdata['2003']>-1) | (tdata['2004']>-1) | (tdata['2005']>-1)]
        
        ensih = tdata[tdata['2003']>-1]['2003'].values
        harv = tdata[tdata['2004']>-1]['2004'].values
        paateh = tdata[tdata['2005']>-1]['2005'].values
        #voisi kayda tsekkaamassa onko taimikonhoito tai varhaisperkaus niin sitten ei muuttaisi puulajia tai ottaisi sen muualta et sailyy kuusena tai mantyna? 
        #tdic={'_2001':'lannoitus','_2002':'kunnostusojitus','_2003':'ensiharvennus','_2004':'harvennus','_2005':'paatehakkuu','_2006':'varhaiperkaus','_2007':'taimikonhoito','_2008':'viljely','_2009':'luontainen uudistuminen'}
        # taimh = tdata[tdata['_2007']>-1]['_2007'].values
        # Pudotetaan viimeisen vuoden tapahtumat pois, koska sekoittavat
        # myÃ¶hemmin tulevan koodin
        ensih = ensih[ensih<100]
        harv = harv[harv<100]
        paateh = paateh[paateh<100]
        
        vuodet = np.unique(np.sort(np.concatenate((ensih, harv, paateh))))
        
        ap = 0
        
        pera = -1
        
        # Tsekataan, onko tapahtumissa perÃ¤kkÃ¤isiÃ¤ vuosia ja poistetaan niistÃ¤ jÃ¤lkimmÃ¤inen
        for v in range(0,len(vuodet)-1):
            if vuodet[v+1]-vuodet[v]==1:
                
                ensih = np.delete(ensih, np.where(ensih == vuodet[v+1]))
                harv = np.delete(harv, np.where(harv == vuodet[v+1]))
                paateh = np.delete(paateh, np.where(paateh == vuodet[v+1]))
                
                vuodet = np.delete(vuodet, [v+1], None)
                
                pera = v
                
                ap = ap + 1
                
            if v+1>=len(vuodet)-1:
                break
        
        if ap>1:
            with open(outfile, "a") as myfile:
                myfile.write('\n\nKaksi kertaa perakkaiset vuodet! stand=' + str(kuviolista[i]) + ', i=' + str(i))
    
        n_mottifiles = len(vuodet) + 1
    
    
        if len(sdata0)==0: # skipataan kuvio, jos dataa ei ole
            continue
    
        if len(sdata0)>21:
            print('sdata0: liikaa rivejÃ¤?!')
        
    
                                    
        # Luodaan pÃ¤Ã¤tehakkuun jÃ¤lkeiset rivit ja korjataan rivin kopioimisesta aiheutuneet
        # luvut
        
        # muokattu 6.8.2021: poistettu tuotos, tuotosKuolleet, tuotosTilavuus comps:sta
        comps = ['DomH','MeanDbh','medianH',\
                          'BA','Kannot_juuretkg',\
                              'tuotosTilavuus',\
                              'TukPL1', 'TukPL2','TukPL3','TukPL4',\
                                  'KuiPL1', 'KuiPL2','KuiPL3','KuiPL4',\
                                      'HukPL1', 'HukPL2','HukPL3','HukPL4',\
                                          'VolPL1', 'VolPL2', 'VolPL3','VolPL4',]  #lisatty 091023 PL4 ositteet mukaan

        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg'] #tama lisatty tahan elokuu 23 jotta namakin tsekataan

        for h in range(0,len(paateh)):
            
            sdata_paateh= sdata0[sdata0['AIKA']==paateh[h]].reset_index(drop=True)
            
            go = False
            
            if len(sdata_paateh)>0:
                if sdata_paateh.loc[0, 'Age']>5:
                    go = True
        
            if (len(sdata_paateh)==0) | (go==True): # Luodaan rivit vain, jos ensimmÃ¤istÃ¤ riviÃ¤ ei ole mukana
                
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>paateh[h]), None)
                
                #sdata0 = sdata0.append([sdata0.iloc[paaterow]], ignore_index=True) # Rivi menee taulukon loppuun
                paaterow_df=pd.DataFrame([sdata0.iloc[paaterow]])
                sdata0 = pd.concat((sdata0,paaterow_df),axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = paateh[h]+1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==paateh[h]+1), None)
                
                # sdata0.loc[paaterow,'Age'] = 1 # TÃ¤mÃ¤ korvattu alla olevalla
                
                if paaterow>0:
                    if sdata0.loc[paaterow-1, 'Age']<5:
                        sdata0.loc[paaterow,'Age'] = sdata0.loc[paaterow-1, 'Age']+1
                    else:
                        sdata0.loc[paaterow,'Age'] = sdata0.loc[paaterow+1, 'Age']-1
                
                sdata0.loc[paaterow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][paaterow]))][0]
                
                
                if sdata0['tuotosKuolleet'][paaterow] < sdata0['TotVol'][paaterow]: #vaihdettiin taman paikkaa
                    sdata0.loc[paaterow,'N'] = (1 + sdata0['tuotosKuolleet'][paaterow]/sdata0['TotVol'][paaterow])*sdata0['N'][paaterow]        

                    if paaterow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][paaterow] > sdata0['N'][paaterow-1]:
                            if paaterow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[paaterow,'N'] = (sdata0['N'][paaterow-1] + sdata0['N'][paaterow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[paaterow,'N'] = 0.99*sdata0['N'][paaterow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.


                for c in range(0,len(comps)):
                    if sdata0.loc[paaterow,comps[c]]>0: #muutos
                        sdata0.loc[paaterow,comps[c]] = (sdata0['TotVol'][paaterow]/sdata0['N'][paaterow]*sdata0['Age'][paaterow])/(sdata0['TotVol'][paaterow+1]/sdata0['N'][paaterow+1]*sdata0['Age'][paaterow+1])*sdata0.loc[paaterow+1, comps[c]]
                        print(kuviolista[i],comps[c], "rivi 303 ongelma")
                    
                            
                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[paaterow,biom_comp[c]]>0: 
                        sdata0.loc[paaterow,biom_comp[c]] = (sdata0['TotVol'][paaterow]/sdata0['N'][paaterow]*sdata0['Age'][paaterow])/(sdata0['TotVol'][paaterow+1]/sdata0['N'][paaterow+1]*sdata0['Age'][paaterow+1])*sdata0.loc[paaterow+1, biom_comp[c]] ####
                
            
        # Luodaan ensiharvennuksen jÃ¤lkeiset rivit (jos vuosi ei osu kohdallaan muuten)
        # EsimerkissÃ¤ nÃ¤ytti siltÃ¤, ettÃ¤ ensiharvennusvuoden riville tulostuu harventamaton tilavuus! Siksi tÃ¤ssÃ¤ kÃ¤sitelty riviÃ¤+1
        # Toisessa esimerkissÃ¤ tulostui harvennettu tilavuus! Joten vaihtelee.
        
        for h in range(0,len(ensih)):
    
    
            ensihrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==ensih[h]+1), None)
            ensihrow2 = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>ensih[h]), None)
            
            if ensihrow is None:
                
                #sdata0 = sdata0.append([sdata0.iloc[ensihrow2]], ignore_index=True) # Rivi menee taulukon loppuun
                ensihrow2_df=pd.DataFrame([sdata0.iloc[ensihrow2]])
                sdata0 = pd.concat((sdata0,ensihrow2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = ensih[h]+1
                
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                ensihrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==ensih[h]+1), None)
                
                # Jos ikÃ¤ lÃ¶ytyy riviÃ¤ alaempaa:
                if sdata0.loc[ensihrow+1,'Age']>0:
                    sdata0.loc[ensihrow,'Age'] = sdata0.loc[ensihrow+1,'Age']-(sdata0['AIKA'][ensihrow+1]-sdata0['AIKA'][ensihrow])
                # Muussa tapauksessa ikÃ¤ lÃ¶ytyy riviÃ¤ ylempÃ¤Ã¤:
                else:
                    if ensihrow>0:
                        sdata0.loc[ensihrow,'Age'] = sdata0.loc[ensihrow-1,'Age']+(sdata0['AIKA'][ensihrow]-sdata0['AIKA'][ensihrow-1])
                    else:
                        # Jos ei lÃ¶ydy kummaltakaan puolelta:
                        with open(outfile, "a") as myfile:
                            myfile.write('\nEnsiharvennuksen jÃ¤lkeistÃ¤ ikÃ¤Ã¤ ei lÃ¶ydy! kuvio=' + str(kuviolista[i]) + ' i:'+ str(i))
                        
                    
                sdata0.loc[ensihrow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][ensihrow]))][0]
                
                # PÃ¤ivitetÃ¤Ã¤n N
                if sdata0['tuotosKuolleet'][ensihrow] < sdata0['TotVol'][ensihrow]:
                    sdata0.loc[ensihrow,'N'] = (1 + sdata0['tuotosKuolleet'][ensihrow]/sdata0['TotVol'][ensihrow])*sdata0['N'][ensihrow]        
    
                    if ensihrow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][ensihrow] > sdata0['N'][ensihrow-1]:
                            if ensihrow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[ensihrow,'N'] = (sdata0['N'][ensihrow-1] + sdata0['N'][ensihrow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[ensihrow,'N'] = 0.99*sdata0['N'][ensihrow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.

                
                for c in range(0,len(comps)):
                    if sdata0.loc[ensihrow,comps[c]]>0:
                        sdata0.loc[ensihrow,comps[c]] = (sdata0['TotVol'][ensihrow]/sdata0['N'][ensihrow]*sdata0['Age'][ensihrow])/(sdata0['TotVol'][ensihrow+1]/sdata0['N'][ensihrow+1]*sdata0['Age'][ensihrow+1])*sdata0.loc[ensihrow+1, comps[c]]
                
                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[ensihrow,biom_comp[c]]>0: 
                        sdata0.loc[ensihrow,biom_comp[c]] = (sdata0['TotVol'][ensihrow]/sdata0['N'][ensihrow]*sdata0['Age'][ensihrow])/(sdata0['TotVol'][ensihrow+1]/sdata0['N'][ensihrow+1]*sdata0['Age'][ensihrow+1])*sdata0.loc[ensihrow+1, biom_comp[c]] ####
                            
        
        # Luodaan harvennuksen jÃ¤lkeiset rivit (jos vuosi ei osu kohdalleen muuten)
        # Harvennusvuoden rivillÃ¤ on harventamaton tilavuus. Siksi kÃ¤sitelty riviÃ¤+1.
        
        for h in range(0,len(harv)):
    
    
            hrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==harv[h]+1), None)
            hrow2 = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>harv[h]), None)
            
            if hrow is None:
                
                #sdata0 = sdata0.append([sdata0.iloc[hrow2]], ignore_index=True) # Rivi menee taulukon loppuun
                hrow2_df=pd.DataFrame([sdata0.iloc[hrow2]])
                sdata0 = pd.concat((sdata0,hrow2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = harv[h]+1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                hrow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x==harv[h]+1), None)
                
                # Jos ikÃ¤ lÃ¶ytyy riviÃ¤ alaempaa:
                if sdata0.loc[hrow+1,'Age']>0:
                    sdata0.loc[hrow,'Age'] = sdata0.loc[hrow+1,'Age']-(sdata0['AIKA'][hrow+1]-sdata0['AIKA'][hrow])
                # Muussa tapauksessa ikÃ¤ lÃ¶ytyy riviÃ¤ ylempÃ¤Ã¤:
                else:
                    if hrow>0:
                        sdata0.loc[hrow,'Age'] = sdata0.loc[hrow-1,'Age']+(sdata0['AIKA'][hrow]-sdata0['AIKA'][hrow-1])
                    else:
                    # Jos ei lÃ¶ydy kummaltakaan puolelta:
                        with open(outfile, "a") as myfile:
                            myfile.write('\nHarvennuksen jÃ¤lkeistÃ¤ ikÃ¤Ã¤ ei lÃ¶ydy! i=' + str(i) + ' kuvio:'+ str(kuviolista[i]))
                
                    
                sdata0.loc[hrow,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][hrow]))][0]
                
                # Korjataan N:
                if sdata0['tuotosKuolleet'][hrow] < sdata0['TotVol'][hrow]:
                    sdata0.loc[hrow,'N'] = (1 + sdata0['tuotosKuolleet'][hrow]/sdata0['TotVol'][hrow])*sdata0['N'][hrow]
                    
                    if hrow>0: # Tarkistetaan, ettei korjattu N ole suurempi kuin edellisen rivin N
                        if sdata0['N'][hrow] > sdata0['N'][hrow-1]:
                            if hrow < len(sdata0)-1: # Jos hrow ei ole viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[hrow,'N'] = (sdata0['N'][hrow-1] + sdata0['N'][hrow-1])/2. # N = edellisen ja seuraavan keskiarvo
                            else: # Jos hrow on viimeisellÃ¤ rivillÃ¤...
                                sdata0.loc[hrow,'N'] = 0.99*sdata0['N'][hrow-1] # Arvioidaan 1% vÃ¤hennys runkolukuun. TÃ¤mÃ¤ ei ole hyvÃ¤ ratkaisu, mutta parempi (lÃ¤hempÃ¤nÃ¤ todellisuutta) kuitenkin kuin ilman mitÃ¤Ã¤n vÃ¤hennystÃ¤.
                            
                
                for c in range(0,len(comps)):
                    if sdata0.loc[hrow,comps[c]]>0:
                        sdata0.loc[hrow,comps[c]] = (sdata0['TotVol'][hrow]/sdata0['N'][hrow]*sdata0['Age'][hrow])/(sdata0['TotVol'][hrow+1]/sdata0['N'][hrow+1]*sdata0['Age'][hrow+1])*sdata0.loc[hrow+1, comps[c]]

                        
                for c in range(0,len(biom_comp)): #taa lisatty 09102023
                    if sdata0.loc[hrow,biom_comp[c]]>0: 
                        sdata0.loc[hrow,biom_comp[c]] = (sdata0['TotVol'][hrow]/sdata0['N'][hrow]*sdata0['Age'][hrow])/(sdata0['TotVol'][hrow+1]/sdata0['N'][hrow+1]*sdata0['Age'][hrow+1])*sdata0.loc[hrow+1, biom_comp[c]] ####
    
        # Tsekataan, ettei ole ekan rivin TotVol pyÃ¶ristyksen kanssa ongelmaa:
        
        if (sdata0['TotVol'][0]>10.) & (-0.1 < sdata0['TotVol'][1]-sdata0['TotVol'][0]<0):
            sdata0.loc[0,'TotVol'] = sdata0['TotVol'][1]
                    
          
        # Jos yllÃ¤ olevien korjausten jÃ¤lkeekin aika on alussa > 1,
        # luodaan aika=0 -rivi.
    
        
        repair = 0
        repairB = 0
        
        if (sdata0['AIKA'][0]>1) & (sdata0['Age'][0]>=sdata0['AIKA'][0]): # Jos aika=0 tai aika=1, lisÃ¤tÃ¤Ã¤n alkuun rivi, joka kopioitu ensimmÃ¤iseltÃ¤ olemassa olevalta riviltÃ¤, kÃ¤ytÃ¤nnÃ¶ssÃ¤ usein varhaisperkauksen jÃ¤lkeinen tilanne
            # print('Aika alussa > 1, i=' + str(i))
            with open(outfile, "a") as myfile:
                myfile.write('\nAika alussa > 1, i=' + str(i) + ' kuvio: '+ str(kuviolista[i]))
             
            #sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
            nolla_df=pd.DataFrame([sdata0.iloc[0]])
            sdata0 = pd.concat((sdata0,nolla_df),axis=0, ignore_index=True) # Rivi menee taulukon loppuun
            sdata0['AIKA'][len(sdata0)-1] = 0 # tehdaan alkuun 0 vuosi
            
            sdata0 = sdata0.sort_values(by=['AIKA']) #sorattaan jotta 0-vuosi alkuun
            sdata0 = sdata0.reset_index(drop=True) #resetoidaan indeksi
            sdata0['Age'][0] = sdata0['Age'][0]-(sdata0['AIKA'][1] - sdata0['AIKA'][0]) #tehdaan n-vuodelle oikea ika
    
            
            # if (sdata0['Age'][0]<0) | (sdata0['Age'][1]<0):
            #     sdata0['Age'][0]=80 # Jos ikÃ¤Ã¤ ei saada puustot-tiedostosta, arvotaan 80 v (tÃ¤llÃ¤ ei pitÃ¤isi olla juurikaan merkitystÃ¤ mihinkÃ¤Ã¤n??)
    
                # # print('Second row age = 0 !!! i = ' +str(i))
                # with open(outfile, "a") as myfile:
                #     myfile.write('\nSecond row age = 0 !!! i = ' +str(i) + 'kuvio: '+str(kuviolista[i]))
            
            
            sdata0['TotVol'][0] = tulotdata['volume0'][0] #luodulle aika=0 riville tilavuus tulotdatan vuosittaisita ilavuusista
            
            # PÃ¤ivitetÃ¤Ã¤n N
            if sdata0['tuotosKuolleet'][0] < sdata0['TotVol'][0]:
                sdata0.loc[0,'N'] = (1 + sdata0['tuotosKuolleet'][0]/sdata0['TotVol'][0])*sdata0['N'][0]       ##miten tama on laskettu 

            # PÃ¤ivitetÃ¤Ã¤n muut komponentit
            for c in range(0,len(comps)):
                if sdata0.loc[0,comps[c]]>0: #jos arvo on yli 0
                    sdata0.loc[0,comps[c]] = (sdata0['TotVol'][0]/sdata0['N'][0]*sdata0['Age'][0])/(sdata0['TotVol'][1]/sdata0['N'][1]*sdata0['Age'][1])*sdata0.loc[1, comps[c]] ####

            #/users/asalmiva/hiilipolku_leena_silvi_mottitiedostojenteko.py                                 
            for c in range(0,len(biom_comp)):
                if sdata0.loc[0,biom_comp[c]]>0: #jos arvo on yli 0
                    sdata0.loc[0,biom_comp[c]] = (sdata0['TotVol'][0]/sdata0['N'][0]*sdata0['Age'][0])/(sdata0['TotVol'][1]/sdata0['N'][1]*sdata0['Age'][1])*sdata0.loc[1, biom_comp[c]] ####runkohukka pitäisi olla ianakin ennen päätehhakkuta eripäin tuo suhde kuin muilla jten pitäisi skaalata eri tavoin... kasvaa 35 vuotaan saakka ehkä ja sitten laksee ainkain lehtipuilla? tai kuuso+lehti
    
            repair = 1
    
        # # Jos harvennus tulee jo ennen Motti-tulosteen toista riviÃ¤, luodaan AIKA=1 -rivi,
        # # jotta pystytÃ¤Ã¤n tehdÃ¤ Motti-tiedosto, josta saadaan interpolointifunktiot
        # # parille ensimmÃ¤iselle vuodelle
    
        # if len(vuodet)>0:
        #     if vuodet[0]>1: # LisÃ¤ys tehdÃ¤Ã¤n vain, jos hakkuu tulee vuonna 2 tai myÃ¶hemmin
        #         if sdata0['AIKA'][1]>=vuodet[0]: # Jos aika=0 tai aika=1, lisÃ¤tÃ¤Ã¤n alkuun rivi, joka kopioitu ensimmÃ¤iseltÃ¤ olemassa olevalta riviltÃ¤, kÃ¤ytÃ¤nnÃ¶ssÃ¤ varhaisperkauksen jÃ¤lkeinen tilanne
        
        #             with open(outfile, "a") as myfile:
        #                 myfile.write('\nAika[1]>vuodet[0]!!! ---- HUOM!!!, i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
                    
        #             sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
                    
        #             sdata0['AIKA'][len(sdata0)-1] = 1
                    
        #             sdata0 = sdata0.sort_values(by=['AIKA'])
        #             sdata0 = sdata0.reset_index(drop=True)
        #             sdata0['Age'][1] = sdata0['Age'][0]+1
        
                    
        #             # if (sdata0['Age'][1]<0) | (sdata0['Age'][1]<0):
        #             #     sdata0['Age'][1]=81 # Jos ikÃ¤Ã¤ ei saada puustot-tiedostosta, arvotaan 80 v (tÃ¤llÃ¤ ei pitÃ¤isi olla juurikaan merkitystÃ¤ mihinkÃ¤Ã¤n??)
        
        #                 # # print('Second row age = 0 !!! i = ' +str(i))
        #                 # with open(outfile, "a") as myfile:
        #                 #     myfile.write('\nSecond row age = 0 !!! i = ' +str(i))
                    
                    
        #             sdata0['TotVol'][1] = tulotdata['volume1'][0] 
            
        #             repairB = 1
    
    
        # Tsekataan, onko TotVol=0 -rivejÃ¤
        
        a = 0
        zeroRow = next((m for m, x in enumerate(sdata0['TotVol'][a:]) if x==0), None)
        
        if zeroRow is not None:
            
            zeroRow = a + zeroRow
                    
            while (zeroRow is not None) & (zeroRow<len(sdata0)-1): # Ei tarkisteta viimeistÃ¤ riviÃ¤ tÃ¤ssÃ¤, koska tuottaa ongelmia. Tsekataan alempana.
                
                new_vol = tulotdata['volume'+str(int(sdata0['AIKA'][zeroRow]))][0]
                
                if new_vol > 5: # PÃ¤Ã¤tehakkuun tapauksessa 0-rivi on Tulot-tiedostossa vasta kohdassa aika+1
                    new_vol = tulotdata['volume'+str(int(sdata0['AIKA'][zeroRow+1]))][0]
                
                sdata0.loc[zeroRow,'TotVol'] = new_vol
                
                # Haetaan seuraava 0-rivi:
                    
                a = zeroRow + 1    
                zeroRow = next((m for m, x in enumerate(sdata0['TotVol'][a:]) if x==0), None)
                
                if zeroRow is not None:
                    zeroRow = a + zeroRow
                else:
                    break
    
    
    
        # Jos 0-biomassoja, tarkistetaan myÃ¶s 0-iÃ¤t ja korjataan ikÃ¤ vastaamaan AIKA-jaksoja
    
        a = 0
        while a < len(sdata0):
            zeroRow = next((m for m, x in enumerate(sdata0['Age'][a:]) if x==0), None)
            if zeroRow is not None:
                zeroRow = zeroRow+a
                if sdata0['Neulasetkg'][zeroRow]==0: # Muokataan vain, jos neulasissa 0-rivejÃ¤
                    nonZeroRow = next((m for m, x in enumerate(sdata0['Age'][zeroRow:]) if x!=0), None)
                    if nonZeroRow is not None: # Jos 0-ikÃ¤ viimeisellÃ¤ rivillÃ¤, ei tehdÃ¤ tarkistuksia/korjauksia
                        nonZeroRow = nonZeroRow + zeroRow
                        if sdata0['AIKA'][nonZeroRow]-sdata0['AIKA'][zeroRow] != sdata0['Age'][nonZeroRow]-sdata0['Age'][zeroRow]:
    
                            sdata0.loc[zeroRow, 'Age'] = int(sdata0['Age'][nonZeroRow] - (sdata0['AIKA'][nonZeroRow]-sdata0['AIKA'][zeroRow]))
                            with open(outfile, "a") as myfile:
                                myfile.write('\n0-age corrected to: ' + str(sdata0.loc[zeroRow, 'Age']) + ', i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
    
                    else:
                        with open(outfile, "a") as myfile:
                            myfile.write('\n0-age detected at last row!! i=' + str(i)+ ' kuvio:'+ str(kuviolista[i]))
    
                a = zeroRow + 1
            else:
                a = len(sdata0)+1     
                
                
    
    
        # Tsekataan TotVol ja DomH nollat ja muokataan
        # TÃ¤hÃ¤n tarvii ehkÃ¤ suuremman alarajan?
        # PitÃ¤Ã¤ lisÃ¤tÃ¤ myÃ¶s tyhjien tsekkaus!
        
        for k in range(0,len(sdata0)):
            
            if (sdata0.loc[k,'TotVol']<0.5) | (np.isnan(sdata0.loc[k,'TotVol']) == True): #voi olla 0.01 tilavuuksiakin
                
                if k < len(sdata0)-1:
                    sdata0.loc[k,'TotVol'] = max(0.01, 0.01*sdata0.loc[k+1,'TotVol']) # LisÃ¤tÃ¤Ã¤n totvol nollien tilalle muu pieni luku (0.001 oli liian pieni, kasvu ei lÃ¤htenyt kÃ¤yntiin, joten kasvatettu 0.2:een) paitsi et tulee 0.5
                else:
                    sdata0.loc[k,'TotVol'] = 0.01
                    
            if (sdata0.loc[k,'DomH']<0.3) | (np.isnan(sdata0.loc[k,'DomH']) == True):
                 
                 if k < len(sdata0)-1:
                     sdata0.loc[k,'DomH'] = max(0.3, 0.3*sdata0.loc[k+1,'DomH']) # LisÃ¤tÃ¤Ã¤n domh nollien tilalle muu pieni luku
                 else:
                     sdata0.loc[k,'DomH'] = 0.3       
       
        
        # Adding values to "missing" biomass components with similar increase-% per year as TotVol
        # Neulaset viimeisenÃ¤, jotta sitÃ¤ voidaan kÃ¤yttÃ¤Ã¤ arvioimaan, tarviiko biomassoja tÃ¤ydentÃ¤Ã¤.
        # Poistettu tÃ¤stÃ¤ listasta 'tuotosKuolleet', koska se on kumulatiivinen..
        
        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg', \
                     'TukPL1', 'TukPL2', 'TukPL3', 'TukPL4',\
                         'KuiPL1', 'KuiPL2', 'KuiPL3', 'KuiPL4',\
                         'HukPL1', 'HukPL2', 'HukPL3', 'HukPL4',\
                             'tuotosTilavuus',\
                                 'Runkopuukg', 'Hukkapuukg',\
                                     'MeanDbh','medianH', 'BA', 'tuotos']  
        #comps = ['DomH','MeanDbh','medianH',\
        #          'BA','Kannot_juuretkg',\
        #                  'tuotosTilavuus',\
        #                      'TukPL1', 'TukPL2','TukPL3',\
        #                          'KuiPL1', 'KuiPL2','KuiPL3',\
        #                              'HukPL1', 'HukPL2','HukPL3',\
        #                                  'VolPL1', 'VolPL2', 'VolPL3']  

        #biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg'] #tama lisatty tahan elokuu 23 jotta namakin tsekataan

        # if repairB==1: # Jos lisÃ¤tty 1-rivi  kopioimalla, korjataan sen komponentit ensin
        #     for j in range(0,len(biom_comp)):
        #         comp = biom_comp[j]
        #         sdata0.loc[1, comp] = sdata0['TotVol'][1]/sdata0['TotVol'][2]*sdata0.loc[2, comp]    
            
        # if repair==1: # Jos lisÃ¤tty 0-rivi kopioimalla, korjataan sen komponentit ensin
        #     for j in range(0,len(biom_comp)):
        #         comp = biom_comp[j]
        #         sdata0.loc[0, comp] = sdata0['TotVol'][0]/sdata0['TotVol'][1]*sdata0.loc[1, comp]
                
                            
    
                
        biom_comp = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg']
    
        for j in range(0,len(biom_comp)):
        
            comp = biom_comp[j]
            k=0
    
            while k < len(sdata0)-1: 
                
                zeroRow = next((m for m, x in enumerate(sdata0[comp][k:]) if x==0), None)
    
                if zeroRow is not None: # If zeroRow exists
                    
                    zeroRow = zeroRow + k
                    
                    # if sdata0['Neulasetkg'][zeroRow]==0: # Muiden komponenttien numeroita muokatan vain, jos neulasissa on nollaa. TÃ¤mÃ¤ poistettu, koska nÃ¤issÃ¤ komponenteissa voi olla nollaa, vaikka neulasmassa ei olisikaan nollaa..
                    nonZeroRow = next((m for m, x in enumerate(sdata0[comp][zeroRow:]) if x!=0), None)
                    
                    if nonZeroRow is not None:
                        
                        nonZeroRow = nonZeroRow + zeroRow
                        sdata0.loc[zeroRow, comp] = (sdata0['TotVol'][zeroRow]/sdata0['N'][zeroRow]*sdata0['Age'][zeroRow])/(sdata0['TotVol'][nonZeroRow]/sdata0['N'][nonZeroRow]*sdata0['Age'][nonZeroRow])*sdata0.loc[nonZeroRow, comp] # tää vaikuttaa siihen myös et 0 riville tulee 0.5 vol ja seuraava on jotain
                    
                        k=k+1
                        
                    else:
                        
                        k=k+1                          
    
                else:
                    k=k+1
                            
    
                        
                        # else:
                            
                            # Jotta ei toistu jokaisella komponentilla erikseen, 
                            # nÃ¤mÃ¤ korjataan vasta lopuksi: 
                            # kopioidaan loppuun rivit, jotka vastaavat pÃ¤Ã¤hakkuun jÃ¤lkeistÃ¤ kehitystÃ¤
                            
                            # EtsitÃ¤Ã¤n toisen pÃ¤Ã¤hakkuun jÃ¤lkeinen rivi
            
    
             
                
             
    
    
    
        # TehdÃ¤Ã¤n sarake 'treesp'=pÃ¤Ã¤puulaji volumen perusteella
        
        sdata0['treesp'] = pd.Series([], dtype='int64') #        sdata0['treesp'] = pd.Series([])

        
        for t in range(0,len(sdata0)):
            
            
            volPine = sdata0['VolPL1'][t]
            volSpruce = sdata0['VolPL2'][t]
            volBirch = sdata0['VolPL3'][t] + sdata0['VolPL4'][t]
        
            spe = 3
        
            if volPine >= max(volSpruce, volBirch):
                spe = 1
            
            if volSpruce > max(volPine, volBirch):
                spe = 2    
    
            sdata0.loc[t, 'treesp']=spe
    
                            
                            
        # Tsekataan, lÃ¶ytyykÃ¶ viimeisen pÃ¤Ã¤tehakkuun jÃ¤lkeen pelkkiÃ¤ 0-rivejÃ¤. Jos nÃ¤in,
        # kopioidaan aiemmat pÃ¤Ã¤tehakkuun jÃ¤lkeiset rivit loppuun.
        # Korjataan sen jÃ¤lkeen myÃ¶s AIKA-sarake
    
        if len(paateh)>0:
            paate1 =  next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>paateh[0]), None) 
            nonZeroRow = next((m for m, x in enumerate(sdata0['Neulasetkg'][paate1:]) if x!=0), None)
            
            if len(paateh)>1: # Jos useampi kuin yksi pÃ¤Ã¤tehakkuu
                paaterow = next((m for m, x in enumerate(sdata0['AIKA'][0:]) if x>=paateh[1]), None) # tÃ¤ssÃ¤ aiemmin x>paateh[1] !!
                nonZeroRow = next((m for m, x in enumerate(sdata0['Neulasetkg'][paaterow:]) if x!=0), None)
            
                if nonZeroRow is None:
                    
                    paate2aika = sdata0['AIKA'][paaterow]
                    
                    b = 0
                    
                    for p in range(paaterow,len(sdata0)):
                        print(p, paate2aika, "paatehakkuun jalkeiset nollarivit")
    
                        currow = paate1 + p - paaterow -1 # LisÃ¤tty tÃ¤hÃ¤n -1, koska rivillÃ¤ 617 lisÃ¤tty x>=paateh[1]
    
                        if currow<b: # LisÃ¤ys 17.10: Joissakin tapauksissa yllÃ¤ toimii -1, toisissa 0, joten muokataan currow tarvittaessa...
                            currow=b
                            
                        sdata0.loc[p]=sdata0.loc[currow]
                        
                        if b==0:
                            sdata0.loc[p, 'AIKA'] = int(paate2aika)
                            sdata0.loc[p, 'tuotosKuolleet'] = sdata0.loc[p-1, 'tuotosKuolleet']
                            b=b+1
                        else:
                            sdata0.loc[p, 'AIKA'] = int(paate2aika + sdata0.loc[currow, 'AIKA'] - sdata0.loc[currow-b, 'AIKA'])
                            sdata0.loc[p, 'tuotosKuolleet'] = sdata0.loc[p-1, 'tuotosKuolleet']
                            b=b+1
                
            else:
                # Jos vain yksi pÃ¤Ã¤tehakkuu, ei voida kopioida rivejÃ¤ aiemman pÃ¤Ã¤tehakkuun jÃ¤lkeen..
                print('MitÃ¤s tÃ¤hÃ¤n? kuvio:', str(kuviolista[i]),paateh)   
                if nonZeroRow is None:
                    with open(outfile, "a") as myfile:
                        myfile.write('\nPÃ¤Ã¤tehakkuun jÃ¤lkeisiÃ¤ 0-rivejÃ¤ ei saada korjattua! i=' + str(i)+ ' kuvio:'+ str(kuviolista[i])+ ' paateh: '+str(paateh)) #Ã¤Ã¤ 
    
        
        # Tsekataan, lÃ¶ytyykÃ¶ viimeiseltÃ¤ riviltÃ¤ (aika=100), nollariviÃ¤.
        # Jos lÃ¶ytyy, poistetaan.
        
    
        s100 = sdata0[(sdata0['AIKA']==100) & (sdata0['Age']==0)] # Vai vol=0?
        
        if len(s100)>0:
            sdata0 = sdata0[sdata0['AIKA']!=100]
            
        # -------------------------------------------------------------------
        # LisÃ¤ys 13.9.2021
        # Tarkistetaan, olisiko ensimmÃ¤inen Motti-tiedosto 1-rivin mittainen:
        
        vuodet_p = vuodet[vuodet>0]
        
        if (len(vuodet_p)>0):
            
            sdata = sdata0[sdata0['AIKA'] < vuodet_p[0]]
            
            if len(sdata)==1:
                #sdata0 = sdata0.append([sdata0.iloc[0]], ignore_index=True) # Rivi menee taulukon loppuun
                nolla2_df=pd.DataFrame([sdata0.iloc[0]])
                sdata0 = pd.concat((sdata0,nolla2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                sdata0['AIKA'][len(sdata0)-1] = vuodet_p[0]-1
        
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                sdata0.loc[1,'TotVol'] = tulotdata['volume'+str(int(sdata0['AIKA'][1]))][0]
                sdata0.loc[1,'Age'] = sdata0['Age'][0]+(sdata0['AIKA'][1]-sdata0['AIKA'][0])
    
                for c in range(0,len(comps)):
                    if sdata0.loc[0,comps[c]]>0:
                        # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
                        # sitÃ¤ voisi supistaa..??
                        sdata0.loc[1,comps[c]] = (sdata0['TotVol'][1]/sdata0['N'][1]*sdata0['Age'][1])/(sdata0['TotVol'][0]/sdata0['N'][0]*sdata0['Age'][0])*sdata0.loc[0, comps[c]]
    
                for c in range(0,len(biom_comp)):
                    if sdata0.loc[0,biom_comp[c]]>0:
                        # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
                        # sitÃ¤ voisi supistaa..??
                        sdata0.loc[1,biom_comp[c]] = (sdata0['TotVol'][1]/sdata0['N'][1]*sdata0['Age'][1])/(sdata0['TotVol'][0]/sdata0['N'][0]*sdata0['Age'][0])*sdata0.loc[0, biom_comp[c]]
    
    
        # -------------------------------------------------------------------
        # LisÃ¤ys 9.3.2022
        # Tarkistetaan, olisiko *toinen* Motti-tiedosto 1-rivin mittainen:
        
        #if (len(vuodet_p)>0):
        #    
        #    if len(vuodet_p)==1:
        #        end_year=101
        #    else:
        #        end_year = vuodet_p[1]
        #        
        #    sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[0]) & (sdata0['AIKA'] < end_year)] # Huom. tÃ¤ssÃ¤ ei saa resetoida indeksiÃ¤, jotta jatko toimii
        #    
        #    if len(sdata)==1:
        #        sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # LisÃ¤tÃ¤Ã¤n rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
        #        sdata0['AIKA'][len(sdata0)-1] = end_year - 1
        # 
        #         sdata0 = sdata0.sort_values(by=['AIKA'])
        #         sdata0 = sdata0.reset_index(drop=True)
        #         
        #        sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
        #        
        #        sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
        #        
        #        if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa riviÃ¤
        #            sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
        #
        #            for c in range(0,len(comps)):
        #                if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
        #                    # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
        #                    # sitÃ¤ voisi supistaa..??
        #                    sdata0.loc[sdatax.index[0],comps[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, comps[c]]
        # 
        #            for c in range(0,len(biom_comp)):
        #                if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
        #                    # Alla BA pÃ¤ivittyy liian isoksi, mutta ei ole perusteita, joilla
        #                    # sitÃ¤ voisi supistaa..??
        #                    sdata0.loc[sdatax.index[0],biom_comp[c]] = sdata0['TotVol'][sdatax.index[0]]/sdata0['TotVol'][sdatax.index[0]-1]*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
        
        
        # -------------------------------------------------------------------
        # Update 6.3.2023: checking all 1-row-motti-files except the last one
        
        for v in range(1, len(vuodet_p)):
            
            if (len(vuodet_p)>0):
                
                if len(vuodet_p)==1:
                    end_year=101
                else:
                    end_year = vuodet_p[v]                
                
                sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[v-1]) & (sdata0['AIKA'] < end_year)] # Huom. tÃ¤ssÃ¤ ei saa resetoida indeksiÃ¤ jotta jatko toimii
                
                if len(sdata)==1:
                    #sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # LisÃ¤tÃ¤Ã¤n rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
                    sdataind_df=pd.DataFrame([sdata0.iloc[sdata.index[0]]])
                    sdata0 = pd.concat((sdata0,sdataind_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                    sdata0['AIKA'][len(sdata0)-1] = end_year - 1
                    
                    sdata0 = sdata0.sort_values(by=['AIKA'])
                    sdata0 = sdata0.reset_index(drop=True)
                    
                    sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
                    
                    sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
                    
                    if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa rivia
                        sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
                        
                        for c in range(0,len(comps)):
                            if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
                                # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                                # sita voisi supistaa..??
                                sdata0.loc[sdatax.index[0],comps[c]] = (sdata0['TotVol'][sdatax.index[0]]/sdata0['N'][sdatax.index[0]]*sdata0['Age'][sdatax.index[0]])/(sdata0['TotVol'][sdatax.index[0]-1]/sdata0['N'][sdatax.index[0]-1]*sdata0['Age'][sdatax.index[0]-1])*sdata0.loc[sdatax.index[0]-1, comps[c]]
                        
                        for c in range(0,len(biom_comp)):
                            if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
                                # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                                # sita voisi supistaa..??
                                sdata0.loc[sdatax.index[0],biom_comp[c]] = (sdata0['TotVol'][sdatax.index[0]]/sdata0['N'][sdatax.index[0]]*sdata0['Age'][sdatax.index[0]])/(sdata0['TotVol'][sdatax.index[0]-1]/sdata0['N'][sdatax.index[0]-1]*sdata0['Age'][sdatax.index[0]-1])*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
        
        # Update 6.3.2023: Checking if the last motti-file would be a 1-row-file
        
        if (len(vuodet_p)>1):
            
            sdata = sdata0[(sdata0['AIKA'] >= vuodet_p[len(vuodet_p)-1])] # Huom. tassa ei saa resetoida indeksia, jotta jatko toimii
            end_year = 101
            
            if (len(sdata)==1) & (sdata0['AIKA'][len(sdata0)-1]<100):
                
                #sdata0 = sdata0.append([sdata0.iloc[sdata.index[0]]], ignore_index=True) # Lisataan rivi sdata:n index-numeron perusteella, Rivi menee taulukon loppuun
                sdataind2_df=pd.DataFrame([sdata0.iloc[sdata.index[0]]])
                sdata0 = pd.concat((sdata0, sdataind2_df), axis=0, ignore_index=True) # Rivi menee taulukon loppuun
                    
                sdata0['AIKA'][len(sdata0)-1] = end_year - 1
                
                sdata0 = sdata0.sort_values(by=['AIKA'])
                sdata0 = sdata0.reset_index(drop=True)
                
                sdatax = sdata0[(sdata0['AIKA'] == end_year - 1)]
                
                sdata0.loc[sdatax.index[0],'TotVol'] = tulotdata['volume'+str(int(end_year-1))][0]
                
                if sdatax.index[0]>0: # Muussa tapauksessa tulee kaksi samaa rivia
                    sdata0.loc[sdatax.index[0],'Age'] = sdata0['Age'][sdatax.index[0]]+(sdata0['AIKA'][sdatax.index[0]]-sdata0['AIKA'][sdatax.index[0]-1])
                    
                    for c in range(0,len(comps)):
                        if sdata0.loc[sdatax.index[0]-1,comps[c]]>0:
                            # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                            # sita voisi supistaa..??
                            sdata0.loc[sdatax.index[0],comps[c]] = (sdata0['TotVol'][sdatax.index[0]]/sdata0['N'][sdatax.index[0]]*sdata0['Age'][sdatax.index[0]])/(sdata0['TotVol'][sdatax.index[0]-1]/sdata0['N'][sdatax.index[0]-1]*sdata0['Age'][sdatax.index[0]-1])*sdata0.loc[sdatax.index[0]-1, comps[c]]
                    
                    for c in range(0,len(biom_comp)):
                        if sdata0.loc[sdatax.index[0]-1,biom_comp[c]]>0:
                            # Alla BA paivittyy liian isoksi, mutta ei ole perusteita, joilla
                            # sita voisi supistaa..??
                            sdata0.loc[sdatax.index[0],biom_comp[c]] = (sdata0['TotVol'][sdatax.index[0]]/sdata0['N'][sdatax.index[0]]*sdata0['Age'][sdatax.index[0]])/(sdata0['TotVol'][sdatax.index[0]-1]/sdata0['N'][sdatax.index[0]-1]*sdata0['Age'][sdatax.index[0]-1])*sdata0.loc[sdatax.index[0]-1, biom_comp[c]]
                            
                
        
        # -------------------------------------------------------------------
        # LisÃ¤ys 9.3.2022
        # Tarkistetaan, tuleeko perÃ¤kkÃ¤isten harvennusten takia liian paljon
        # vÃ¤heneviÃ¤ tilavuuksia 
        
        comps_notk = ['Runkopuukg', 'Hukkapuukg', 'Oksatkg', 'kannot', '>2mm_juuretkg', 'h_juuretkg', 'DomH','Neulasetkg', \
                             'tuotosTilavuus','MeanDbh','medianH', 'BA', 'tuotos', 'TotVol']  
        
        
        if pera >= 0:
            
            vuodet_p = vuodet[vuodet>=0] # LisÃ¤ys 17.10.2022: vuodet>=0 (jos tÃ¤ssÃ¤ olisi vuodet>0, vuodet_p saattaisi jÃ¤Ã¤dÃ¤ tyhjÃ¤ksi -> virhe seuraavalla rivillÃ¤)
            
            sdata = sdata0[sdata0['AIKA'] < vuodet_p[0]]
            
            if (len(sdata)<=2) & (len(sdata)>0): # LisÃ¤ys 17.10.2022: len(sdata)>0, koska jos sdata on tyhjÃ¤ -> error
                
                if (sdata['TotVol'][0] - sdata['TotVol'][1]) > 0:
                    
                    sdata0.loc[0]=sdata0.loc[1]
                    sdata0.loc[0, 'AIKA'] = sdata0['AIKA'][1] -1
                    sdata0.loc[0, 'Age'] = sdata0['Age'][1] -1
                    
                    
                    for c in range(0,len(comps_notk)):
                        
                        comp = comps_notk[c]
                        sdata0.loc[0, comp] = sdata0.loc[0, comp] * 0.97 # Arvio, jotta 
                                            
    
    
        
        
        # Muutetaan Age=-1 --> Age=0, nÃ¤itÃ¤ lÃ¶ytyi suoraan Motti-datasta
        
        sdata0['Age'] = np.clip(sdata0['Age'], 0, np.nan) 
        
        # -------------------------------------------------------------------
        # Lisays 28.3.2023:
        # Tarkistetaan, loytyyko iasta perakkaisia nollia
        
        a = 0
        while a < len(sdata0):
            zeroRow = next((m for m, x in enumerate(sdata0['Age'][a:]) if x==0), None)
            
            if zeroRow is not None:
                zeroRow = zeroRow+a
                try:
                    if sdata0['Age'][zeroRow+1]==0:
                        sdata0.loc[zeroRow, 'Age'] = 1
                        sdata0.loc[zeroRow+1, 'Age'] = 1 + int(sdata0['AIKA'][zeroRow+1]-sdata0['AIKA'][zeroRow])
                        with open(outfile, "a") as myfile:
                            myfile.write('\nDouble zero ages!! i=' + str(i))
                        
                    a = zeroRow + 1
                except:
                    a = zeroRow + 1
                        
            else:
                a = len(sdata0)+1     
        
        
        
        
        # Varsinaisten Motti-tiedostojen teko alkaa tÃ¤stÃ¤
        
        if len(sdata0)>0:
            
            if len(vuodet)>0:
                
                for j in range(0,n_mottifiles):
                
                    if j == 0:
                        sdata = sdata0[sdata0['AIKA'] < vuodet[j]]
                        if len(sdata)==0:
                            continue
                    
                    if (j>0) & (j<n_mottifiles-1):
                        sdata = sdata0[(sdata0['AIKA'] >= vuodet[j-1]) & (sdata0['AIKA'] < vuodet[j])]
                        sdata = sdata.reset_index(drop=True)
                        if len(sdata)==0:
                            continue
                    
                    if j == n_mottifiles-1:
                        sdata = sdata0[sdata0['AIKA'] >= vuodet[j-1]]
                        sdata = sdata.reset_index(drop=True)
                        if len(sdata)==0:
                            continue
                    
                    kasvatus = pd.Series(np.ones(len(sdata)))
                    vuosi = sdata['AIKA']
                    ika = sdata['Age']
                    N = sdata['N']
                    PPA = sdata['BA']
                    Hg = sdata['medianH']
                    Dg = sdata['MeanDbh']
                    Hdom = sdata['DomH']
                    Tilavuus = sdata['TotVol']
                    Tukki = sdata['TukPL1']+sdata['TukPL2']+sdata['TukPL3']+sdata['TukPL4']
                    Kuitu = sdata['KuiPL1']+sdata['KuiPL2']+sdata['KuiPL3']+sdata['KuiPL4']
                    Hukka = sdata['HukPL1']+sdata['HukPL2']+sdata['HukPL3']+sdata['HukPL4']
                    Tuotos = sdata['tuotosKuolleet'] + sdata['tuotosTilavuus']
                    Kuolleisuus = Tuotos*0. # Ei kÃ¤ytetÃ¤ Susissa
                    runko_aines = sdata['Runkopuukg']/1000. # YksikkÃ¶ tn
                    runko_hukka = sdata['Hukkapuukg']/1000.
                    elavat_oksat = sdata['Oksatkg']*0.8/1000. # Oksatkg = elÃ¤vÃ¤t ja kuolleet oksat yhteensÃ¤
                    kuolleet_oksat = sdata['Oksatkg']*0.2/1000.
                    lehdet = sdata['Neulasetkg']/1000.
                    kannot = sdata['kannot']/1000.
                    juuret = sdata['>2mm_juuretkg']/1000.
                    hienojuuret = sdata['h_juuretkg']/1000.
                    
                    # Tsekkaa, ettÃ¤ tÃ¤mÃ¤ on ok!!!
                    for t in range(0,len(sdata)):
                        if t==0:
                            Tuotos[0]=sdata['tuotosTilavuus'][0]
                        
                        else:
                            Tuotos[t] = sdata['tuotosTilavuus'][t] + (sdata['tuotosKuolleet'][t] - sdata['tuotosKuolleet'][t-1])
                    
                    # ---tÃ¤hÃ¤n asti
                    
                    data = pd.DataFrame(np.transpose(np.array([kasvatus, vuosi, ika, N, PPA, Hg, Dg, Hdom, Tilavuus, Tukki, Kuitu, \
                                         Hukka, Tuotos, Kuolleisuus, runko_aines, runko_hukka, \
                                         elavat_oksat, kuolleet_oksat, lehdet, kannot, juuret, hienojuuret])), \
                                         columns=['Kasvatus', 'Vuosi', 'Ikä', 'N', 'PPA', 'Hg', 'Dg', 'Hdom', 'Tilavuus', 'Tukki', 'Kuitu',\
                                                  'Hukka', 'Tuotos', 'Kuolleisuus', 'runko(aines)', 'runko(hukka)', \
                                                  'elävät oksat', 'kuolleet oksat', 'lehdet', 'Kannot', 'Juuret >2mm', 'Hienojuuret'])
                    
                    standInfo = kuviot[kuviot['KUVIO']==kuviolista[i]]
                    standInfo = standInfo.reset_index(drop=True)
                    
                    """ TSEKKAA, ETTÃ„ PUULAJI TULEE OIKEIN! """
                    
                    
                    
                    if (project == 'life') | (project == 'hiilipolku'):
                        spe = sdata['treesp'][len(sdata)-1]
                    
                    else: 
                        if (project =='vmi') | (project=='suo'):
                            
                            spe = int(sdata['treesp'][len(sdata)-1]) # viimeisen rivin puulaji
                            # treesp=kdata['Selite'][0][28:33]
                            # spe = 3
                            
                            # if treesp=='MÃ¤nty':
                            #     spe = 1
                            # if treesp=='Kuusi':
                            #     spe = 2       
                                
                            
                        else:
                            
                            treesp = pdata2['Selite'][0][28:33]
                            
                            spe = 3
                            
                            if treesp=='Mänty': #onk tÃ¤ a va Ã¤#
                                spe = 1
                            if treesp=='Kuusi':
                                spe = 2 
                    
                    
                    # print(kuviolista[i], ': puulaji', spe)
                    with open(outfile, "a") as myfile:
                        myfile.write('\n' + str(kuviolista[i]) +  ': puulaji ' + str(spe)+ ' n'+str(j))
                    
                    
                    fout = mottifolder + '/' + str(kuviolista[i]) + '_n' + str(j) + '.xls'
                    
                    cols = pd.DataFrame([['','','','',spe,'','','','','','','','','']],columns=['Kasvatus',	'Vuosi',	'id Harvennus',	'Harvennus',	'id Puulaji',	'Puulaji',	'Tukki[m³/ha]',	'Pikkutukki[m³/ha]',	'Kuitu[m³/ha]',	'Energiapuu, runko(aines)[m³/ha]',	'Energiapuu, runko(hukka)[m³/ha]',	'Energiapuu, oksat(e)[m³/ha]',	'Energiapuu, oksat(k)[m³/ha]',	'Energiapuu, kannot ja juuret[m³/ha]'])
                    
                    writer = pd.ExcelWriter(fout, engine = 'xlsxwriter')
                    data.to_excel(writer, index=False, sheet_name = 'Puustotunnukset')
                    cols.to_excel(writer, index=False, sheet_name = r'Kertymät')
                    #writer.save()
                    writer.close()
            else:
                # print(str(kuviolista[i]) + ': poistumat puuttuvat! i=' + str(i))
                with open(outfile, "a") as myfile:
                    myfile.write('\n' + str(kuviolista[i]) + ': poistumat puuttuvat! i=' + str(i))
                
                sdata = sdata0.copy()
                
                kasvatus = pd.Series(np.ones(len(sdata)))
                vuosi = sdata['AIKA']
                ika = sdata['Age']
                N = sdata['N']
                PPA = sdata['BA']
                Hg = sdata['medianH']
                Dg = sdata['MeanDbh']
                Hdom = sdata['DomH']
                Tilavuus = sdata['TotVol']
                Tukki = sdata['TukPL1']+sdata['TukPL2']+sdata['TukPL3']+sdata['TukPL4']
                Kuitu = sdata['KuiPL1']+sdata['KuiPL2']+sdata['KuiPL3']+sdata['KuiPL4']
                Hukka = sdata['HukPL1']+sdata['HukPL2']+sdata['HukPL3']+sdata['HukPL4']
                Tuotos = sdata['tuotos']
                Kuolleisuus = Tuotos*0. # Ei kÃ¤ytetÃ¤ Susissa
                runko_aines = sdata['Runkopuukg']/1000. # YksikkÃ¶ tn
                runko_hukka = sdata['Hukkapuukg']/1000.
                elavat_oksat = sdata['Oksatkg']*0.8/1000. # Oksatkg = elÃ¤vÃ¤t ja kuolleet oksat yhteensÃ¤
                kuolleet_oksat = sdata['Oksatkg']*0.2/1000.
                lehdet = sdata['Neulasetkg']/1000.
                kannot = sdata['kannot']/1000.
                juuret = sdata['>2mm_juuretkg']/1000.
                hienojuuret = sdata['h_juuretkg']/1000.
                
                data = pd.DataFrame(np.transpose(np.array([kasvatus, vuosi, ika, N, PPA, Hg, Dg, Hdom, Tilavuus, Tukki, Kuitu, \
                                     Hukka, Tuotos, Kuolleisuus, runko_aines, runko_hukka, \
                                     elavat_oksat, kuolleet_oksat, lehdet, kannot, juuret, hienojuuret])), \
                                     columns=['Kasvatus', 'Vuosi', 'Ikä', 'N', 'PPA', 'Hg', 'Dg', 'Hdom', 'Tilavuus', 'Tukki', 'Kuitu',\
                                              'Hukka', 'Tuotos', 'Kuolleisuus', 'runko(aines)', 'runko(hukka)', \
                                              'elävät oksat', 'kuolleet oksat', 'lehdet', 'Kannot', 'Juuret >2mm', 'Hienojuuret'])
                
                standInfo = kuviot[kuviot['KUVIO']==kuviolista[i]]
                standInfo = standInfo.reset_index(drop=True)
                
                
                """ TSEKKAA, ETTÃ„ PUULAJI TULEE OIKEIN! """
                if (project == 'life') | (project == 'hiilipolku'):
                    spe = sdata['treesp'][len(sdata)-1]
                
                else:
                    
                    if project =='vmi':
                        
                        # spe = int(sdata['treesp'][len(sdata)-1]) # viimeisen rivin puulaji
                       
                        treesp=kdata['Selite'][0][28:33]
                        spe = 3
                        
                        if treesp=='Mänty': #Ã¤
                            spe = 1
                        if treesp=='Kuusi':
                            spe = 2    
                    else:
                        treesp = pdata2['Selite'][0][28:33]
                        
                        spe = 3
                        
                        if treesp=='Mänty':
                            spe = 1
                        if treesp=='Kuusi':
                            spe = 2 
                    
                # print(kuviolista[i], ': puulaji ', spe)
                with open(outfile, "a") as myfile:
                    myfile.write('\n' + str(kuviolista[i]) +  ':' + str(spe)+ ' n0')
                
                
                
                fout = mottifolder + '/' + str(kuviolista[i]) + '_n0.xls'
                
                cols = pd.DataFrame([['','','','',spe,'','','','','','','','','']],columns=['Kasvatus',	'Vuosi',	'id Harvennus',	'Harvennus',	'id Puulaji',	'Puulaji',	'Tukki[m³/ha]',	'Pikkutukki[m³/ha]',	'Kuitu[m³/ha]',	'Energiapuu, runko(aines)[m³/ha]',	'Energiapuu, runko(hukka)[m³/ha]',	'Energiapuu, oksat(e)[m³/ha]',	'Energiapuu, oksat(k)[m³/ha]',	'Energiapuu, kannot ja juuret[m³/ha]'])
                
                writer = pd.ExcelWriter(fout, engine = 'xlsxwriter')
                data.to_excel(writer, index=False, sheet_name = 'Puustotunnukset')
                cols.to_excel(writer, index=False, sheet_name = r'Kertymät')
                #writer.save()
                writer.close()

###########################

scenarios=['Halvanjoki_BAU_A', 'Halvanjoki_BAU_B', 'Halvanjoki_BIO_A', 'Halvanjoki_BIO_B', 'Halvanjoki_HII_A', 'Halvanjoki_HII_B']
create_motti_files_silvi_lp4('hiilipolku',scenarios[0])
create_motti_files_silvi_lp4('hiilipolku',scenarios[1])
create_motti_files_silvi_lp4('hiilipolku',scenarios[2])                
create_motti_files_silvi_lp4('hiilipolku',scenarios[3])                
create_motti_files_silvi_lp4('hiilipolku',scenarios[4])                
create_motti_files_silvi_lp4('hiilipolku',scenarios[5])                

scenarios=['Halvanjoki_BAU_A', 'Halvanjoki_BAU_B', 'Halvanjoki_BIO_A', 'Halvanjoki_BIO_B', 'Halvanjoki_HII_A', 'Halvanjoki_HII_B']
create_motti_files_silvi_lp4_vol01('hiilipolku',scenarios[0])
create_motti_files_silvi_lp4_vol01('hiilipolku',scenarios[1])
create_motti_files_silvi_lp4_vol01('hiilipolku',scenarios[2])                
create_motti_files_silvi_lp4_vol01('hiilipolku',scenarios[3])                
create_motti_files_silvi_lp4_vol01('hiilipolku',scenarios[4])                
create_motti_files_silvi_lp4_vol01('hiilipolku',scenarios[5])                

create_motti_files_silvi_lp4_vol001('hiilipolku',scenarios[0])
create_motti_files_silvi_lp4_vol001('hiilipolku',scenarios[1])
create_motti_files_silvi_lp4_vol001('hiilipolku',scenarios[2])                
create_motti_files_silvi_lp4_vol001('hiilipolku',scenarios[3])                
create_motti_files_silvi_lp4_vol001('hiilipolku',scenarios[4])                
create_motti_files_silvi_lp4_vol001('hiilipolku',scenarios[5])                

create_motti_files_silvi_lp4_volnage('hiilipolku',scenarios[0])
create_motti_files_silvi_lp4_volnage('hiilipolku',scenarios[1])
create_motti_files_silvi_lp4_volnage('hiilipolku',scenarios[2])                
create_motti_files_silvi_lp4_volnage('hiilipolku',scenarios[3])                
create_motti_files_silvi_lp4_volnage('hiilipolku',scenarios[4])                
create_motti_files_silvi_lp4_volnage('hiilipolku',scenarios[5])                



scenarios=['Kuonanjoki_BAU_A', 'Kuonanjoki_BAU_B', 'Kuonanjoki_BIO_A', 'Kuonanjoki_BIO_B', 'Kuonanjoki_HII_A', 'Kuonanjoki_HII_B']
create_motti_files_silvi('hiilipolku',scenarios[0])
#create_motti_files_silvi('hiilipolku',scenarios[1])
#create_motti_files_silvi('hiilipolku',scenarios[2])                
#create_motti_files_silvi('hiilipolku',scenarios[3])                
#create_motti_files_silvi('hiilipolku',scenarios[4])                
#create_motti_files_silvi('hiilipolku',scenarios[5])                

scenarios=['Halvanjoki_BAU_A', 'Halvanjoki_BAU_B', 'Halvanjoki_BIO_A', 'Halvanjoki_BIO_B', 'Halvanjoki_HII_A', 'Halvanjoki_HII_B']
create_motti_files_silvi('hiilipolku',scenarios[0])
create_motti_files_silvi('hiilipolku',scenarios[1])
create_motti_files_silvi('hiilipolku',scenarios[2])                
create_motti_files_silvi('hiilipolku',scenarios[3])                
create_motti_files_silvi('hiilipolku',scenarios[4])                
create_motti_files_silvi('hiilipolku',scenarios[5])                
               
                
                