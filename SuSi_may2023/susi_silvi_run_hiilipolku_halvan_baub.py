# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 12:25:55 2021

@author: Leena Stenberg / LUKE

"""



import numpy as np
import pandas as pd
import datetime
from susi.susi_utils import  get_motti, read_FMI_weather, get_mese_input
#from inputs.susi_para import get_susi_para
from susi.susi_para import get_susi_para
import susi.susi_io
from susi.susi_main import Susi
import os




''' 
Tässä mukana:
    - simulointi katkeaa break_vol perusteella tai mikäli simuloitava aika tulee täyteen
    - ageSim ei voi olla nolla, muutetaan 0.001


How to use:
    1) run this file
    2) run call_local_susi_motti_silvi_list()

'''







''' Define input paths and input data: ------------------------------------'''
    
def get_paths():
    
    susiPath = r'/scratch/project_2002470/SUSI_HIILIPOLKU/susi/' # susi python files path
    outPath=r'/scratch/project_2002470/SUSI_HIILIPOLKU/outputs/Halvanjoki_BAU_B/' # outputs path
    wpath = r'/scratch/project_2002470/SUSI_HIILIPOLKU/inputs/' # weather data path
    #mottifolder = r'/scratch/project_2002470/SUSI_HIILIPOLKU/Halvanjoki_BAU_B/'
    mottipath = r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/Halvanjoki_BAU_B/'
    outfile = r'/scratch/project_2002470/SUSI_HIILIPOLKU/outputs/Halvanjoki_BAU_B/Halvanjoki_BAU_B.txt'    
    if not os.path.exists(outPath):
        os.mkdir(outPath)    
    return susiPath, wpath, mottipath, outPath, outfile

def get_motti_data():
        
    kuviot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/Halvanjoki_BAU_B/motti/Halvanjoki_BAU_B_kuviot.csv', encoding='latin1', sep=';')
    puustot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/Halvanjoki_BAU_B/motti/Halvanjoki_BAU_B_puustot.csv', encoding='latin1', sep=';')
    poistumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/Halvanjoki_BAU_B/motti/Halvanjoki_BAU_B_poistumat.csv', encoding='latin1', sep=';')
    tapahtumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/Halvanjoki_BAU_B/motti/Halvanjoki_BAU_B_tapahtumat.csv', encoding='latin1', sep=';', index_col=False)
    tulot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Halvanjoki/Halvanjoki_BAU_B/motti/Halvanjoki_BAU_B_kasvut.csv', encoding='latin1', sep=';')
    
    return kuviot, puustot, poistumat, tapahtumat, tulot

def get_inputs():
        
    wdata='halvansuo_weathercont.csv' #'halvansuo_weather.csv' #'puruvesi_saa_19810101-2020062.csv' #'halvansuo_weather.csv' #'kuonanjoki_weather.csv'# # # Weather data file name
    sarka = 40. # Ditch spacing
    n_ditch_scens = 4 # Number of ditch scenarios 
    start_yr_ini = 1991 # Initial starting year
    totSimYears=50 # Total lenght of simulations, in years
    
    return wdata, sarka, n_ditch_scens, start_yr_ini, totSimYears

''' -----------------------------------------------------------------------'''





def get_sfc_spec(sfc, puulaji):
    
    ''' Yksinkertainen pikaversio. Tätä voisi parantaa tarkemmilla puulajisuhteilla.
    Puulajisuhteet pitää ottaa sellaisesta vaiheesta, jossa puusto on kuviolla jo reilusti, esim. 80 m3/ha (tähän ei joka kuviolla välttämättä päästä)    
    Soili: Mtkg on helppo: I on kuusikko ja II männikkö (hiestä voi olla). Ptkg on vaikeampi: I on yleensä puhdas tai lähes puhdas männikkö, II:lla on yleensä koivua suht reilusti sekapuuna (ellei hakattu pois). 
    '''
    
    sfc_spec = 2 # Oletus
    
    if (sfc==3) & (puulaji=='Kuusi'):
        sfc_spec = 1
        
    if (sfc==3) & (puulaji=='Mänty'):
        sfc_spec = 2
        
    if (sfc==4) & (puulaji=='Mänty'):
        sfc_spec = 1
        
    return sfc_spec


def get_motti_silvi(ifile, return_spe = False):
    #---read the Motti-simulation to be used as a basis for the Susi-simulation

    cnames=['yr', 'age', 'N', 'BA', 'Hg', 'Dg', 'hdom', 'vol', 'logs', 'pulp', 'loss', 'yield','mortality', 
            'stempulp', 'stemloss', 'branch_living', 'branch_dead', 'leaves', 'stump', 'roots_coarse', 'roots_fine']
    df = pd.read_excel(ifile, sheet_name=0, usecols=range(22), skiprows=1, header=None )
    df = df.drop([0], axis=1)
    df.columns=cnames

    cname=['idSpe']
    df2 = pd.read_excel(ifile, sheet_name=1, usecols=[4], skiprows=1, header=None )
    df2.columns = cname
    
    #---- find thinnings and add a small time to lines with the age to enable interpolation---------
    # df = df[df['age']!=0]    
   
    steps = np.array(np.diff(df['age']), dtype=float)
    idx = np.ravel(np.argwhere(steps<1.))+1
    df['age'][idx]=df['age'][idx]+5./365.
    if return_spe: 
        return df, df2['idSpe'][0]
    else:
        return df


def check_end(n_ditch_scens, new_end_yr, end_yr):
    
    end_reached = True
    
    # Jos yhdessäkin ojasyvyysskenaariossa on simulointivuosia jäljellä, jatketaan simulointia    
    
    for e in range(0, n_ditch_scens): 
        if new_end_yr[e] < end_yr:
            if new_end_yr[e] > 0:
                end_reached = False
            break
        
    return end_reached


def call_local_susi_silvi(wdata, mottifile, sarkaSim, ageSim, sfc_mean, sfc_spec, site, break_vol, start_yr_arr, start_yr_ini, end_yr, vol_aft, susi_outfile, ash_year):

    susiPath, wpath, mottipath, folderName, outfile = get_paths()
    
    start_date = datetime.datetime(start_yr_ini,1,1) # for forc only
    end_date=datetime.datetime(end_yr,12,31)
    
    n = int(sarkaSim / 2)                       # number of computation nodes within the strip
    
    ageSimD = {'dominant': ageSim*np.ones(n),
              'subdominant': 0*np.ones(n),
              'under': 0*np.ones(n)}            # age of the stand in each node
    
    sfc =  np.ones(n, dtype=int)*sfc_mean       # site fertility class, assuming same value (sfc_mean) in every node
    
    forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
                
    wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                              folderName=folderName, hdomSim=None,  
                                                                              ageSim=ageSimD, sarkaSim=sarkaSim, sfc=sfc, 
                                                                              susiPath=None,
                                                                              n=n, ash_year=ash_year, sfc_spec=sfc_spec)
    outpara['netcdf']=susi_outfile   
    new_susi = Susi()        
    
    new_end_yr = new_susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr_arr, end_yr, wlocation = 'undefined', 
                            mottifile=mottifile, peat= 'other', photosite='All data', 
                            folderName=folderName,ageSim=ageSimD, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath, break_vol=break_vol)       
                 
    #os.rename(folderName + "susi"+str(i)+".nc", folderName + susi_outfile)
    
    return new_end_yr
    


''' -----------------------------------------------------------------------'''
# Jos hajautetaan CSC:llä:
def call_local_susi_motti_silvi_list(i):
# def call_local_susi_motti_silvi_list():
    
    # Get inputs:
    kuviot, puustot, poistumat, tapahtumat, tulot = get_motti_data()    
    wdata, sarka, n_ditch_scens, start_yr_ini, totSimYears = get_inputs() 
    susiPath, wpath, mottipath, folderName, outfile = get_paths()
    
    kuviot = kuviot.rename(columns={'kuvio':'KUVIO'})
    puustot = puustot.rename(columns={'kuvio':'KUVIO'})
    poistumat = poistumat.rename(columns={'kuvio':'KUVIO'})
    tapahtumat = tapahtumat.rename(columns={'kuvio':'KUVIO',
                                  '_2003': '2003',
                                  '_2004': '2004',
                                  '_2005': '2005',
                                  '_2001': '2001'})
    tulot = tulot.rename(columns={'kuvio':'KUVIO'})

    end_yr = start_yr_ini + totSimYears - 1
    start_yr = start_yr_ini*np.ones(n_ditch_scens, dtype=int)
    
    kuviolista = list(set(kuviot['KUVIO'][kuviot['kuivatustilanne']>=3])) # id number list for peatland stands
    
    # jos hajautetaan CSC:llä, poistetaan tästä for-lause ja i annetaan funktiolle inputtina
    # for i in range(0,len(kuviolista)):
     
    with open(outfile, "a") as myfile:
        myfile.write('\n\nStand=' + str(kuviolista[i]) + ', i=' + str(i))  
        
    stand = kuviolista[i]

    pdata = poistumat[(poistumat['KUVIO']==stand) & (poistumat['HARV']>1)]       
    pdata = pdata.reset_index(drop=True)
    
    pdata2 = poistumat[(poistumat['KUVIO']==stand)]
    pdata2 = pdata2.reset_index(drop=True)
    
    tdata = tapahtumat[(tapahtumat['KUVIO']==stand)]
    tdata = tdata.reset_index(drop=True)       
    
    tulotdata = tulot[(tulot['KUVIO']==stand)]
    tulotdata = tulotdata.reset_index(drop=True)
    
    ensih = tdata[tdata['2003']>-1]['2003'].values # ensiharvennus
    harv = tdata[tdata['2004']>-1]['2004'].values # harvennus
    paateh = tdata[tdata['2005']>-1]['2005'].values # päätehakkuu
    # varhaisp = tdata[tdata['2006']>-1]['2006'].values # varhaisperkaus
    # taimikonh = tdata[tdata['2007']>-1]['2007'].values # taimikonhoito
    
    '''
    Tuhkalannoitukselle pitää keksiä parempi ratkaisu. Tässä tulee suoraan Motti-datan vuoden mukaan.
    Parempi olisi sinä vuonna, jolloin vastaava ba tai vol saavutetaan. 
    '''
    
    ash = tdata[tdata['2001']>-1]['2001'].values
    # ash = np.array([1], dtype='int64') # Jos pakotetaan tuhkalannoitus tietylle vuodelle
    
    if len(ash)>0:
        ash_year = ash[0] + start_yr_ini # Oletetaan, että tuhkalannoitus vain kerran simulointiaikana
    else:
        ash_year=10000 # Jos ei tuhkalannoitusta -> niin iso luku, ettei sitä saavuteta
    
    vuodet = np.unique(np.sort(np.concatenate((ensih, harv, paateh))))
    
    # Otetaan perkkäisten vuosien poistettava vuosi talteen
    tuplat = pd.Series([], dtype=int)
    tup = 0
    
    # Tsekataan, onko tapahtumissa peräkkäisiä vuosia ja poistetaan niistä jälkimmäinen
    for v in range(0,len(vuodet)-1):
        if vuodet[v+1]-vuodet[v]==1:
            tuplat[tup]=vuodet[v+1]
            tup = tup+1
            vuodet = np.delete(vuodet, [v+1], None)
        if v+1>=len(vuodet)-1:
            break
    
    vol_before = pd.Series([], dtype=float) # Tämä ei sisällä tilavuuksia, jotka edeltävät vuosina 0 tai 1 tapahtumia! ..koska näitä vuosia ei simuloida eikä volumeja siten tavoitella
    vol_after = pd.Series([], dtype=float) # Tämä ei sisällä volumeja, jotka tulevat vuosien 0 tai 1 tapahtumien jälkeen!         
    
    tup=0
    
    for y in range(0,len(vuodet)):
        vol_before[y]=tulotdata['volume'+str(int(vuodet[y]))][0]
        if vuodet[y]<99: # aiemmin 100, 
            vol_after[y]=tulotdata['volume'+str(int(vuodet[y]+1))][0]
            if len(tuplat)>0:
                if tulotdata['volume'+str(int(vuodet[y]+2))][0] < tulotdata['volume'+str(int(vuodet[y]+1))][0]:
                    if tuplat[tup]-1==vuodet[y]:
                        tup=tup+1
                        vol_after[y]=tulotdata['volume'+str(int(vuodet[y]+2))][0]
        else:
            vol_after[y]=np.nan
                
    n_mottifiles = len(vuodet) + 1

    simYearSum = 0*np.ones(n_ditch_scens, dtype=int)
    poist0 = 0 # Ennen varsinaista simulointia hakattu puusto (oletus 0)
    
    standInfo = kuviot[kuviot['KUVIO']==stand]
    standInfo = standInfo.reset_index(drop=True)
    sfc = standInfo['kasvup'][0]
    puulaji = standInfo['puulaji'][0] # Huom. puulaji voi olla eri kuin yksittäisessä Susi-yhteensopivassa Motti-tiedostossa..
    
    sfc_spec = get_sfc_spec(sfc, puulaji) # kasvupaikkaluokan lisämääre (Mtkg, Ptkg I tai II), oletus II
 
    vol_ini_0 = np.nan
           
    # Assuming sfc between 2-5 (that is available in Susi)
    
    if sfc<2:
        sfc=2
        
    if sfc>5:
        sfc=5
    
    # If there are no cuttings:
    if len(pdata)==0:
        break_vol=None
        j = 0
        new_end_yr, motti_found = call_susi_help(kuviot, stand, j, end_yr, start_yr, wdata, sarka, break_vol, start_yr_ini, n_ditch_scens, np.nan, sfc_spec, ash_year=ash_year, vol_ini_0=vol_ini_0)
        
    
    if len(pdata)>0:
    
        j = 0
        
        while j < n_mottifiles:
                        
            break_vol = vol_before[j]
            
            if j==0:
                
                if vuodet[0]<5: # Jos ensimmäinen tapahtuma ennen v.5, vuosille 0-5 ei yleensä pysty muodostaa omaa järkevää motti-tiedostoa
                        poist0 = vol_before[0]-vol_after[0]
                        new_end_yr = start_yr
                        vol_ini_0 = vol_before[0]
                        j=j+1 # Hypätään suoraan seuraavaan Motti-tiedostoon
                        
                        if len(vuodet)>1: # Jos useampia kuin yksi tapahtuma...
                            if vuodet[1]<5: #.. ja toinenkin tapahtuma ennen v.5:
                                poist0 = poist0 + vol_before[1]-vol_after[1]
                                j=j+1 # Hypätään suoraan seuraavaan (kolmanteen) Motti-tiedostoon
                            
                else:   

                    new_end_yr, motti_found = call_susi_help(kuviot, stand, j, end_yr, start_yr, wdata, sarka, break_vol, start_yr_ini, n_ditch_scens, vol_after[j], sfc_spec, ash_year=ash_year, vol_ini_0=vol_ini_0)


                    with open(outfile, "a") as myfile:
                        myfile.write('\n' + 'j=0: ' + str(start_yr) + '...' + str(new_end_yr))  
                        
                    end_reached = check_end(n_ditch_scens, new_end_yr, end_yr)
                    
                    if end_reached == True:
                        break
                    
                    else:
                        
                        j=j+1
                        
                        if motti_found==True:
                            simYearSum = new_end_yr - start_yr + 1
                        else:
                            simYearSum=0
        
                        
                
            if j==n_mottifiles-1: # Jos tullaan viimeiseen Motti-tiedostoon...                  
                break_vol = None # Ei enää keskeytetä laskentaa, koska Motin mukaan tapahtumia ei tule
            else:
                break_vol = vol_before[j] # Muussa tapauksessa break_vol normaalisti
                
            
            if (j>0) & (j<n_mottifiles-1):
                
                start_yr_arr = start_yr + simYearSum # uusi aloitusvuosi
                
                new_end_yr, motti_found = call_susi_help(kuviot, stand, j, end_yr, start_yr_arr, wdata, sarka, break_vol, start_yr_ini, n_ditch_scens, vol_after[j], sfc_spec, ash_year=ash_year, vol_ini_0=vol_ini_0)    
                end_reached = check_end(n_ditch_scens, new_end_yr, end_yr)
                
                if end_reached == True: # Tsekataan, päästiinkö jo kokonaissimulointiajan loppuun
                    with open(outfile, "a") as myfile:
                        myfile.write('\n' + 'j=' + str(j) + ': ' + str(start_yr_arr) + '...' + str(new_end_yr))  
                    break
                else:
                    if motti_found==True:
                        
                        simYearSum = new_end_yr - start_yr + 1
                        
                        with open(outfile, "a") as myfile:
                            myfile.write('\n' + 'j=' + str(j) + ': ' + str(start_yr_arr) + '...' + str(new_end_yr) + ', ' + str(simYearSum) + ' years')  

                    else:

                        with open(outfile, "a") as myfile:
                            myfile.write('\n' + 'j=' + str(j) + ': ' + str(start_yr_arr) + '...' + str(new_end_yr) + ', ' + str(simYearSum) + ' years, Motti not found')  

                    j=j+1
            
            if j==n_mottifiles-1: # Jos tullaan viimeiseen Motti-tiedostoon...                  
                break_vol = None # Ei enää keskeytetä laskentaa, koska Motin mukaan tapahtumia ei tule
            else:
                break_vol = vol_before[j] # Muussa tapauksessa break_vol normaalisti
            
            if n_mottifiles > 2:
                if j == n_mottifiles-1: # Jos tullaan viimeiseen motti-tiedostoon...
                    
                    start_yr_arr = start_yr + simYearSum

                    new_end_yr, motti_found = call_susi_help(kuviot, stand, j, end_yr, start_yr_arr, wdata, sarka, break_vol, start_yr_ini, n_ditch_scens, vol_after[j], sfc_spec, ash_year=ash_year, vol_ini_0=vol_ini_0)    

                    with open(outfile, "a") as myfile:
                        myfile.write('\n' + 'j=' + str(j) + ': ' + str(start_yr_arr) + '...' + str(new_end_yr))  
                    
                        
                    end_reached = check_end(n_ditch_scens, new_end_yr, end_yr)
                    
                    if end_reached == True: # Tsekataan, päästiinkö jo kokonaissimulointiajan loppuun
                        with open(outfile, "a") as myfile:
                            myfile.write('\n End OK.')  
                        break
                    else:
                        with open(outfile, "a") as myfile:
                            myfile.write('\n Total simulation years not reached!')  
                        break

            else:   # Jos motti-tiedostoja 1 tai 2 kpl
                if j==1:
                    start_yr_arr = start_yr + simYearSum
                    
                    if j==n_mottifiles-1: # Jos viimeinen Motti-tiedosto...    
                        new_end_yr, motti_found = call_susi_help(kuviot, stand, j, end_yr, start_yr_arr, wdata, sarka, break_vol, start_yr_ini, n_ditch_scens, np.nan, sfc_spec, ash_year=ash_year, vol_ini_0=vol_ini_0)    

                    else:
                        new_end_yr, motti_found = call_susi_help(kuviot, stand, j, end_yr, start_yr_arr, wdata, sarka, break_vol, start_yr_ini, n_ditch_scens, vol_after[j], sfc_spec, ash_year=ash_year, vol_ini_0=vol_ini_0)    
                    
                    with open(outfile, "a") as myfile:
                        myfile.write('\n' + 'j=' + str(j) + ': ' + str(start_yr_arr) + '...' + str(new_end_yr))  
                                            
                    end_reached = check_end(n_ditch_scens, new_end_yr, end_yr)
                    
                    if end_reached == True: # Tsekataan, päästiinkö jo kokonaissimulointiajan loppuun
                        with open(outfile, "a") as myfile:
                            myfile.write('\n End OK.')  
                        break
                    else:
                        with open(outfile, "a") as myfile:
                            myfile.write('\n Total simulation years not reached!')  
                        break



def call_susi_help(kuviot, stand, j, end_yr, start_yr_arr, wdata, sarka, break_vol, start_yr_ini, n_ditch_scens, vol_aft, sfc_spec, maku=None, ash_year=None, vol_ini_0=None):

    susiPath, wpath, mottipath, folderName, outfile = get_paths()

    mottifile = {'path':mottipath,
                  'dominant':{1: str(stand) + '_n' + str(j) + '.xls'},
                  'subdominant':{0:'susi_motti_input_lyr_1.xlsx'}, # if key=0, file name can be anything
                  'under':{0:'susi_motti_input_lyr_2.xlsx'}}  # if key=0, file name can be anything
    
    susi_outfile = str(stand) + '_n' + str(j) + '.nc'
    
    standInfo = kuviot[kuviot['KUVIO']==stand]
    standInfo = standInfo.reset_index(drop=True)
    
    treesp = pd.read_excel(mottifile['path'] + '/' + mottifile['dominant'][1], sheet_name='Kertymät')['id Puulaji'][0] # Puulajin pitäisi olla Motti-tiedostossa oikein, vaikka puulaji vaihtuu
    
    try:
        df0 = pd.read_excel(mottifile['path'] + '/'+ mottifile['dominant'][1], sheet_name=0)
        motti_found = True
    except:
        with open(outfile, "a") as myfile:
            myfile.write('\nMotti file missing:', mottifile['path'] + '/'+ mottifile['dominant'][1])
        
        new_end_yr=np.nan*np.ones(n_ditch_scens)
        motti_found = False
        return new_end_yr, motti_found            
    
    if len(df0)>0:
        
        df = get_motti_silvi(mottifile['path'] + '/'+ mottifile['dominant'][1])
            
        ageSim = df['age'][0]
        if ageSim==0.:
            ageSim=0.001
        
        if (treesp==1) | (treesp==2) | (treesp==3) | (treesp==4): # Simulaatio toteutuu vain, jos pÃ¤Ã¤puulaji tiedossa
            
            sfc = standInfo['kasvup'][0]
            subgroup = standInfo['alaryh'][0]
            sarkaSim = sarka
            
            if (subgroup>=2) & (subgroup<=3) & (sfc<=5): # Simulaatio toteutuu vain, jos subgroup 2 tai 3, ja sfc<=5
                # Korvet
                if (subgroup==2) & (treesp==1):
                    site = 'pine_S'
                if (subgroup==2) & (treesp==2):
                    site = 'spruce_S'
                    
                # Rämeet
                if (subgroup==3) & (treesp==1) & (sfc<5):
                    site = 'pine_A'
                if (subgroup==3) & (treesp==1) & (sfc==5):
                    site = 'pine_S'
                if (subgroup==3) & (treesp==2) & (sfc<5):
                    site = 'spruce_A'
                if (subgroup==3) & (treesp==2) & (sfc==5):
                    site = 'spruce_S'
        
                if treesp>2:
                    site = 'birch_A'
                
                # Nevat ja letot -> onko 
                if (subgroup>=4) & (treesp==1):
                    site = 'pine_S'

                if (subgroup>=4) & (treesp==1):
                    site = 'pine_S'
                    
                new_end_yr = call_local_susi_silvi(wdata, mottifile, sarkaSim, ageSim, sfc, sfc_spec, site, break_vol, start_yr_arr, start_yr_ini, end_yr, vol_aft, susi_outfile, ash_year)
             
                return new_end_yr, motti_found
            
            else:
                with open(outfile, "a") as myfile:
                    myfile.write('\nStand ' + str(stand) +'j=0: skipped simulation because subgroup not 2 or 3, or sfc>5') 

        else:
            with open(outfile, "a") as myfile:
                myfile.write('\nError in stand number ' + str(stand)) 

    else:
        with open(outfile, "a") as myfile:
            myfile.write('\nMottifile empty, stand number ' + str(stand))     
        
        
        

    
    
    
    
    