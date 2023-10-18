# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 09:38:08 2020

"""

import sys
import os
import pandas as pd

i = int(sys.argv[1])

from susi_silvi_run_hiilipolku_ import call_local_susi_motti_silvi_list
#from susi_silvi_run_hiilipolku_halvan_baua import call_local_susi_motti_silvi_list

def get_paths():
    
    susiPath = r'/scratch/project_2002470/SUSI_HIILIPOLKU/susi/' # susi python files path
    outPath=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/Neuvonta_lp4/Halvanjoki_HII_B/' # outputs path
    wpath = r'/scratch/project_2002470/SUSI_HIILIPOLKU/inputs/' # weather data path
    #mottifolder = r'/scratch/project_2002470/SUSI_HIILIPOLKU/Halvanjoki_HII_B/'
    mottipath = r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan_lp4/Halvanjoki/Halvanjoki_HII_B/'
    outfile = r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/Neuvonta_lp4/Halvanjoki_HII_B/Halvanjoki_HII_B.txt'    
    outfol=r'/scratch/project_2002470/SUSI_HIILIPOLKU_outputs/Neuvonta_lp4/' # outputs path
    if not os.path.exists(outfol):
        os.mkdir(outfol)    
    if not os.path.exists(outPath):
        os.mkdir(outPath)    
    return susiPath, wpath, mottipath, outPath, outfile

def get_motti_data():
        
    kuviot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/Halvanjoki/Halvanjoki_HII_B/motti/Halvanjoki_HII_B_kuviot.csv', encoding='latin1', sep=';')
    puustot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/Halvanjoki/Halvanjoki_HII_B/motti/Halvanjoki_HII_B_puustot.csv', encoding='latin1', sep=';')
    poistumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/Halvanjoki/Halvanjoki_HII_B/motti/Halvanjoki_HII_B_poistumat.csv', encoding='latin1', sep=';')
    tapahtumat = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/Halvanjoki/Halvanjoki_HII_B/motti/Halvanjoki_HII_B_tapahtumat.csv', encoding='latin1', sep=';', index_col=False)
    tulot = pd.read_csv(r'/scratch/project_2002470/HIILIPOLKU_data/Neuvontaan/Halvanjoki/Halvanjoki_HII_B/motti/Halvanjoki_HII_B_kasvut.csv', encoding='latin1', sep=';')
    
    return kuviot, puustot, poistumat, tapahtumat, tulot

def get_inputs():
        
    wdata='halvansuo_weathercont.csv' #'halvansuo_weather.csv' #'puruvesi_saa_19810101-2020062.csv' #'halvansuo_weather.csv' #'kuonanjoki_weather.csv'# # # Weather data file name
    sarka = 40. # Ditch spacing
    n_ditch_scens = 4 # Number of ditch scenarios 
    start_yr_ini = 1991 # Initial starting year
    totSimYears=50 # Total lenght of simulations, in years
    
    return wdata, sarka, n_ditch_scens, start_yr_ini, totSimYears

kuviot, puustot, poistumat, tapahtumat, tulot = get_motti_data()    
wdata, sarka, n_ditch_scens, start_yr_ini, totSimYears = get_inputs() 
susiPath, wpath, mottipath, folderName, outfile = get_paths()
    
call_local_susi_motti_silvi_list(i,kuviot, puustot, poistumat, tapahtumat, tulot,wdata, sarka, n_ditch_scens, start_yr_ini, totSimYears,susiPath, wpath, mottipath, folderName, outfile) # i=0-99
#call_local_susi_motti_silvi_list(i+100,kuviot, puustot, poistumat, tapahtumat, tulot,wdata, sarka, n_ditch_scens, start_yr_ini, totSimYears,susiPath, wpath, mottipath, folderName, outfile) # i=100-199
#call_local_susi_motti_silvi_list(i+200,kuviot, puustot, poistumat, tapahtumat, tulot,wdata, sarka, n_ditch_scens, start_yr_ini, totSimYears,susiPath, wpath, mottipath, folderName, outfile) # i=200-299
#call_local_susi_motti_silvi_list(i+300,kuviot, puustot, poistumat, tapahtumat, tulot,wdata, sarka, n_ditch_scens, start_yr_ini, totSimYears,susiPath, wpath, mottipath, folderName, outfile) # i=300-399

