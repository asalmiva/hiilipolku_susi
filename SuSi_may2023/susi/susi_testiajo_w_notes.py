# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 12:13:58 2022
Created on Tue Mar  1 19:52:35 2022

@author: 03081194
@author: alauren
otettu githubista master branchista päivitetty 17.8.2022,c815bf5, joudutiin muokkaa kyllä että toimi omalla koneella

"""

import numpy as np
import pandas as pd
import datetime
from susi_utils import  get_motti, read_FMI_weather, get_mese_input
from susi_para import get_susi_para
# from susi_main import run_susi
import susi_io
from susi_main import Susi
import os

#***************** local call for SUSI*****************************************************
susiPath = r'C:/Users/03081194/OneDrive - Valtion/SUSI/susi_2022-master/susi/' # susi python files path
folderName=r'C:/Users/03081194/OneDrive - Valtion/SUSI/output_test/' # outputs path lisaa tama args inputtina
wpath = r'C:/Users/03081194/OneDrive - Valtion/SUSI/susi_2022-master/inputs/' # weather data path lisaa tama args inputtina
wdata='parkano_weather.csv' #tarkista miten säädatan pitää olla, pitääkö olla kaikki päivät jne

mottifile = {'path':r'C:/Users/03081194/OneDrive - Valtion/SUSI/susi_2022-master/inputs/', #tassa ok olla xslx mutta olis hyv muuttaa myös csv:ksi jotta saadaan haettua hakkuita, ja muita site ja ikä tietoja
              'dominant':{1: 'susi_motti_input_lyr_0.xlsx'},#taa perusmotti input file, palastelut pitää olla tehtynä eka eli jos tulee hakkkuu
              'subdominant':{0:'susi_motti_input_lyr_1.xlsx'}, # if key=0, file name can be anything
              'under':{0:'susi_motti_input_lyr_2.xlsx'}}  # if key=0, file name can be anything

filename = 'test4.nc'#argumentti inputtina, pitaa huomioida et tekee ensin yhen nimisen ja lopuks muuttaa, jos moniajoa niin pitää hoitaa

start_date = datetime.datetime(2000,1,1) #argumenttina ehkä kans
end_date=datetime.datetime(2005,12,31) #argumenttina ehkä kans
start_yr = start_date.year 
end_yr = end_date.year
yrs = (end_date - start_date).days/365.25

sarkaSim = 40.                                  # strip width, ie distance between ditches, m, 
n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip

ageSim = {'dominant': 40.*np.ones(n),       #tahan vois ottaa motista iän
          'subdominant': 0*np.ones(n),
          'under': 0*np.ones(n)}                                          # age of the stand in each node

sfc =  np.ones(n, dtype=int)*3              #tähän kans se tulis motista tia muualta input sfc, leenalta koodit                            # site fertility class

site = 'develop_scens'    #susi_para.py tiedostossa soil and stand parameters (spara) dictionary elementti, johon voi muutttaa esim ojien syvyyden muutokset, nyt ei muutosta, nyt männylle vaan, mutta Leenalta tulee Kuuselle ja koivulle

forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
           
wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=None,  
                                                                          ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          susiPath=None,
                                                                          n=n)
        
susi_test = Susi()    #susi_main.py, täytyy miettiä miten lannoitusvuo saadan input argumenttina sparan kautta, ja miten hakkuuvuosi ja hakkuun tapa annetaan (to_ba :n kautta ehkä), stand.py:ssä periaatteessa , miten motti hoitaa hakkuut???
#   tuleeko --- Organic matter decomposition and nutrient release--------------- lakettua vaikka stand.pyssä tyhjää
#hakkuun tekeminen stand.pyssä? jatkuvapeittenen eri tavalla pitäisi selvittää toimiiko ja onnistuuko niin ettei pilkota mottitiedotsotja? ja jos pilkottaisiinkin niin miten niin että päivittyy spara/uudlel atiedostolla...
#joka tapauksessa output samaan yhteen fileen, olis hyvä useista tiedostoista, toissalta netdcf
#output/kuvio/reseptiseknaario/n1,n2,n3 hakkuiden mukaan -> yhteenveto



susi_test.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                        mottifile=mottifile, peat= 'other', photosite='All data', 
                        folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath)          

             
os.rename(folderName + 'susi.nc', folderName + filename)
