# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 09:38:08 2020

"""

import sys

i = int(sys.argv[1])

from susi_silvi_run_hiilipolku_kuonan_hiib import call_local_susi_motti_silvi_list

call_local_susi_motti_silvi_list(i) # i=0-99
call_local_susi_motti_silvi_list(i+100) # i=100-199
call_local_susi_motti_silvi_list(i+200) # i=200-299
call_local_susi_motti_silvi_list(i+300) # i=300-399
call_local_susi_motti_silvi_list(i+400) # i=400-499
call_local_susi_motti_silvi_list(i+500) # i=500-599
call_local_susi_motti_silvi_list(i+600) # i=600-699
call_local_susi_motti_silvi_list(i+700) # i=700-799
