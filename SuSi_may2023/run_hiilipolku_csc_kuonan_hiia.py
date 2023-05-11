# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 09:38:08 2020

"""

import sys

i = int(sys.argv[1])

from susi_silvi_run_hiilipolku_kuonan_hiia import call_local_susi_motti_silvi_list

call_local_susi_motti_silvi_list(i) # i=0-99
call_local_susi_motti_silvi_list(i+100) # i=100-199
call_local_susi_motti_silvi_list(i+200) # i=200-299