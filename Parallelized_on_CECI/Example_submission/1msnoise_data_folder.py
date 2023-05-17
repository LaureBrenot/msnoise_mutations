# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 14:20:21 2022

@author: breno
"""
import matplotlib
from msnoise import api
from msnoise import s000installer
import os

db = api.connect() #connect to database

#api.reset_jobs(db, 'CC', alljobs=True)
#api.reset_jobs(db, 'MWCS', alljobs=True)
#api.reset_jobs(db, 'DTT', alljobs=True)

volcano = "Pavlof"
data_location = "/scratch/ulb/gtime/lbrenot/Alaska_volcanoes/"+volcano+"/data_"+volcano
#data_location ="/CECI/trsf/ulb/gtime/lbrenot/Askja_trsf/"
#data_location = "/scratch/ulb/gtime/lbrenot/Aso/data/"

api.update_config(db, 'data_folder', data_location) #set data location

