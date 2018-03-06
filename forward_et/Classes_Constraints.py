# -*- coding: utf-8 -*-
"""
Created on Thu Feb 01 12:33:32 2018

@author: gmgo
"""
#Calculate vegetation constraints
import numpy as np


class VegConstraints(object):
    Tmean=0

    def __init__(self, Dictionary):
        self.Tmax = Dictionary['DASEMON_Tmax'][:]
        self.Tmin = Dictionary['DASEMON_Tmin'][:]
        self.NDVI = Dictionary['NDVI_16'][:]
        self.FPAR = Dictionary['fPAR'][:]

    def ft(self):
        Tmean=(self.Tmax + self.Tmin)/2
        #This function calculates the temperature constrain. It accounts for the reduction in the photosintet efficiency when plants grow.
        #Read further in Garcia et al under the methods section.  Check table 2 for eqaution implemented. and assumptions.
        Topt=25
        ft=1.1814*(1+np.exp(0.3*(Topt-10-Tmean)))**-1
#        ft[np.isinf(ft)]=np.nan
        return ft
    def fg(self):
        import Functions as fn
        FIPAR=fn.calFIPAR(self.NDVI )
        fc=FIPAR
        fg=self.FPAR/fc
        fg[fg<0]=0
        fg[fg>1]=1
        return fg
    def fm(self):
        fPARmax=np.load('fPAR_max.npy')
        fm=(self.FPAR)/fPARmax
        return fm 

