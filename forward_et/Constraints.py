# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:51:06 2018

@author: gmgo
"""


def cal_ft(Tmean):
    import numpy as np
    # This function calculates the temperature constrain. It accounts for the reduction in the photosintet efficiency when plants grow.
    # Read further in Garcia et al under the methods section.  Check table 2 for eqaution implemented. and assumptions.
    Topt = 25
    ft = 1.1814*(1+np.exp(0.3*(Topt-10-Tmean)))**-1
    ft[np.isinf(ft)] = np.nan
    return ft


def cal_fg(NDVI, FPAR):
    import Functions as fn
    FIPAR = fn.calFIPAR(NDVI)
    fc = FIPAR
    fg = FPAR/fc
    fg[fg < 0] = 0
    fg[fg > 1] = 1
    return fg


def calVegConstrains_fm(fPAR):
    import os
    import numpy as np
    folder = (os.path.dirname(os.path.realpath(__file__)))
    fPARmax = np.load(folder+'\\Data\\fPAR_max_Andalucia.npy')
    # fPARmax=np.load(folder+'\\fPAR_max.npy')
    fm = fPAR/fPARmax
    return fm


def calEQ_all(NetRadiationSoil, psi, s, Jul, Eqs_all):
    import numpy as np
    As = NetRadiationSoil
    Epsilon = s/psi
    Eq = Epsilon*As/(Epsilon+1)
    Eq = Eq*0.0864  # Converts to m<mjm2<DAY
    Eq = Eq*0.408  # To get mm
    Eqs_all[Jul-1, :, :] = Eq
    return Eqs_all


def calSOILMoistureConstrain_f_Zhang(count, Fzang_all, fdrying_all):
    import numpy as np
    import os
    i = count
    folder = (os.path.dirname(os.path.realpath(__file__)))
    Precip = np.load(folder+'\\Data\\Precipitation_Andalucia.npy')
    RainyDay = np.load(folder+'\\Data\\Rainyday_Andalucia_1100m.npy')
    Precip = np.load(folder+'\\Data\\Precipitation_Andalucia.npy')
    Eq = np.load(folder+'\\Data\\EQ_All_Andalucia_new_.npy')
    if i < 15:
        SumPrec = np.nansum(Precip[0:i+1, :, :], axis=0)
        SumPrec[SumPrec < 0] = 0

        Eqs_sum = np.nansum(Eq[0:i+1, :, :], axis=0)
        Eqs_sum[Eqs_sum == 0] = 0

    else:
        SumPrec = np.nansum(Precip[i-15:i+1, :, :], axis=0)
        SumPrec[SumPrec < 0] = 0
        Eqs_sum = np.nansum(Eq[i-15:i+1, :, :], axis=0)
        Eqs_sum[Eqs_sum == 0] = 0
    MaxRainDay = np.nanmax(RainyDay[:i+1, :, :], axis=0)
    Delta_t = i-MaxRainDay
    alpha_f = 0.5
    MaxPos = np.nanargmax(RainyDay[:i+1, :, :], axis=0)
    Fzang = SumPrec/Eqs_sum
    Fzang[Fzang > 1] = 1
    Fzang[Fzang < 0] = 0
    Fzang_all[i, :, :] = Fzang
#    flp=np.zeros([834,1115])
#    for ii in range(0,834):
#        for iii in range(0,1115):
#            flp[ii,iii]=Fzang_all[(MaxPos[ii,iii]),ii,iii]
    flp = np.zeros([276, 478])
    for ii in range(0, 276):
        for iii in range(0, 478):
            flp[ii, iii] = Fzang_all[(MaxPos[ii, iii]), ii, iii]
    fdrying = flp*np.exp(-alpha_f*Delta_t)
    fdrying[fdrying < 0] = 0
    fdrying[fdrying > 1] = 1
    fdrying_all[:, :, i] = fdrying
    return Fzang_all, fdrying_all
