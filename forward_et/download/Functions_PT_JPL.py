# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:09:27 2018

@author: gmgo
"""


def CheckDataConsistency(fileName):
    import os
    folder = os.walk(fileName)
    return True


def JulDay_JulDay8fromDate(date):
    from jdcal import gcal2jd

    year = str(date.year)
    month = date.month
    if month < 10:
        monthstr = '0'+str(month)
    else:
        monthstr = str(month)

    day = (date.day)
    if day < 10:
        daystr = '0'+str(day)
    else:
        daystr = str(day)
    day_full = year+'-'+monthstr+'-'+daystr

    Jul = gcal2jd(year, 1, 1)
    JulAct = gcal2jd(year, month, day)
    JulianDay = (JulAct[1]+1-Jul[1])
    JulDay8 = ('001', '009', '017', '025', '033', '041', '049', '057', '065', '073', '081',
               '089', '097', '105', '113', '121', '129', '137', '145', '153', '161', '169',
               '177', '185', '193', '201', '209', '217', '225', '233', '241', '249', '257',
               '265', '273', '281', '289', '297', '305', '313', '321', '329', '337', '345',
               '353', '361')
    JulDayIni = ('001', '009', '017', '025', '033', '041', '049', '057', '065', '073', '081',
                 '089', '097', '105', '113', '121', '129', '137', '145', '153', '161', '169',
                 '177', '185', '193', '201', '209', '217', '225', '233', '241', '249', '257',
                 '265', '273', '281', '289', '297', '305', '313', '321', '329', '337', '345',
                 '353')
    JulDayEnd = ('009', '017', '025', '033', '041', '049', '057', '065', '073', '081', '089',
                 '097', '105', '113', '121', '129', '137', '145', '153', '161', '169', '177',
                 '185', '193', '201', '209', '217', '225', '233', '241', '249', '257', '265',
                 '273', '281', '289', '297', '305', '313', '321', '329', '337', '345', '353',
                 '361')
    for i in range(0, len(JulDay8)):
        if JulianDay == int(JulDay8[i]):
            return JulianDay, JulDay8[i]
        if (JulianDay > int(JulDayIni[i])) and (JulianDay < int(JulDayEnd[i])):
            return JulianDay, JulDay8[i]
        if JulianDay >= 361:
            i = 45
            return int(JulianDay), int(JulDay8[i])


def JulDay_JulDay16fromDate(date):
    from jdcal import gcal2jd

    year = str(date.year)
    month = date.month
    if month < 10:
        monthstr = '0'+str(month)
    else:
        monthstr = str(month)
    day = (date.day)
    if day < 10:
        daystr = '0'+str(day)
    else:
        daystr = str(day)
    day_full = year+'-'+monthstr+'-'+daystr

    Jul = gcal2jd(year, 1, 1)
    JulAct = gcal2jd(year, month, day)
    JulianDay = (JulAct[1]+1-Jul[1])
    JulDay16 = ('001', '017', '033', '049', '065', '081',
                '097', '113', '129', '145', '161', '177',
                '193', '209', '225', '241', '257', '273',
                '289', '305', '321', '337', '353')
    JulDayIni = ('001', '017', '033', '049', '065', '081',
                 '097', '113', '129', '145', '161', '177',
                 '193', '209', '225', '241', '257', '273',
                 '289', '305', '321', '337', '353')
    JulDayEnd = ('017', '033', '049', '065', '081', '097',
                 '113', '129', '145', '161', '177', '193',
                 '209', '225', '241', '257', '273', '289',
                 '305', '321', '337', '366')
    for i in range(0, len(JulDay16)):
        if JulianDay == int(JulDay16[i]):
            return JulDay16[i]
        if (JulianDay > int(JulDayIni[i])) and (JulianDay < int(JulDayEnd[i])):
            return JulDay16[i]
        if JulianDay >= 353:
            i = 22
            return int(JulDay16[i])
