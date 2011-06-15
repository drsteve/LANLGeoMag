#!/usr/bin/env python

'''Script to verify the Lstar calculation of LanlGeoMag

Uses the centred dipole to get the expected value'''

from __future__ import division
import matplotlib.pyplot as plt
import itertools

import lgmpy
import spacepy.time as spt
import spacepy.datamodel as dm
import spacepy.toolbox as tb

nvals = 5
dates = spt.tickrange('2001-04-22T12:00:00', '2001-04-23T23:00:00', 1/24)
pos = dm.dmarray([-5, 0, 0], attrs={'sys': 'SM'})
alpha = 45
Lstar = [[]]*nvals
for date, qual in itertools.product(dates.UTC, range(nvals)):
    data = lgmpy.get_Lstar(pos, date, alpha=alpha, coord_system='SM', Bfield='Lgm_B_cdip', LstarThresh=20.0, extended_out=True, LstarQuality=qual)
    try:
        Lstar[qual].extend(data[alpha]['Lstar'])
    except TypeError:
        Lstar[qual].append(data[alpha]['Lstar'].tolist())

tb.pmm(Lstar)
plt.boxplot(Lstar)
plt.ylim([4.96, 5.02])
plt.show()
