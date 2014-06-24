#!/usr/bin/env python

'''Script to verify the Lstar calculation of LanlGeoMag

Uses the T89 model to check spread and convergence of L* calcs'''

from __future__ import division
import matplotlib.pyplot as plt
import itertools
import datetime as dt
import lgmpy, copy
import spacepy.time as spt
import spacepy.datamodel as dm
import spacepy.toolbox as tb

#make figure
fig = plt.figure(figsize=(17, 21))

#set up test runs
nvals, rvals = 5, [3,4,5,6]
dates = spt.tickrange('2001-04-22T12:00:00', '2001-04-22T12:01:00', dt.timedelta(seconds=5))
alpha = 45
print('T89: alpha=%s, date range [%s, %s]' % (alpha, dates.UTC[0], dates.UTC[-1]))
#loop over radial positions, dates and quality flags
for rval in rvals:
    print('Radial Distance: %s' % rval)
    Lstar = [[]]*nvals
    for date, qual in itertools.product(dates.UTC, range(nvals)):
        print('%s, Quality=%d' % (date, qual))
        pos = dm.dmarray([-1*rval, 0, 0], attrs={'sys': 'GSM'})
        data = lgmpy.get_Lstar(pos, date, alpha=alpha, coord_system='GSM', Bfield='Lgm_B_T89c', LstarThresh=20.0, extended_out=True, LstarQuality=qual)
        try:
            Lstar[qual].extend(data[alpha]['Lstar'])
        except TypeError:
            Lstar[qual].append(data[alpha]['Lstar'].tolist())

    #make plots
    fstr = '%d1%d' % (len(rvals), rval-rvals[0]+1)
    ax = fig.add_subplot(fstr)
    ax.boxplot(Lstar)
    ax.set_title('T89 [-%d, 0, 0]$_{GSM}$; PA=%d$^{o}$' % (rval, alpha))
    ax.set_ylabel('L* (LanlGeoMag)')
    ax.set_xlabel('Quality Flag')
    ax.set_xticklabels([str(n) for n in range(5)])
    
fig.savefig('lgm_T89_%d_zoom_5sec.png' % alpha, dpi=300)
