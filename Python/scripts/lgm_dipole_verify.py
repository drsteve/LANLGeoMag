#!/usr/bin/env python

'''Script to verify the Lstar calculation of LanlGeoMag

Uses the centred dipole to get the expected value'''

from __future__ import division
import matplotlib.pyplot as plt
import itertools

import lgmpy, copy
import spacepy.time as spt
import spacepy.datamodel as dm
import spacepy.toolbox as tb

#make figure
fig = plt.figure(figsize=(17, 21))

#set up test runs
nvals, rvals = 5, [2, 3, 4, 5, 6]
dates = spt.tickrange('2001-04-22T12:00:00', '2001-04-23T23:00:00', 1/24)
alpha = 45

#loop over radial positions, dates and quality flags
for rval in rvals:
    Lstar = [[]]*nvals
    for date, qual in itertools.product(dates.UTC, range(nvals)):
        pos = dm.dmarray([-1*rval, 0, 0], attrs={'sys': 'SM'})
        data = lgmpy.get_Lstar(pos, date, alpha=alpha, coord_system='SM', Bfield='Lgm_B_cdip', LstarThresh=20.0, extended_out=True, LstarQuality=qual)
        try:
            Lstar[qual].extend(data[alpha]['Lstar'])
        except TypeError:
            Lstar[qual].append(data[alpha]['Lstar'].tolist())
    print('Did [-%d,0,0] for all qualities' % rval)
    #make plots
    fstr = '%d1%d' % (len(rvals), rval-rvals[0]+1)
    ax = fig.add_subplot(fstr)
    ax.boxplot(Lstar)
    #ax = plt.gca()
    ax.set_title('LGM - Centred dipole [-%d, 0, 0]$_{SM}$; PA=%d$^{o}$' % (rval, alpha))
    ax.set_ylabel('L* (LanlGeoMag)')
    ax.set_xlabel('Quality Flag')
    ax.set_xticklabels([str(n) for n in range(5)])
    
    tb.savepickle('lgm_cdip_lstar%d_alpha%d.pkl' % (rval, alpha), {'Lstar': Lstar})

fig.savefig('lgm_cdip_verify%d_zoom.png' % alpha, dpi=300)
    ##plt.ylim([4.96, 5.02])
    ##plt.savefig('lgm_cdip_verify%d.png' % alpha, dpi=300)
