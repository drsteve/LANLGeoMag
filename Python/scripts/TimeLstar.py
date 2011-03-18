# -*- coding: utf-8 -*-
from __future__ import division

import datetime
import LstarVersusPA
import spacepy.toolbox as tb
from pylab import *
import socket

now = datetime.datetime.now

pos = [-4,0,0]
date = datetime.datetime(2010, 12, 1)


results = {}

for i in range(9):
    results[i] = {}
    results[i]['Lsimple']=[]
    results[i]['time']=[]
    results[i]['Lstar']=[]


for qual in range(9):
    for n in range(20):
        print qual, n
        t0 = now()
        ans = LstarVersusPA.LstarVersusPA(pos, date, LstarQuality=qual, Bfield = 'Lgm_B_T89')
        t1 = now()
        results[qual]['time'].append((t1-t0).seconds + (t1-t0).microseconds/1000000)
        results[qual]['Lsimple'].append(ans[90]['Lsimple'])
        try:
            results[qual]['Lstar'].append(ans[90]['Lstar'])
        except KeyError:
            pass

colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g')
data = []
for i in range(9):
    for j, val in enumerate(results[i]['time']):
        plot(j, val, color=colors[i], marker='o', markersize=10)
    data.append(results[i]['time'])






figure()
boxplot(data, notch=1, positions=range(9))
ax = gca()
ax.set_xlabel('Quality number')
ax.set_ylabel('Run time')
ax.set_title(socket.gethostname() + ' LstarVersusPA Calcs ' +
             str(datetime.datetime.now().month) + '-' +
             str(datetime.datetime.now().day) +
             '-' + str(datetime.datetime.now().year))
#draw()
dstr = datetime.datetime.now().strftime('%d%b%Y')
savefig('TimeLstar_'+socket.gethostname()+'_'+dstr+'.png')
