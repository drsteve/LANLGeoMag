# -*- coding: utf-8 -*-
from __future__ import division

import datetime
import LstarVersusPA
import spacepy.toolbox as tb
from pylab import *
import socket
import scipy.stats

now = datetime.datetime.now

pos = [-4,0,0]
date = datetime.datetime(2010, 12, 1)

print('*'*50)
print("Something is wrong here with the plotting, don't run me")
print('*'*50)

results = {}

for i in range(9*2):
    results[i] = {}
    results[i]['Lgm_B_T89'] = {}
    results[i]['Lgm_B_OP77'] = {}

    results[i]['Lgm_B_T89']['LMcIlwain']=[]
    results[i]['Lgm_B_T89']['time']=[]
    results[i]['Lgm_B_T89']['Lstar']=[]

    results[i]['Lgm_B_OP77']['LMcIlwain']=[]
    results[i]['Lgm_B_OP77']['time']=[]
    results[i]['Lgm_B_OP77']['Lstar']=[]

nQual = 9
nTrials = 20

for qual in range(nQual):
    for n in range(nTrials): # get 20 samples to compare
        print qual, n
        t0 = now()
        ans = LstarVersusPA.LstarVersusPA(pos, date, LstarQuality=qual, Bfield = 'Lgm_B_OP77')
        t1 = now()
        t00 = now()
        ans = LstarVersusPA.LstarVersusPA(pos, date, LstarQuality=qual, Bfield = 'Lgm_B_T89')
        t11 = now()
        results[qual]['Lgm_B_T89']['time'].append((t1-t0).seconds + (t1-t0).microseconds/1000000)
        results[qual]['Lgm_B_T89']['LMcIlwain'].append(ans[90]['LMcIlwain'])
        try:
            results[qual]['Lgm_B_T89']['Lstar'].append(ans[90]['Lstar'])
        except KeyError:
            pass
        results[qual]['Lgm_B_OP77']['time'].append((t11-t00).seconds + (t11-t00).microseconds/1000000)
        results[qual]['Lgm_B_OP77']['LMcIlwain'].append(ans[90]['LMcIlwain'])
        try:
            results[qual]['Lgm_B_OP77']['Lstar'].append(ans[90]['Lstar'])
        except KeyError:
            pass

dataT89 = []
for i in range(nQual):
    dataT89.append(results[i]['Lgm_B_T89']['time'])
dataT89 = array(dataT89)
dataOP77 = []
for i in range(nQual):
    dataOP77.append(results[i]['Lgm_B_OP77']['time'])
dataOP77  = array(dataOP77)

data = hstack([dataT89.transpose(), dataOP77.transpose()])
data = hstack([data[0:,0::2], data[0:,1::2]]) # T89 and OP77 alternate now

# compute statistics on if the medians are different in the two runs
stats = []
for val in range(nQual):
    # for each pair of data do the test
    stats.append(scipy.stats.mannwhitneyu(data[:,val*2], data[:,val*2+1])[1]*2) # get p value *2 for 2 sided
for i, val in enumerate(stats):
    if val < 0.05:
        stats[i] = 'Different'
    else:
        stats[i] = 'Same'

figure()
boxplot(data, notch=True, positions=range(nQual*2))
ax = gca()
ax.set_xlabel('Quality number')
ax.set_ylabel('Run time')
ax.set_title(socket.gethostname() + ' LstarVersusPA Calcs ' +
             str(datetime.datetime.now().month) + '-' +
             str(datetime.datetime.now().day) +
             '-' + str(datetime.datetime.now().year))
ax.set_xticks(range(nQual*2))
labels = [ [str(val)+'-T89\n' + stats[val], str(val)+'-OP77'] for val in  range(nQual)]
labels = list(flatten(labels))
ax.set_xticklabels(labels)

#draw()
dstr = datetime.datetime.now().strftime('%d%b%Y')
savefig('TimeLstar_T89_OP77'+socket.gethostname()+'_'+dstr+'.png')
