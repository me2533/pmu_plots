import sys
import os
import numpy as np
import time
import math
import csv
import copy
from utils import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter
import matplotlib.dates as md
import datetime as dt
from matplotlib import cm

class ScalarFormatterForceFormat(ScalarFormatter):
  def _set_format(self):  # Override function that finds format to use.
    self.format = "%1.2f"  # Give format here                          



def analyzefiles(log, all_data, strvar):

  nbuses = all_data['nbuses']
  initialtime = dt.datetime.strptime('%s %s'%(all_data['datestr'],all_data['firsthourstr']),'%Y-%m-%d %H.%M')

  data = np.zeros((nbuses,27000))
  if strvar == 'Frequency':
    idxvar = 2
  elif strvar == 'VoltageAng':
    idxvar = 3
  elif strvar == 'VoltageMag':
    idxvar = 4


  T = all_data['T']
  S = all_data['S'] 
  Q = 30*60

  for bus in range(1,nbuses):
    filename = all_data['datafolder']+'/bus%03d.csv'%bus
    if os.path.isfile(filename)==0:
      continue
    
    log.joint('reading ' + filename + '\n')
    csvfile = open(filename,'rb')
    reader = csv.reader(csvfile)
    count = 0

    for row in reader:
      data[bus,int(round(30*float(row[0])))] = float(row[idxvar])
      count += 1

    if count!=27000:
      data[bus,:] = 0
      continue


  neigs = all_data['neigs']
  eig = {}
  for i in range(neigs+1):
    eig[i] = []

  rawdata = data
  data = 0*data

  for bus in range(nbuses):
    if sum(rawdata[bus,:]==0)==27000:
      continue

    for i in range(S,27000):
      mean = np.mean(rawdata[bus,(i-S):i])
      std = np.std(rawdata[bus,(i-S):i])
      data[bus,i] = (rawdata[bus,i]-mean)/std


  all_data['plfr'] = plfr = 10
  plidx = 0
  plcount = 0

  ac_fn = np.zeros(nbuses)
  firstacbus = all_data['firstacbus']

  for i in np.arange(4*60*30,14*60*30,plfr):
    thisdata = data[:,(i-T):i]

    for bus in range(nbuses):
      if sum(thisdata[bus,:]==0)==T:
        continue
    
      fn = np.zeros(1800)
      busdatat = thisdata[bus,-Q:]
      for ss in range(1800):
        if ss==0:
          busdatas = thisdata[bus,(-Q):]
        else:
          busdatas = thisdata[bus,(-Q-ss):(-ss)]
        fn[ss] = np.dot(busdatat,busdatas)/Q
      fn = fn/fn[0]

      ac_fn[bus] = np.mean(np.absolute(fn[30:]))


    if firstacbus[0]<0:
      sortedacinds = np.argsort(ac_fn)[::-1]
      firstacbus = sortedacinds[0:5]#sorted(sortedacinds[0:5])

    ac_firstbuses = np.zeros((5,1800))
    for ii in range(5):
      busdatat = thisdata[firstacbus[ii],-Q:]
      for ss in range(1800):
        if ss==0:
          busdatas = thisdata[firstacbus[ii],(-Q):]
        else:
          busdatas = thisdata[firstacbus[ii],(-Q-ss):(-ss)]
        ac_firstbuses[ii,ss] = np.dot(busdatat,busdatas)/Q
      ac_firstbuses[ii,:] = ac_firstbuses[ii,:]/ac_firstbuses[ii,0]


    actualtime  = initialtime + dt.timedelta(seconds=i/30.)
    eig = doplot(log,all_data,strvar,actualtime,thisdata,plcount,eig,ac_fn,firstacbus,ac_firstbuses)

    plcount += 1

        

          

def doplot(log, all_data, strvar, actualtime, data, plcount, eig, ac_fn, firstacbus, ac_firstbuses):

  nbuses = all_data['nbuses']
  T = all_data['T']
  neigs = all_data['neigs']
  coords = all_data['coords']

  deletebuses = []
  drawbuses = []

  # repeat last value when does not exist
  for bus in range(nbuses):
    if sum(data[bus,:]==0)==T:
      deletebuses.append(bus)
    else:
      drawbuses.append(bus)


  data = np.delete(data,deletebuses,0)
  corrmat = np.matmul(data,data.transpose())/T

  eigvals,eigvects = np.linalg.eigh(corrmat)
  sortedeigvals = np.sort(eigvals)[::-1]
  sortedeiginds = np.argsort(eigvals)[::-1]

  vec = {}
  for i in range(neigs):
    eig[i].append(sortedeigvals[i])
    vec[i] = eigvects[:,sortedeiginds[i]]

  dlim = [0,0]
  thistime = actualtime
  dlim[0] = thistime + dt.timedelta(seconds=-60)
  dlim[1] = thistime
  eig[neigs].append(thistime)

  dm,dM = [],[]
  for i in range(61):
    thistime = dlim[0] + dt.timedelta(seconds=i)
    if thistime.second in [0,10,20,30,40,50]:
      dm.append(thistime)
    if thistime.second in [0,20,40] and i>7:
      dM.append(thistime)
      
    
  fig = plt.figure(tight_layout=True,figsize=(20,3))
  gs = matplotlib.gridspec.GridSpec(1, 7, width_ratios=[1,2.3,1,1,1,1,1], height_ratios=[1] )

  axvals = plt.subplot(gs[0])
  axvals.tick_params(axis='both', which='major', labelsize=8, direction='in')

  xfmt = md.DateFormatter('%H:%M:%S')
  axvals.scatter(eig[neigs],eig[0],s=3,c='b')
  axvals.scatter(eig[neigs],eig[1],s=3,c='g')
  axvals.scatter(eig[neigs],eig[2],s=3,c='r')
  axvals.scatter(eig[neigs],eig[3],s=3,c='c')
  axvals.xaxis.set_major_formatter(xfmt)
  axvals.set_xlim([dlim[0],dlim[1]])
  ylim = [0,1.2*max(eig[0][-60:])]
  axvals.set_ylim([ylim[0],ylim[-1]])
  yexp = math.floor(math.log10(ylim[-1]))
  axvals.text(0,1.01,'1e%d'%(yexp),size=8,transform=axvals.transAxes)
  axvals.set_yticks([ylim[0],(ylim[0]+ylim[-1])/2.0,ylim[-1]])
  ylim = [ylim[0]*10**(-yexp),ylim[-1]*10**(-yexp)]
  axvals.set_yticklabels(['0','%.2f'%((ylim[0]+ylim[-1])/2.0),'%.2f'%ylim[-1]])
  axvals.text(1,1.01,'Largest 4 eigvals of '+strvar,fontsize=8,ha='right',va='bottom',transform=axvals.transAxes)
  axvals.text(0.99,0.94,r'$T$ = %d [s]'%(all_data['T']/(30)),size=8,ha='right',va='bottom',transform=axvals.transAxes)
  #axvals.set_title('Largest 4 eigvalues of '+strvar,fontsize=10)
  #axvals.text(0.99,0.93,actualtime,size=8,ha='right',va='bottom',transform=axvals.transAxes)
  axvals.set_xticks(dM)
  axvals.set_xticks(dm, minor=True)
  axvals.grid(which='minor', alpha=0.3)
  axvals.grid(which='major', alpha=0.7)
  plt.setp( axvals.xaxis.get_majorticklabels(), rotation=0, ha='right')


  deletebuses = []
  drawbuses = []

  # repeat last value when does not exist
  for bus in range(nbuses):
    if ac_fn[bus]==0:
      deletebuses.append(bus)
    else:
      drawbuses.append(bus)

  #autocorrelation indicator
  cmap2 = plt.get_cmap('rainbow')
  axdata = plt.subplot(gs[1])
  x,y,c,sz = [130,130],[30,30],[0,1],[0,0]
  for i in range(len(drawbuses)):
    x.append(coords[drawbuses[i]][0])
    y.append(coords[drawbuses[i]][1])
    c.append(ac_fn[drawbuses[i]]**.5)
    sz.append(10+ac_fn[drawbuses[i]]*100)
  #cdata = axdata.scatter(x,y,30,c,alpha=1,cmap=cm.plasma_r,edgecolor='0.4')
  cdata = axdata.scatter(x,y,s=sz,c=c,alpha=1,cmap=cm.plasma_r,edgecolor='1',linewidth='0')
  axdata.set_xlim(141.4,146.6)
  axdata.set_ylim(25,28.6)
  axdata.text(0.99,0.88,'size = '+str(len(drawbuses)),size=8,ha='right',va='bottom',transform=axdata.transAxes)
  #axdata.text(0.99,0.88,r'$R$ = 60 [s]',size=8,ha='right',va='bottom',transform=axdata.transAxes)
  axdata.text(0.99,0.94,r'$t$ = '+ actualtime.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3],size=8,ha='right',va='bottom',transform=axdata.transAxes)
  axdata.text(1,1.01,r'$\rho([\Delta_{min},\Delta_{max}];t)$: Autocorrelation residual map of '+strvar,fontsize=8,ha='right',va='bottom',transform=axdata.transAxes)
  #axdata.set_title(actualtime, fontsize=8)
  axdata.tick_params(axis='both',which='both',bottom='off',left='off',labelbottom='off',labelleft='off')

  cbarticks = [0,.1**.5,.2**.5,.3**.5,.4**.5,.5**.5,.8**.5,1]
  cbar = fig.colorbar(cdata,ticks=[cbarticks],orientation='vertical')
  cbar.ax.set_yticklabels(['%.1f'%(cbt**2) for cbt in cbarticks])
  cbar.ax.tick_params(labelsize=8)
  #cbar.ax.set_title('Autocorrelation indicator of ' + strvar,fontsize=8)

  clr = {firstacbus[0]:'xkcd:royal blue',firstacbus[1]:'xkcd:hunter green',firstacbus[2]:'xkcd:brick',firstacbus[3]:'xkcd:dull green',firstacbus[4]:'xkcd:grape'}
  mrk = {firstacbus[0]:'o',firstacbus[1]:'s',firstacbus[2]:'X',firstacbus[3]:'D',firstacbus[4]:'p'}
  for i in firstacbus:
    axdata.scatter(coords[i][0],coords[i][1],s=(70+ac_fn[i]*100),facecolors='none',edgecolors=clr[i],marker=mrk[i])


  axac = {}
  for i in range(5):
    axac[i] = plt.subplot(gs[2+i])
    axac[i].text(1,1.01,r'${\cal A}_i(\Delta;t)$: Autocorrelation of '+strvar,fontsize=8,ha='right',va='bottom',transform=axac[i].transAxes)
    axac[i].tick_params(axis='both', which='major', labelsize=8, direction='in')
    #axac[i].plot(np.arange(300)/30.0,fr_all['autocorr']/fr_all['autocorr'][0])
    axac[i].plot(np.arange(1800)/30.0,ac_firstbuses[i,:],color=clr[firstacbus[i]],linewidth=0.7)
    axac[i].set_xlim([-0.5,10.5])          
    axac[i].set_ylim([-1.3,1.3])
    #axac[i].set_ylim([-1,1])    
    axac[i].set_yticks([-1,0,1])
    axac[i].set_yticklabels(['-1','0','1'])
    axac[i].text(0.99,0.90,'bus i = %d'%firstacbus[i],size=8,ha='right',va='bottom',transform=axac[i].transAxes)
    axac[i].text(0.99,0.94,r'$T$ = %d [s]'%(60),size=8,ha='right',va='bottom',transform=axac[i].transAxes)
    #axac[i].text(0.99,0.94,'Data normalized by mov.avg.',size=6,ha='right',va='bottom',transform=axac[i].transAxes)
    axac[i].text(0.99,0.04,r'$\Delta$ [s]',size=8,ha='right',va='bottom',transform=axac[i].transAxes)
    axac[i].grid(which='major', alpha=0.3)
    axac[i].scatter(0.06,0.94,s=(70+ac_fn[firstacbus[i]]*60),facecolors='none',edgecolors=clr[firstacbus[i]],marker=mrk[firstacbus[i]],transform=axac[i].transAxes)
    

  if strvar=='VoltageMag':
    axac[4].text(1.01,-0.12,all_data['sysargv'][0]+'  '+all_data['videostamp'],size=6,ha='right',va='bottom',transform=axac[4].transAxes)
    axdata.text(0,-0.12,r'$S = 1[s]$, $\Delta_{min} = 1[s]$, $\Delta_{max} = 60[s]$',size=8,ha='left',va='bottom',transform=axdata.transAxes)



  log.joint(all_data['foldername']+'/plot_' + strvar + '_%05d.png'%plcount + '\n')
  fig.savefig(all_data['foldername']+'/plot_' + strvar + '_%05d.png'%plcount)
  plt.clf()
  plt.cla()
  plt.close(fig)


  return eig


if __name__ == '__main__':
  if len(sys.argv) > 3 or len(sys.argv) < 2:
    print 'Usage: plot_autocorr.py file.config [logfile]'
    exit(0)
  t0 = time.time()

  mylogfile = "plot_autocorr.log"

  if len(sys.argv) == 3:
    mylogfile = sys.argv[2]

  log = Logger(mylogfile)

  all_data = read_configfile(log,sys.argv[1],sys.argv[0].split('/')[-1][:-3])
  all_data = read_nominalvoltsfile(log,all_data)
  all_data = read_coordfile(log,all_data)
  all_data['sysargv'] = sys.argv

  if os.path.isdir(all_data['foldername'])==0:
    os.makedirs(all_data['foldername'])

  if all_data['variable2plot']=='all':
    all_data['videostamp'] = createstamp()
    log2 = Logger(all_data['foldername']+'/'+all_data['videostamp']+'.txt')
    log2.closelog()

    analyzefiles(log,all_data,'Frequency')
    analyzefiles(log,all_data,'VoltageAng')
    analyzefiles(log,all_data,'VoltageMag')


  else:
    if all_data['variable2plot']=='VoltageMag':
      all_data['videostamp'] = createstamp()
      log2 = Logger(all_data['foldername']+'/'+all_data['videostamp']+'.txt')
      log2.closelog()

    analyzefiles(log,all_data,all_data['variable2plot'])


  log.closelog()


          
      


