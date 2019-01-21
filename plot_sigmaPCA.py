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

  all_data['lastvec'] = {}
  ac_bus = all_data['autocorrbus']
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


  for i in np.arange(4*60*30,14*60*30,plfr):
    thisdata = data[:,(i-T):i]
    
    fn = np.zeros(1800)
    busdatat = thisdata[ac_bus,-Q:]
    for ss in range(1800):
      if ss==0:
        busdatas = thisdata[ac_bus,(-Q):]
      else:
        busdatas = thisdata[ac_bus,(-Q-ss):(-ss)]
      fn[ss] = np.dot(busdatat,busdatas)/Q
    fn = fn/fn[0]
    all_data['ac_fn'] = fn

    actualtime  = initialtime + dt.timedelta(seconds=i/30.)
    eig = doplot(log,all_data,strvar,actualtime,thisdata,plcount,eig)

    plcount += 1


        


def doplot(log, all_data, strvar, actualtime, data, plcount, eig):

  nbuses = all_data['nbuses']
  T = all_data['T']
  neigs = all_data['neigs']
  coords = all_data['coords']

  deletebuses = []
  drawbuses = []

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
    newvec = np.zeros(nbuses)
    newvec[drawbuses] = vec[i]
    if len(all_data['lastvec'])==neigs:
      if np.linalg.norm(newvec+all_data['lastvec'][i])<np.linalg.norm(newvec-all_data['lastvec'][i]):
        vec[i] = -vec[i]
    newvec[drawbuses] = vec[i]
    all_data['lastvec'][i] = newvec

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
    if thistime.second in [0,20,40] and i>3:
      dM.append(thistime)
      
  #plt.rcParams['font.family'] = 'Times New Roman'
  fig = plt.figure(tight_layout=True,figsize=(18,3))
  gs = matplotlib.gridspec.GridSpec(2, 7, width_ratios=[1,1,1,1,1,1,1], height_ratios=[1, 3] )

  #plotting raw data
  cmap2 = plt.get_cmap('rainbow')
  axdata = plt.subplot(gs[0:8:7])
  x,y,c = [130,130],[30,30],[-6,6]
  for i in range(len(drawbuses)):
    x.append(coords[drawbuses[i]][0])
    y.append(coords[drawbuses[i]][1])
    c.append(data[i,-1])
  cdata = axdata.scatter(x,y,15,c,alpha=1,cmap=cmap2,edgecolor='0.7')
  axdata.set_xlim(141.4,146.6)
  axdata.set_ylim(25,28.6)
  #axdata.text(146.55, 28.4, actualtime, fontsize=8, ha='right')
  axdata.set_title(r'$t$ = '+actualtime.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3], fontsize=9)
  axdata.tick_params(axis='both',which='both',bottom='off',left='off',labelbottom='off',labelleft='off')

  cbarticks = [-5,-2,0,2,5]
  cbar = fig.colorbar(cdata,ticks=[cbarticks],orientation='horizontal')
  cbar.ax.set_xticklabels(['%d'%cbt for cbt in cbarticks])
  cbar.ax.tick_params(labelsize=8)
  cbar.ax.set_title('Normalized ' + strvar,fontsize=8)

  
  #plotting all eigvals
  axallvals = plt.subplot(gs[1:9:7])
  axallvals.bar(1+np.arange(len(sortedeigvals)),sortedeigvals)
  axallvals.set_xlim([0,40])
  axallvals.tick_params(axis='both', which='major', labelsize=8, direction='out')
  axallvals.text(0.98,0.94,'size = '+str(len(eigvals)),size=8,ha='right',va='bottom',transform=axallvals.transAxes)
  axallvals.text(0.98,0.88,r'$T$' + ' = %d [s]'%(T/30.),size=8,ha='right',va='bottom',transform=axallvals.transAxes)
  axallvals.text(0.98,0.82,r'$S$' + ' = %d [s]'%(all_data['S']/30.),size=8,ha='right',va='bottom',transform=axallvals.transAxes)
  axallvals.text(1,1.02,'First 40 eigvals of '+strvar,fontsize=8,ha='right',va='bottom',transform=axallvals.transAxes)
  #axallvals.yaxis.get_major_formatter().set_powerlimits((-2, 2))    
  ymin,ymax = axallvals.get_ylim()
  ymax = ymax*1.05
  axallvals.set_ylim([ymin,ymax])
  yexp = math.floor(math.log10(ymax))
  axallvals.set_yticks([ymin,ymax/4.0,ymax/2.0,ymax*3.0/4,ymax])
  if ymax*10**(-yexp)>9.99:
    yexp += 1
  ymin,ymax = 0,ymax*10**(-yexp)
  axallvals.set_yticklabels(['0','%.2f'%(ymax/4.0),'%.2f'%(ymax/2.0),'%.2f'%(ymax*3.0/4.0),'%.2f'%ymax])
  axallvals.text(0,1.01,'1e%d'%(yexp),size=8,transform=axallvals.transAxes)
  plt.setp( axallvals.xaxis.get_majorticklabels(), ha='right' )

      
  # plotting first 4 eigvals/vects
  axval,axvec = {},{}
  axval[0] = plt.subplot(gs[2])
  axval[1] = plt.subplot(gs[3])
  axval[2] = plt.subplot(gs[4])
  axval[3] = plt.subplot(gs[5])

  axvec[0] = plt.subplot(gs[9])
  axvec[1] = plt.subplot(gs[10])
  axvec[2] = plt.subplot(gs[11])
  axvec[3] = plt.subplot(gs[12])

  strvec = ['1st','2nd','3rd','4th']
  latval = [r'$\lambda_1(t)$',r'$\lambda_2(t)$',r'$\lambda_3(t)$',r'$\lambda_4(t)$']
  latvec = [r'$\xi_1(t)$',r'$\xi_2(t)$',r'$\xi_3(t)$',r'$\xi_4(t)$']
  for i in axval.keys():
    axval[i].tick_params(axis='both', which='major', labelsize=6, direction='in')
    axvec[i].text(146.55, 28.3,latvec[i] + ': '+ strvec[i]+' eigvector', fontsize=8, ha='right')
    axvec[i].set_xlim(141.4,146.6)
    axvec[i].set_ylim(25,28.6)
    axvec[i].tick_params(axis='both',which='both',bottom='off',left='off',labelbottom='off',labelleft='off')

  xfmt = md.DateFormatter('%H:%M:%S')

  for i in [0,1,2,3]:
    axval[i].scatter(eig[neigs],eig[i],s=1)
    axval[i].xaxis.set_major_formatter(xfmt)
    axval[i].set_xlim([dlim[0],dlim[1]])
    #ylim = [0,1.1*max(fr_eig[i])/2,1.1*max(fr_eig[i])]   
    if len(eig[i])>60*30/all_data['plfr']:
      plotedeig = eig[i][-(60*30/all_data['plfr']):]
    else:
      plotedeig = eig[i]
    ylim = [0.9*min(plotedeig),1.1*max(plotedeig)]
    axval[i].set_ylim([ylim[0],ylim[-1]])
    yexp = math.floor(math.log10(ylim[-1]))
    axval[i].text(0,1.05,'1e%d'%(yexp),size=8,transform=axval[i].transAxes)
    axval[i].set_yticks([ylim[0],(ylim[0]+ylim[-1])/2.0,ylim[-1]])
    ylim = [ylim[0]*10**(-yexp),ylim[-1]*10**(-yexp)]
    axval[i].set_yticklabels(['%.2f'%ylim[0],'%.2f'%((ylim[0]+ylim[-1])/2.0),'%.2f'%ylim[-1]])
    axval[i].text(1,1.02,latval[i]+': '+strvec[i]+' eigval of '+strvar,fontsize=8,ha='right',va='bottom',transform=axval[i].transAxes)
    axval[i].text(1.01,0,r'$t$',fontsize=7,ha='left',va='center',transform=axval[i].transAxes)
    axval[i].set_xticks(dM)
    axval[i].set_xticks(dm, minor=True)
    axval[i].grid(which='minor', alpha=0.3)
    axval[i].grid(which='major', alpha=0.7)
    plt.setp( axval[i].xaxis.get_majorticklabels(), rotation=0, ha='right')

  cmap = plt.get_cmap('seismic')
  c,cax = {},{}
  x,y = [130,130],[30,30]
  for i in range(neigs):
    c[i] = [-1,1]

  for i in range(len(drawbuses)):
    c_x,c_y = coords[drawbuses[i]][0],coords[drawbuses[i]][1]
    x.append(c_x)
    y.append(c_y)
    for j in range(neigs):
      c[j].append(vec[j][i])
  for i in range(neigs):
    cax[i] = axvec[i].scatter(x,y,18,c[i],alpha=1,cmap=cmap,edgecolor='0.5',linewidth='0.2')


  #plotting autocorr

  axautocorr = plt.subplot(gs[6:14:7])
  axautocorr.text(1,1.02,'Autocorrel. of '+strvar,fontsize=8,ha='right',va='bottom',transform=axautocorr.transAxes)
  axautocorr.tick_params(axis='both', which='major', labelsize=8, direction='in')
  #axautocorr.plot(np.arange(300)/30.0,fr_all['autocorr']/fr_all['autocorr'][0])
  axautocorr.plot(np.arange(1800)/30.0,all_data['ac_fn'],linewidth=0.7)
  axautocorr.set_xlim([-0.5,10.5])          
  axautocorr.set_ylim([-1.5,1.5])
  #axautocorr.set_ylim([-1,1])    
  axautocorr.set_yticks([-1,0,1])
  axautocorr.set_yticklabels(['-1','0','1'])
  axautocorr.text(0.98,0.88,'bus %d'%all_data['autocorrbus'],size=8,ha='right',va='bottom',transform=axautocorr.transAxes)
  axautocorr.text(0.98,0.93,r'$T$' + ' = %d [s]'%(1800/30.),size=8,ha='right',va='bottom',transform=axautocorr.transAxes)
  #axautocorr.text(0.98,0.94,'Data normalized by mov.avg.',size=6,ha='right',va='bottom',transform=axautocorr.transAxes)
  axautocorr.text(0.98,0.04,r'$\Delta$'+' [s]',size=8,ha='right',va='bottom',transform=axautocorr.transAxes)
  axautocorr.grid(which='major', alpha=0.3)


  if strvar=='VoltageMag':
    axautocorr.text(1.01,-0.12,all_data['sysargv'][0]+'  '+all_data['videostamp'],size=6,ha='right',va='bottom',transform=axautocorr.transAxes)



  log.joint(all_data['foldername']+'/plot_' + strvar + '_%05d.png'%plcount + '\n')
  fig.savefig(all_data['foldername']+'/plot_' + strvar + '_%05d.png'%plcount)
  plt.clf()
  plt.cla()
  plt.close(fig)


  return eig



if __name__ == '__main__':
  if len(sys.argv) > 3 or len(sys.argv) < 2:
    print 'Usage: plot_sigmaPCA.py file.config [logfile]'
    exit(0)
  t0 = time.time()

  mylogfile = "plot_sigmaPCA.log"

  if len(sys.argv) == 3:
    mylogfile = sys.argv[2]

  log = Logger(mylogfile)

  all_data = read_configfile(log,sys.argv[1],sys.argv[0].split('/')[-1][:-3])
  all_data = read_nominalvoltsfile(log,all_data)
  all_data = read_coordfile(log,all_data)
  all_data['autocorrbus'] = 1
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


          
      


