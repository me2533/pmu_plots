#!/usr/bin/python

import sys
import math
import cmath
import csv
import time
import datetime as dt


'''
LOGGING CLASS
'''
class Logger:
  def __init__(self, logfilename):
    self.logfile = open(logfilename,"w")
    self.screen = 1
    self.log = 1
    localtime= time.asctime(time.localtime(time.time()))
    self.logfile.write("starting log: %s\n" % localtime)
    from socket import gethostname; self.logfile.write("running on: " + gethostname() + "\n\n")

  def closelog(self):
    localtime= time.asctime(time.localtime(time.time()))
    self.logfile.write("\nclosing log: %s\n" % localtime)
    self.logfile.close()

  def joint(self, mystring):
    if self.log:
      self.logfile.write(mystring)
      self.logfile.flush()
    if self.screen:
      print(mystring),

  def screen_on(self):
    self.screen = 1

  def screen_off(self):
    self.screen = 0

  def log_on(self):
    self.log = 1

  def log_off(self):
    self.log = 0

  def stateandquit(self, stuff):
    self.log = self.screen = 1
    self.joint(stuff + "\n")
    self.closelog()
    sys.exit("\nquitting\n")    



'''
USEFUL FUNCTIONS
'''
def breakexit(foo):
    stuff = raw_input("("+foo+") break> ")
    if stuff == 'x' or stuff == 'q':
        sys.exit("bye")


def createstamp():
    date = dt.datetime.now()
    datestamp = '%4d%02d%02d_'%(date.year,date.month,date.day)
    timestamp = '%02d%02d%02d'%(date.hour,date.minute,date.second)

    return datestamp+timestamp


'''
FILE READERS
'''
def read_configfile(log,filename,pyfunction):
  log.joint('reading config file ' + filename + '\n')

  try:
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
  except IOError as (errno, strerror):
    log.stateandquit("cannot open file "+ filename)
    sys.exit("failure")

  all_data = {}

  coordsfilename = 'coords.csv'
  nominalvoltsfilename = 'nominalvolts.csv'
  T = 2*30*60
  S = 30
  variable2plot = 'Frequency'
  datafolder = 'data'
  filterfr = 4
  
  alpha = 0.05

  firstacbus = [-1, -1, -1, -1, -1] 

  linenum = 0
  while linenum < len(lines):
    thisline = lines[linenum]
    thisline = thisline[0:thisline.find('#')]
    thisline = thisline.split()
    if len(thisline) > 0:
      if thisline[0] == 'END':
        break
        
      elif thisline[0] == 'coordsfilename':
        coordsfilename = thisline[1]
        
      elif thisline[0] == 'nominalvoltsfilename':
        nominalvoltsfilename = thisline[1]

      elif thisline[0] == 'T':
        if len(thisline)==2:
          T = int(thisline[1])
        elif len(thisline)==3:
          if thisline[2]=='sec':
            T = 30*int(thisline[1])
          elif thisline[2]=='min':
            T = 30*60*int(thisline[1])

      elif thisline[0] == 'S':
        if len(thisline)==2:
          S = int(thisline[1])
        elif len(thisline)==3:
          if thisline[2]=='sec':
            S = 30*int(thisline[1])
          elif thisline[2]=='min':
            S = 30*60*int(thisline[1])

      elif thisline[0] == 'variable2plot':
        variable2plot = thisline[1]

      elif thisline[0] == 'alpha':
        alpha = float(thisline[1])

      elif thisline[0] == 'filterfr':
        filterfr = int(thisline[1])

      elif thisline[0] == 'datafolder':
        datafolder = thisline[1]

      elif thisline[0] == 'firstacbus':
        if len(thisline)!=6:
          log.stateandquit("wrong line size: firstacbus")
        for k in range(1,6):
          firstacbus[k-1] = int(thisline[k])

      elif thisline[0] == 'DATESTRING':
        log.joint(' '.join(thisline) + '\n')
        all_data['datestr'] = datestr = thisline[1]
        firsthour  = int(thisline[3])
        firstmin = int(thisline[5])
        secondhour  = int(thisline[7])
        secondmin = int(thisline[9])
        all_data['firsthourstr'] = '%02d.%02d'%(firsthour,firstmin)

        all_data['timestr'] = '%02d.%02d_to_%02d.%02d'%(firsthour,firstmin,secondhour,secondmin)
        all_data['endtime'] = datestr + ' %02d:%02d:00.000'%(secondhour,secondmin)
        all_data['logfilename'] = datestr + '_%02d.%02d'%(firsthour,firstmin) + '_to_' + '%02d.%02d'%(secondhour,secondmin) + '.log'
        all_data['foldername'] = datestr + '_%02d%02d'%(firsthour,firstmin) + '_' + pyfunction

      else:
        print 'illegal input', thisline[0]
        sys.exit('bye')

    linenum += 1

  all_data['nbuses'] = 240
  all_data['neigs'] = 4
  all_data['alpha'] = alpha
  all_data['coordsfilename'] = coordsfilename
  all_data['nominalvoltsfilename'] = nominalvoltsfilename
  all_data['T'] = T
  all_data['S'] = S
  all_data['variable2plot'] = variable2plot
  all_data['datafolder'] = datafolder
  all_data['filterfr'] = filterfr
  all_data['firstacbus'] = firstacbus
  
  return all_data



def read_coordfile(log,all_data):
  filename = all_data['coordsfilename']
  log.joint('reading coord file ' + filename + '\n')

  try:
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
  except IOError as (errno, strerror):
    log.stateandquit("cannot open file "+ filename)
    sys.exit("failure")

  coords = {}
  linenum = 0

  max_x,max_y = -1e10,-1e10
  min_x,min_y = 1e10,1e10

  while linenum < len(lines):
    thisline = lines[linenum]
    thisline = thisline[0:thisline.find('#')]
    thisline = thisline.split()
    if len(thisline) == 3:
      x,y = float(thisline[0]), float(thisline[1])
      coords[int(thisline[2])] = [ x, y ]
      max_x = max(x,max_x)
      max_y = max(y,max_y)
      min_x = min(x,min_x)
      min_y = min(y,min_y)
    linenum += 1

  xdif = max_x-min_x
  ydif = max_y-min_y
  all_data['coords'] = coords
  all_data['coords_limits'] = [min_x-xdif/10.,max_x+xdif/10.,min_y-ydif/10.,max_y-ydif/10.]
  return all_data



def read_nominalvoltsfile(log,all_data):
  filename = all_data['nominalvoltsfilename']
  log.joint('reading coord file ' + filename + '\n')

  csvfile = open(filename,'rb')
  reader = csv.DictReader(csvfile)

  nominalvolt = {}

  for row in reader:
    bus = int(row['TermID'])
    nominalvolt[bus] = float(row['Nominal_Voltage'])

  all_data['nominalvolt']  = nominalvolt

  return all_data


 


