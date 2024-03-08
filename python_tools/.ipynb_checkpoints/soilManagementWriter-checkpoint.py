import numpy as np
import numpy as np
from netCDF4 import Dataset
import os,time,sys,argparse
from datetime import datetime
from array import array
import stringTools as strtool
import warnings
warnings.filterwarnings("ignore")


config_lake_dict={
'case':'lake',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/lake/',
'outdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/lake/',    
'ntopu':'2'
}


config_sample_dict={
'case':'sample',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/sample/',
'outdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/sample/',    
'ntopu':'1'
}

config_dry_dict={
'case':'dryland',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/dryland_maize/',
'outdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/dryland_maize/',    
'ntopu':'1'
}


config_dict=config_lake_dict


def write_soil_mgmt(config_dict):
  print('generate soil managment data for '+config_dict['case'])

  mdir=config_dict['mdir']
  case=config_dict['case']

  ntopu=int(config_dict['ntopu'])

  current_dateTime = datetime.now()
  nc_f=config_dict['mdir']+case+'_soilmgmt_%4d%02d%02d.nc'%(current_dateTime.year,current_dateTime.month,current_dateTime.day)

  print("check file %s"%nc_f)

  nc_fid = Dataset(nc_f, 'w')

  nc_fid.description='soil managment data created on %4d/%02d/%02d/%02d:%02d:%02d'% \
    (current_dateTime.year,current_dateTime.month,current_dateTime.day, \
    current_dateTime.hour,current_dateTime.minute,current_dateTime.second)


  nc_fid.createDimension('ntopou', ntopu)
  nc_fid.createDimension('nchar1', 10)
  nc_fid.createDimension('nchart', 24)
  nc_fid.createDimension('ntill', 12)
  nc_fid.createDimension('nfert', 12)
  nc_fid.createDimension('ncharf', 128)
  nc_fid.createDimension('year', None)

  w_nc_var = nc_fid.createVariable('year', 'i4', ('year'))
  w_nc_var.long_name='year AD'

  w_nc_var = nc_fid.createVariable('NH1', 'i4', ('ntopou'))
  w_nc_var.long_name='Starting column from the west for a topo unit'
  w_nc_var.units='None'

  w_nc_var = nc_fid.createVariable('NV1', 'i4', ('ntopou'))
  w_nc_var.long_name='Ending column at the east for a topo unit'
  w_nc_var.units='None'

  w_nc_var = nc_fid.createVariable('NH2', 'i4', ('ntopou'))
  w_nc_var.long_name='Starting row from the north  for a topo unit'
  w_nc_var.units='None'

  w_nc_var = nc_fid.createVariable('NV2', 'i4', ('ntopou'))
  w_nc_var.long_name='Ending row at the south  for a topo unit'
  w_nc_var.units='None'

  w_nc_var = nc_fid.createVariable('fertf', 'S1', ('year','ntopou','nchar1'))
  w_nc_var.long_name='Fertilization info for a topo unit'

  w_nc_var = nc_fid.createVariable('tillf', 'S1', ('year','ntopou','nchar1'))
  w_nc_var.long_name='Tillage info for a topo unit'

  w_nc_var = nc_fid.createVariable('irrigf', 'S1', ('year','ntopou','nchar1'))
  w_nc_var.long_name='Irrigation info for a topo unit'


  if case=='dryland':
    w_nc_var=nc_fid.createVariable('me01t', 'S1', ('ntill','nchart'))
    w_nc_var.long_name='Tillage file'

    nc_fid.variables['NH1'][:]=1
    nc_fid.variables['NV1'][:]=1
    nc_fid.variables['NH2'][:]=1
    nc_fid.variables['NV2'][:]=1

    mfname='me01t'
    ifile=mdir+mfname
    with open(ifile,"r") as infile:
      line=infile.readline().strip()
      for k1 in range(12):
        nc_fid.variables[mfname][k1,:]=' '*24
      k1=0
      while line:
        sarr1=strtool.string2arr(line)
        len1=len(sarr1)
        nc_fid.variables[mfname][k1,0:len1]=sarr1
        line=infile.readline().strip()
        k1=k1+1

    w_nc_var=nc_fid.createVariable('me01f', 'S1', ('nfert','ncharf'))
    w_nc_var.long_name='fertilization file'
    mfname='me01f'
    ifile=mdir+mfname
    with open(ifile,"r") as infile:
      line=infile.readline().strip()
      for k1 in range(12):
        nc_fid.variables[mfname][k1,:]=' '*24

      k1=0
      while line:
        sarr1=strtool.string2arr(line)
        len1=len(sarr1)
        nc_fid.variables[mfname][k1,:]=' '*128
        nc_fid.variables[mfname][k1,0:len1]=sarr1
        line=infile.readline().strip()
        k1=k1+1

    w_nc_var=nc_fid.createVariable('me03f', 'S1', ('nfert','ncharf'))
    w_nc_var.long_name='fertilization file'
    mfname='me03f'
    ifile=mdir+mfname
    with open(ifile,"r") as infile:
      line=infile.readline().strip()
      for k1 in range(12):
        nc_fid.variables[mfname][k1,:]=' '*24
      k1=0
      while line:
        sarr1=strtool.string2arr(line)
        len1=len(sarr1)
        nc_fid.variables[mfname][k1,:]=' '*128
        nc_fid.variables[mfname][k1,0:len1]=sarr1
        line=infile.readline().strip()
        k1=k1+1

    w_nc_var=nc_fid.createVariable('me05f', 'S1', ('nfert','ncharf'))
    w_nc_var.long_name='fertilization file'
    mfname='me05f'
    ifile=mdir+mfname
    with open(ifile,"r") as infile:
      line=infile.readline().strip()
      for k1 in range(12):
        nc_fid.variables[mfname][k1,:]=' '*24
      k1=0
      while line:
        sarr1=strtool.string2arr(line)
        len1=len(sarr1)
        nc_fid.variables[mfname][k1,:]=' '*128
        nc_fid.variables[mfname][k1,0:len1]=sarr1
        line=infile.readline().strip()
        k1=k1+1

    nc_fid.variables['year'][:]=range(2001,2007)
    for jj in range(2007-2001):
      nc_fid.variables['fertf'][jj,0,:]=' '*10
      nc_fid.variables['tillf'][jj,0,:]=' '*10
      nc_fid.variables['irrigf'][jj,0,:]=' '*10

    nc_fid.variables['tillf'][0,0,0:5]=['m','e','0','1','t']
    nc_fid.variables['fertf'][0,0,0:5]=['m','e','0','1','f']
    nc_fid.variables['irrigf'][0,0,0:2]=['N','O']

    nc_fid.variables['tillf'][1,0,0:2]=['N','O']
    nc_fid.variables['fertf'][1,0,0:2]=['N','O']
    nc_fid.variables['irrigf'][1,0,0:2]=['N','O']

    nc_fid.variables['tillf'][2,0,0:2]=['N','O']
    nc_fid.variables['fertf'][2,0,0:5]=['m','e','0','3','f']
    nc_fid.variables['irrigf'][2,0,0:2]=['N','O']

    nc_fid.variables['tillf'][3,0,0:2]=['N','O']
    nc_fid.variables['fertf'][3,0,0:2]=['N','O']
    nc_fid.variables['irrigf'][3,0,0:2]=['N','O']

    nc_fid.variables['tillf'][4,0,0:2]=['N','O']
    nc_fid.variables['fertf'][4,0,0:5]=['m','e','0','5','f']
    nc_fid.variables['irrigf'][4,0,0:2]=['N','O']

    nc_fid.variables['tillf'][5,0,0:2]=['N','O']
    nc_fid.variables['fertf'][5,0,0:2]=['N','O']
    nc_fid.variables['irrigf'][5,0,0:2]=['N','O']

  elif case=='sample':
    nc_fid.variables['year'][:]=range(1979,2011)
    for jj in range(2011-1979):
      nc_fid.variables['fertf'][jj,0,:]=' '*10
      nc_fid.variables['tillf'][jj,0,:]=' '*10
      nc_fid.variables['irrigf'][jj,0,:]=' '*10
      nc_fid.variables['fertf'][jj,0,0:2]=['N','O']
      nc_fid.variables['tillf'][jj,0,0:2]=['N','O']
      nc_fid.variables['irrigf'][jj,0,0:2]=['N','O']
  elif case=='lake':
    nc_fid.variables['NH1'][:]=[1,2]
    nc_fid.variables['NV1'][:]=[1,1]
    nc_fid.variables['NH2'][:]=[1,2]
    nc_fid.variables['NV2'][:]=[1,1]
    nc_fid.variables['year'][:]=range(2001,2009)
    for jj in range(2009-2001):
      nc_fid.variables['fertf'][jj,:,:]=' '*10
      nc_fid.variables['tillf'][jj,:,:]=' '*10
      nc_fid.variables['irrigf'][jj,:,:]=' '*10
      nc_fid.variables['fertf'][jj,:,0:2]=[['N','O'],['N','O']]
      nc_fid.variables['tillf'][jj,:,0:2]=[['N','O'],['N','O']]
      nc_fid.variables['irrigf'][jj,:,0:2]=[['N','O'],['N','O']]

  nc_fid.close()  # close the new file
