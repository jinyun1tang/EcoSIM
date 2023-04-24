import numpy as np
import numpy as np
from netCDF4 import Dataset
import os,time,sys,argparse
from datetime import datetime
from array import array
import warnings
warnings.filterwarnings("ignore")


"""
define auxillary functions
"""
def split_var(var_info):
  x=var_info.split(':')
  return x

def string2arr(string,len1=0):
  """
  conver a string into array of chars
  """
  if len1 == 0:
    to_array = []
    for x in string:
      to_array.extend(x)
  else:
    to_array=np.array([' '] * len1, dtype='S1')
    ll=0
    for x in string:
      to_array[ll]=x
      ll=ll+1

  return to_array


def strpack(strarr):
  """
  pack the string array, to remove empty elements
  """
  to_array = []
  for x in strarr:
    if x:
      to_array.append(x)
  return to_array

#writing pft management data, annual time series


topo_unit_info={
'NH1':'Starting column from the west:none:i1:0d',
'NH2':'Ending column at the east:none:i1:0d',
'NV1':'Starting row from the north:none:i1:0d',
'NV2':'Ending row at the south:none:i1:0d'
}


#        READ(12,*,END=540)DY,PPI(NZ,NY,NX),SDPTHI(NZ,NY,NX)
#        READ(12,*,END=540)DY,ICUT,JCUT,HCUT,PCUT &
#          ,ECUT11,ECUT12,ECUT13,ECUT14,ECUT21,ECUT22,ECUT23,ECUT24

pft_mgmt_info={
'PLTDY':'Planting date DDMMYYYY:none:i32',      #PLTDY(NZ,NY,NX)
'PPI':'Initial planting density:100> for grass, ~1 for trees:m-2:f4', #PPI(NZ,NY,NX)
'SDPTHI':'Seeding depth:m:f4',           #SDPTHI(NZ,NY,NX)
'MGMTDY':'Management date DDMMYYYY:could be harvest, grazing etc:none:i32',   #MGMTDY(ntimes,NZ,NY,NX)
'ICUT':'Harvest type:0=none,1=grain,2=all above-ground,3=pruning,4=grazing,5=fire,6=herbivory:none:i1', # (ntimes,NZ,NY,NX)
'JCUT':'Flag for PFT termination:0=no,1=yes,2=yes,but reseed:none:i1', #(ntimes,NZ,NY,NX)
'PCUT':'Disturbance rate:if ICUT>=0, =cutting height, else =fraction' \
      +' of LAI removed; if grazing, grazer consumption rate:g DW g FW-1 d-1:f4',#(ntimes,NZ,NY,NX)
'ECUT11':'Fraction of leaf removed from PFT:none:f4', #(ntimes,NZ,NY,NX)
'ECUT12':'Fraction of non-foliar removed from PFT:none:f4', #(ntimes,NZ,NY,NX)
'ECUT13':'Fraction of wood removed from PFT:none:f4', #(ntimes,NZ,NY,NX)
'ECUT14':'Fraction of standing dead removed from PFT:none:f4',#(ntimes,NZ,NY,NX)
'ECUT21':'Fraction of leaf removed from ecosystem:none:f4',#(ntimes,NZ,NY,NX)
'ECUT22':'Fraction of non-foliar removed from ecosystem:none:f4',#(ntimes,NZ,NY,NX)
'ECUT23':'Fraction of woody removed from ecosystem:none:f4',#(ntimes,NZ,NY,NX)
'ECUT24':'Fraction of standing dead removed from ecosystem:none:f4' #(ntimes,NZ,NY,NX)}
}

config_dry_dict={
'case':'dryland',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/dryland_maize/',
'pftf':'me01p:me02p:me03p:me04p:me05p:me06p',
'topf':'metopo',
'ntopu':'1',
'ncol':'1',
'nrow':'1'
}

config_lake_dict={
'case':'lake',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/lake/',
'pftf':'va99p:vaxxp',
'topf':'vatopo',
'ntopu':'2',
'ncol':'2',
'nrow':'1'
}

config_sample_dict={
'case':'sample',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/sample/',
'pftf':'pft_arctic_p:pft_arctic_g',
'topf':'tp022852',
'ntopu':'1',
'ncol':'1',
'nrow':'1'
}

config_dict=config_lake_dict

print('generate pft data for '+config_dict['case'])

mdir=config_dict['mdir']
case=config_dict['case']
pfts=split_var(config_dict['pftf'])

ntopu=int(config_dict['ntopu'])
ncol=int(config_dict['ncol'])
nrow=int(config_dict['nrow'])

current_dateTime = datetime.now()
nc_f=case+'_pft_%4d%02d%02d.nc'%(current_dateTime.year,current_dateTime.month,current_dateTime.day)

print("check file %s"%nc_f)

nc_fid = Dataset(nc_f, 'w')

nc_fid.description='PFT input data created on %4d/%02d/%02d/%02d:%02d:%02d'% \
  (current_dateTime.year,current_dateTime.month,current_dateTime.day, \
  current_dateTime.hour,current_dateTime.minute,current_dateTime.second)+ \
  '\n use READ(tline,*)DY,PPI(NZ,NY,NX),SDPTHI(NZ,NY,NX) to read planting information from pft_pltinfo;' +\
  ' use READ(tline,*)DY,ICUT,JCUT,HCUT,PCUT,ECUT11,ECUT12,ECUT13,ECUT14,ECUT21,ECUT22,ECUT23,ECUT24 to '+\
  'read management information from pft_mgmt'

nc_fid.createDimension('ntopou', ntopu)
nc_fid.createDimension('nchar1', 10)
nc_fid.createDimension('nchar2', 10)
nc_fid.createDimension('ncharmgnt',128)
nc_fid.createDimension('maxpfts', 5)
nc_fid.createDimension('maxpmgt',24) #management is done at most once every two weeks
nc_fid.createDimension('year', None)
nc_fid.createDimension('ncol', ncol)
nc_fid.createDimension('nrow', nrow)

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

w_nc_var = nc_fid.createVariable('NZ', 'i4', ('ntopou'))
w_nc_var.long_name='Number of pfts on a topo unit'
w_nc_var.units='None'

w_nc_var = nc_fid.createVariable('nmgnts', 'i2', ('year','ntopou','maxpfts'))
w_nc_var.long_name='Number of managements for a given pft in in given topo unit in a year'
#w_nc_var = nc_fid.createVariable('pft_mgmt', 'S1', ('ntopu','nchars'))
#w_nc_var = nc_fid.createVariable('pft_mgmt', 'S1', ('ntopu','nchars'))
w_nc_var = nc_fid.createVariable('pft_type', 'S1', ('year','ntopou','maxpfts','nchar1'))

w_nc_var = nc_fid.createVariable('pft_pltinfo', 'S1', ('year','ntopou','maxpfts','ncharmgnt'))
w_nc_var.long_name='string containing planting information'


w_nc_var = nc_fid.createVariable('pft_mgmt', 'S1', ('year','ntopou','maxpfts','maxpmgt','ncharmgnt'))
w_nc_var.long_name='string containing plant management information'

pflag = nc_fid.createVariable('pft_dflag', 'i4', ())

k=0
pft_dflag=-1  #no pft data
if len(pfts)==2:
  pft_dflag=0  #constant pft managment data, only plantation
elif len(pfts)>2:
  pft_dflag=1  #transient pft managment data, including multiple
pflag[:]=pft_dflag
pflag.setncattr('long_name','Flag for plant management data')
pflag.setncattr('flags','-1 no pft data, 0 only plantation information, 1 transient pft data')

#for v in pft_mgmt_info:
#  print(v)
#  ss=split_var(pft_mgmt_info[v])
#  if len(ss)==4:
#    long_name,flags,units,dtype=ss
#  else:
#    flags=''
#    long_name,units,dtype=ss
#  print(long_name)
#  print(units)

k1=0
for pf in pfts:
  pftnm=mdir+pf
  k=0
  with open(pftnm,"r") as pftf:
    while 1:
      line=pftf.readline()
      if not line:
        break
      arr=strpack(line.strip().split(' '))
      NH1=int(arr[0])
      NV1=int(arr[1])
      NH2=int(arr[2])
      NV2=int(arr[3])
      NZ=int(arr[4])

      nc_fid.variables['NH1'][k]=NH1
      nc_fid.variables['NV1'][k]=NV1
      nc_fid.variables['NH2'][k]=NH2
      nc_fid.variables['NV2'][k]=NV2
      nc_fid.variables['NZ'][k]=NZ
      line=pftf.readline()
      arr=strpack(line.strip().split(' '))
      for nz in range(NZ):
        sarr=string2arr(arr[nz*2])
        len1=len(sarr)
        nc_fid.variables['pft_type'][k1,k,nz,0:len1]=sarr
        if arr[nz*2+1] != 'NO':
          print('%s %s\n'%(arr[nz*2],arr[nz*2+1]))
          mgfnm=mdir+arr[nz*2+1]
          with open(mgfnm,"r") as mgf:
            print(mgfnm)
            line1=mgf.readline()
            sarr1=string2arr(line1.strip())
            len1=len(sarr1)
            nc_fid.variables['pft_pltinfo'][k1,k,nz,:]=' '*128
            nc_fid.variables['pft_pltinfo'][k1,k,nz,0:len1]=sarr1
            k2=0
            while 1:
              line1=mgf.readline()
              if not line1:
                break
              sarr1=string2arr(line1.strip())
              len1=len(sarr1)
              nc_fid.variables['pft_mgmt'][k1,k,nz,k2,:]=' '*128
              nc_fid.variables['pft_mgmt'][k1,k,nz,k2,0:len1]=sarr1
              k2=k2+1
              print(line1)
          nc_fid.variables['nmgnts'][k1,k,nz]=k2
      k=k+1
  k1=k1+1
nc_fid.close()  # close the new file
