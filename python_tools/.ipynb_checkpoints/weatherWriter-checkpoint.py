import numpy as np
import numpy as np
from netCDF4 import Dataset
import os,time,sys,argparse
from datetime import datetime
from array import array
import warnings
warnings.filterwarnings("ignore")


case_sample_dict={
'dir':'sample',
'flist':'weth_list',
'ttype':'THWPR',#the five variables: temperature, humidity, wind, precip, solar_radiation
'ctype':'KRSMW',#the units: Kelvin, relative humidity, m/s, mm/m2, W m-2
}
mdir='/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/'

case=case_sample_dict
flist=mdir+case['dir']+'/'+case['flist']

current_dateTime = datetime.now()
nc_f=case['dir']+'_weather_%4d%02d%02d.nc'%(current_dateTime.year,current_dateTime.month,current_dateTime.day)
nc_fid = Dataset(nc_f, 'w')

ngrid=1
nc_fid.createDimension('year', None)
nc_fid.createDimension('day', 366)
nc_fid.createDimension('hour', 24)
nc_fid.createDimension('ngrid', ngrid)


with open(flist,"r") as fnm:
    line=fnm.readline().strip()
    while line:
        print("%s\nProcessing file %s"%('='*50,line))
        fwethnm=mdir+case['dir']+'/'+line
        with open(fwethnm,"r") as fweth:
            line1=fweth.readline().strip()
#READ(3,'(2A1,2I2,50A1)')TTYPE,CTYPE,NI,NN,(IVAR(K),K=1,NI),(VAR(K),K=1,NN)
            print(line1)
            line1=fweth.readline().strip()
#    READ(3,'(50A1)')(TYP(K),K=1,NN)
            print(line1)
            line1=fweth.readline().strip()
            datv=line1.split(',')
            Z0G=float(datv[0])
            IFLGW=int(datav(1))
            ZNOONG=float(datv[2])
            print(line1)
            line1=fweth.readline().strip()
#    READ(3,*)PHRG,CN4RIG,CNORIG,CPORG,CALRG,CFERG,CCARG,CMGRG,CNARG &
#      ,CKARG,CSORG,CCLRG
            print(line1)
            line1=fweth.readline().strip()
            while line1:
                dat=line1.split(',')
                print(dat)
                line1=fweth.readline().strip()

        line=fnm.readline().strip()

nc_fid.close()  # close the new file
