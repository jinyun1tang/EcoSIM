#libraries
import numpy as np
from netCDF4 import Dataset
import subprocess
from pathlib import Path
import json
import os

class ParEditor:
    # contruct the parameter editor ParEditor
    def __init__(self,pftparfile=None,micparfile=None):
        self.mic_recordjsonl=''
        self.pft_recordjsonl=''
        if pftparfile:
            self.pftparfile=pftparfile
            self.pft_recordjsonl=self.__get_jsonf_stem(pftparfile)
        if micparfile:    
            self.micparfile=micparfile        
            self.mic_recordjsonl=self.__get_jsonf_stem(micparfile)
        
    def reset(self,sure=False):
        """
        remove jsonl file if reset
        """
        if not sure:
            prompt='do you want remove the jasonl files that records past changes?'
            user_input = input(prompt).lower().strip()
            if user_input in ['y', 'yes']:
                self.__rmfile(self.mic_recordjsonl)
                self.__rmfile(self.pft_recordjsonl)
                print('Old jsonl files are removed.')
            else:
                print('Old jsonl files are used.')
        else:
            self.__rmfile(self.mic_recordjsonl)
            self.__rmfile(self.pft_recordjsonl)           
            print('Old jsonl files are removed.')

        
    def __rmfile(self,filepath):
        """
        remove filepath if it exists
        """
        filePath=filepath+'.jsonl'
        if os.path.exists(filePath):
            os.remove(filePath)    
            
    
    def __get_jsonf_stem(self,fname):
        """
        extract file stem for the jsonl file
        """
        parts0=fname.split('/')
        parts1=parts0[-1].split('.')
        return parts1[0]
        
                 
    def PlantParamModify(self,pft,pars,iscale=False,verbose=False):
        """
        Modify parameter parnames using parvals for pft on file parfile
        """
        try:
            print(self.pftparfile)
        except NameError:
            print("The pft parameter file is not defined.")
    
        new_dict={}
        with Dataset(self.pftparfile, 'r+') as nc_file:
            variable = nc_file.variables['pfts']        
            pft_loc=0
            #initialize list for recording parameter values
        
            for var in variable:
                result_string=''
                for byte in var:
                    if byte:
                        result_string=''.join([result_string,byte.decode('utf-8')])
                #locate the pft       
                if result_string.strip()==pft:                
                    #locate the variable
                    parvs= [0] * len(pars) 
                    parnames=['']*len(pars)
                    id=0
                    for parnm,parval in pars.items():
                        variable1=nc_file.variables[parnm]
                        long_name = variable1.getncattr('long_name') if 'long_name' in variable1.ncattrs() else 'No long_name attribute'
                    
                        if iscale:
                            if verbose:
                                print("%-80s: %s for %s is %f, and changed to %f"%(long_name,parnm,pft,variable1[pft_loc],parval*variable1[pft_loc]))
                            variable1[pft_loc]=parval*variable1[pft_loc]
                        else:
                            if verbose:
                                print("%-80s: %s for %s is %f, and changed to %f"%(long_name,parnm,pft,variable1[pft_loc],parval))
                            variable1[pft_loc]=parval
                        parvs[id]=float(variable1[pft_loc])
                        parnames[id]=parnm
                        id=id+1                    
                    new_dict={'pft':pft,'parvarnames':parnames,'parvals':parvs}
                    break                                                    
                pft_loc=pft_loc+1
        if new_dict:
            self.RecordPftPars(new_dict,self.pft_recordjsonl)
            

    def MicrobeParamModify(self,pars,iscale=False,verbose=False):
        """
        Modify microbial parameters parnames using parvals on file parfile
        """
        try:
            print(self.micparfile)
        except NameError:
            print("The microbial parameter file is not defined.")
            
        new_dict={}
        with Dataset(self.micparfile, 'r+') as nc_file:
            #initialize list for recording parameter values
            parvs=[0]*len(pars)
            parnames=['']*len(pars)
            id=0
            #locate the variable
            for parnm,parval in pars.items():
                variable1=nc_file.variables[parnm]
            
                long_name = variable1.getncattr('long_name') if 'long_name' in variable1.ncattrs() else 'No long_name attribute'
                if iscale:
                    if verbose:
                        print("%80s: %s is %f, and changed to %f"%(long_name,parnm,variable1[:],parval*variable1[:]))
                    variable1[:]=parval*variable1[:]                
                else:
                    if verbose:
                        print("%80s: %s is %f, and changed to %f"%(long_name,parnm,variable1[:],parval))
                    variable1[:]=parval
                parvs[id]=float(variable1[:])
                parnames[id]=parnm
                id=id+1
        new_dict={'parvarnames':parnames,'parvals':parvs}
        self.RecordPftPars(new_dict,self.mic_recordjsonl)
        

    def RecordPftPars(self,new_dict,RecordFileName):
        #record the old parameters
        # 'a' means append mode
        with open(RecordFileName+'.jsonl', 'a') as file:
            # Convert dict to a JSON string and write it, followed by a newline
            file.write(json.dumps(new_dict) + '\n')

