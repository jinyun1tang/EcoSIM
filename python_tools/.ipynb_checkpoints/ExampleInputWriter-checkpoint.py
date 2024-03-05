from pftMgmtWriter import write_pft_mgmt

config_dry_dict={
'case':'dryland',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/dryland_maize/',
'pftf':'me01p:me02p:me03p:me04p:me05p:me06p',
'year':'2001:2002:2003:2004:2005:2006',
#'topf':'metopo',
'ntopu':'1',
'ncol':'1',
'nrow':'1'
}

config_lake_dict={
'case':'lake',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/lake/',
'pftf':'va99p:vaxxp',
#'topf':'vatopo',
'ntopu':'2',
'ncol':'2',
'nrow':'1'
}

config_sample_dict={
'case':'sample',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/sample/',
'pftf':'pft_arctic_p:pft_arctic_g',
#'topf':'tp022852',
'ntopu':'1',
'ncol':'1',
'nrow':'1'
}

config_DaringLake_dict={
'case':'DaringLake',
'mdir':'/Users/jinyuntang/work/github/ecosys2/ecosys_benchmark_runs_2022_oct_24/DaringLake/',
'pftf':'dlmt00p:dlmtxxp',
#'topf':'dlmtto',
'ntopu':'2'
}

config_Fen_StordIsland_dict={
'case':'FenStordIsland',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/fenrun_shuai/',
'pftf':'sdlxxp2:sdlxxp_no',
#'topf':'dlmtto',
'ntopu':'1',
'ncol':'1',
'nrow':'1'

}


config_MeditPastureCA_dict={
'case':'MeditPastureCA',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/Meditteranean_Pasture_CA/',
'pftf':'va99p:vaxxp',
#'topf':'dlmtto',
'ntopu':'1',
'ncol':'1',
'nrow':'1'

}

config_US_Ton_dict={
'case':'US_Ton',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/US_Ton/',
'pftf':'pft_uston_p:pft_uston_g',    

#'topf':'dlmtto',
'ntopu':'1',
'ncol':'1',
'nrow':'1'

}


config_SemiaridGrassland_dict={
'case':'SemiaridGrassland',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/Semiarid_Grassland_AB/',
'pftf':'lepg:lepx',    
#'topf':'dlmtto',
'ntopu':'1',
'ncol':'1',
'nrow':'1'
}    

config_WarmTempOakTN_dict={
'case':'WarmTempOakTN',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/Warm_Temperate_Oak_TN/',
'pftf':'td60p:tdxxp',    
'ntopu':'1',
'ncol':'1',
'nrow':'1'
}    


config_SnodgrassTransect_dict={
'case':'SnodgrassTransect',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/Snodgrass_transect/',
'pftf':'pft_ert_p:pft_ert_g',
#'topf':'dlmtto',
'ntopu':'6'
}    
case=4

if case==1:
    config_dict=config_Fen_StordIsland_dict
elif case==2:    
    config_dict=config_WarmTempOakTN_dict
elif case==3:    
    config_dict=config_US_Ton_dict
elif case==4:    
    config_dict=config_SemiaridGrassland_dict
elif case==5:    
    config_dict=config_MeditPastureCA_dict
elif case==6:
    config_dict=config_SnodgrassTransect_dict

write_pft_mgmt(config_dict)


#write topgraphy and site data

"""
main body of code
"""

#topo unit file
config_dry_dict={
'case':'dryland',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/dryland_maize/',
'sitef':'mesite',
'topf':'metopo',
'ntopu':'1',
'ncol':'1',
'nrow':'1',
'NHW':'1',
'NVN':'1',
'NHE':'1',
'NVS':'1'
}

config_lake_dict={
'case':'lake',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/lake/',
'sitef':'vasite',
'topf':'vatopo',
'ntopu':'2',
'ncol':'2',
'nrow':'1',
'NHW':'1',
'NVN':'1',
'NHE':'2',
'NVS':'1'
}

config_sample_dict={
'case':'sample',
'mdir':'/Users/jinyuntang/work/github/ecosim2/EcoSIM/examples/inputs/sample/',
'sitef':'st022852',
'topf':'tp022852',
'ntopu':'1',
'ncol':'1',
'nrow':'1',
'NHW':'1',
'NVN':'1',
'NHE':'1',
'NVS':'1'
}


config_DaringLake_dict={
'case':'DaringLake',
'mdir':'/Users/jinyuntang/work/github/ecosys2/ecosys_benchmark_runs_2022_oct_24/DaringLake/',
'sitef':'dlmts',
'topf':'dlmtto',
'ntopu':'6',
'ncol':'1',
'nrow':'6',
'NHW':'1',
'NVN':'1',
'NHE':'1',
'NVS':'6'
}

config_Fen_StordIsland_dict={
'case':'FenStordIsland',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/fenrun_shuai/',
'sitef':'sdlsite',
'topf':'sdltopo',
'ntopu':'1',
'ncol':'1',
'nrow':'1',
'NHW':'1',
'NVN':'1',
'NHE':'1',
'NVS':'1'
}

config_US_Ton_dict={
'case':'US_Ton',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/US_Ton/',
'sitef':'st_uston',
'topf':'tp_uston',
'ntopu':'1',
'ncol':'1',
'nrow':'1',
'NHW':'1',
'NVN':'1',
'NHE':'1',
'NVS':'1'
}

config_MeditPastureCA_dict={
'case':'MeditPastureCA',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/Meditteranean_Pasture_CA/',
'sitef':'vasite',
'topf':'vatopo',
'ntopu':'1',
'ncol':'1',
'nrow':'1',
'NHW':'1',
'NVN':'1',
'NHE':'1',
'NVS':'1'
}

config_SemiaridGrassland_dict={
'case':'SemiaridGrassland',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/Semiarid_Grassland_AB/',
'sitef':'lsite',
'topf':'ltopo',
'ntopu':'1',
'ncol':'1',
'nrow':'1',
'NHW':'1',
'NVN':'1',
'NHE':'1',
'NVS':'1'
}

config_WarmTempOakTN_dict={
'case':'WarmTempOakTN',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/Warm_Temperate_Oak_TN/',
'sitef':'tds',
'topf':'tdtopo',
'ntopu':'1',
'ncol':'1',
'nrow':'1',
'NHW':'1',
'NVN':'1',
'NHE':'1',
'NVS':'1'
}

config_SnodgrassTransect_dict={
'case':'SnodgrassTransect',
'mdir':'/Users/jinyuntang/work/ecosys_benchmark/Snodgrass_transect/',
'sitef':'st074773_e',
'topf':'tp074773_e',    
'ntopu':'6',
'ncol':'6',
'nrow':'1',
'NHW':'1',
'NVN':'1',
'NHE':'6',
'NVS':'1'
}


from SiteTopoWriter import write_site_topo_data

if case==1:
    config_dict=config_Fen_StordIsland_dict
elif case==2:    
    config_dict=config_WarmTempOakTN_dict
elif case==3:    
    config_dict=config_US_Ton_dict
elif case==4:    
    config_dict=config_SemiaridGrassland_dict
elif case==5:    
    config_dict=config_MeditPastureCA_dict
elif case==6:
    config_dict=config_SnodgrassTransect_dict
write_site_topo_data(config_dict)