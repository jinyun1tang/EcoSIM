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


config_dict=config_DaringLake_dict

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


from SiteTopoWriter import write_site_topo_data

config_dict=config_DaringLake_dict
write_site_topo_data(config_dict)