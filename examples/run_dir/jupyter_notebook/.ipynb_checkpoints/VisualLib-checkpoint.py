import numpy as np

def get_hist_yr(histfile):
    """
    get the first year of history file
    """
    index = histfile.find('h0.')
    substring_after_h0 = histfile[index + len('h0.'):]    
    year=int(substring_after_h0[:4])
    return year

def is_leap_year(year):
    if (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0):
        return 1
    else:
        return 0

def get_ts_nyears(ts):
    return int(len(ts)/365)
    
def cumts2ts(cts,year0):
    """
    convert cumulative time series into difference
    """
    nyears=get_ts_nyears(cts)
    print('nyears=%d'%nyears)
    ts=np.zeros((nyears*366,1))
    id0=0
    id1=365
    nmels=len(cts)
    for yr in range(nyears):
        year=yr+year0
        print('year=%d'%year)
        leap=is_leap_year(year)
        id1=min(id1+leap,nmels)        
        ts[id0+1:id1]=cts[id0+1:id1]-cts[id0:id1-1]
        id0=id1
        id1=id1+365
    return ts
