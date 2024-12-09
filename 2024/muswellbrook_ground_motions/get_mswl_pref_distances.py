from mapping_tools import distance

def return_all_au_station_data():
    from datetime import datetime
    from os import getcwd
    
    mswl_station_file = 'mswl_station_data.dat'
    
    lines = open(mswl_station_file).readlines()[1:]
    
    sta_dict = []
    
    for line in lines:
        dat = line.strip().split('\t')
        #print(line)
        
        if int(dat[5]) < 1:
            dat[5] = 1
        if int(dat[7]) < 1:
            dat[7] = 1
            
        tmp = {'sta':dat[0], 'stlo':float(dat[1]), 'stla':float(dat[2]), 
               'startdate':datetime(int(dat[4]),int(dat[5]),1), 
               'enddate':datetime(int(dat[6]),int(dat[7]),1)}
        
        # append to sta_dict
        sta_dict.append(tmp)
        
    return sta_dict
    
sta_dict = return_all_au_station_data()

eqla = -32.347
eqlo = 150.869

for sta in sta_dict:
    rng = distance(eqla, eqlo, sta['stla'], sta['stlo'])[0]
    
    print(','.join((sta['sta'], str(rng))))
    