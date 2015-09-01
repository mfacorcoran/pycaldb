__author__ = 'mcorcora'

from utils import list_ftpdir

def get_caldb_missions(caldb="ftp://heasarc.gsfc.nasa.gov/caldb"):
    """
    get list of all missions listed in $CALDB/data
    """
    ftpurl=caldb+'/data'
    mission=list_ftpdir(ftpurl)
    mission.remove('README')
    mission.remove('Storage')
    mission.remove('pcf')
    test=list(mission) # test=mission won't work, since then they are the same object so removing an element from test removes it from mission
    for item in test:
        if 'chandra' in item:
            mission.remove(item)
    mission.append('chandra')
    mission.sort()
    return  mission

def get_caldb_instruments(mission, caldb="ftp://heasarc.gsfc.nasa.gov/caldb"):
    instdir=caldb+"/data/"+mission.strip().lower()
    instruments=list_ftpdir(instdir)
    try:
        instruments.remove('README')
        instruments.remove('Storage')
    except:
        pass
    return instruments

def get_supported_missions_dict(caldb='ftp://heasarc.gsfc.nasa.gov/caldb'):
    """
    Creates a dictionary for the HEASRC caldb of all supported missions and instruments
    """
    caldb_dict=dict()
    missions=get_caldb_missions()
    for mis in missions:
        inst=get_caldb_instruments(mis)
        caldb_dict[mis]=inst
    return caldb_dict



if __name__ == "__main__":
    caldb_dict=get_supported_missions_dict(caldb="ftp://heasarc.gsfc.nasa.gov/caldb")
    for key in caldb_dict.keys():
        print "Mission={0}, instruments={1}\n".format(key,caldb_dict[key])

