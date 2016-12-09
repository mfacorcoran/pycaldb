def cif2caldbwwwprod (telescop, instrume, version,
                      workdir="/FTP/caldb/local/software/pro/DATA_DELIVERIES/pro/mission_summaries/work",
                      caldbwww="/www.prod/htdocs/docs/heasarc/caldb/data"):
    """
    moves html summary version of calibration index file to caldbwwwprod
    :param telescop:
    :param instrume:
    :param version:
    :return:
    """

    import os
    import shutil

    cifname = "cif_"+telescop.strip().lower()+"_"+instrume.strip().lower()+"_"+version+".html"
    cifwwwdir = caldbwww+"/"+telescop.strip().lower()+"/"+instrume.strip().lower()+'/index'
    shutil.copy(workdir+"/"+cifname, cifwwwdir+"/"+cifname)
    os.symlink(cifwwwdir+"/"+cifname,cifwwwdir+"/index.html")

    return
