def patch_nustar_caldb(version,
                       nupatchserver = "hassif.caltech.edu",
                       nupatchworkdir = "/pub/nustar/pipeline",
                       caldbstage = "/web_chroot/FTP/caldb/staging",
                       ck_file_exists='True'):
    """
    This function creates a caldb_update.txt message as described in
    https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_2003_001/cal_gen_2003_001.html

    Once this update file is created, the caldbmgr needs to run
    /FTP/caldb/local/scripts/DATA_DELIVERIES/nustar_caldbupdate.csh

    in order to complete the update.

    See /web_chroot/FTP/.caldb_staging/data/nustar/README_rev1.txt for more detailed information on how such updates
    were done manually.

    After this function runs successfully, then, to ingest the data (as in other missions),
    manually run /FTP/caldb/local/scripts/DATA_DELIVERIES/nustar_caldbupdate.csh

    :param version: string in 'YYYYMMDD' format, from nustar patch notification (from BG)
    :param nupatchworkdir: URL of directory where patch tar file is located, from BG
    :param caldbstage: location of caldb staging area
    :param ck_file_exists: Checks that the new file is in the new caldb.indx file (Default =  'True')
    :return:
    """

    import os
    import ftputil
    import tarfile
    import datetime
    import shutil
    from heasarc.pycaldb.pycaldb import ck_file_existence
    from astropy.io import fits as pyfits

    #
    # create patch install directory
    #
    versiondir=caldbstage+"/data/nustar/versions/"+version.strip().lower()
    patchinstalldir = versiondir
    if os.path.isdir(versiondir):
        ans = raw_input("{0} exists; Do you want to remove it (Y/n)? ".format(versiondir))
        try:
            if ans.upper().strip()[0] <> "N":
                shutil.rmtree(versiondir)
                print("Creating clean {0} directory".format(patchinstalldir))
                os.makedirs(patchinstalldir)
        except IndexError: # user hit <CR> generating an empty string
            shutil.rmtree(versiondir)
            print("Creating clean {0} directory".format(patchinstalldir))
            os.makedirs(patchinstalldir)
    else:
        print("Creating {0} directory".format(patchinstalldir))
        os.makedirs(patchinstalldir)
    #
    #create a full version of the existing nustar caldb in the
    #    $CALDB/data/nustar in the nustar/versions are via rsync:
    # cdirs=['bcf','cpf','index', 'caldb.indx']
    # for c in cdirs:
    #     cmd = "rsync -avz " + caldb + "/data/nustar/fpm/" + c + " " + patchinstalldir + "/fpm"
    #     if os.path.isdir(patchinstalldir+'/fpm/'+c):
    #         ans = raw_input('\nDirectory '+patchinstalldir+'/fpm/'+c+' Exists.  Remove it and re-download (Y/n)? ')
    #         if ans.lower().strip()[0] == 'y':
    #             shutil.rmtree(patchinstalldir+'/'+c)
    #         else:
    #             cmd = "echo keeping pre-existing directory "+patchinstalldir+'/fpm/'+c
    #     print (cmd)
    #     os.system(cmd)
    #
    # get the patch file
    #
    host = ftputil.FTPHost(nupatchserver,'anonymous','caldbmgr@nasa.gov')
    patchfiledir = nupatchworkdir+'/'+version.strip()
    patchfile = 'NUSTARCALDB_patch'+version.strip()+'.tgz'
    print('\nDownloading {0}...'.format(patchfile))
    try:
        host.download(patchfiledir+'/'+patchfile,patchinstalldir+'/'+patchfile,mode='b') # need mode flag for ftputil<3.0
    except TypeError:
        host.download(patchfiledir + '/' + patchfile, patchinstalldir + '/' + patchfile)  # if mode flag not present (ftputil>=3.0)
    #
    # untar the file
    #
    print('\nExtracting {0} to {1}'.format(patchfile, patchinstalldir))
    ptf = tarfile.open(patchinstalldir+'/'+patchfile)
    ptf.extractall(patchinstalldir)
    tarlist = ptf.getnames()
    print "Tarlist:"
    print tarlist
    newfiles = [x for x in tarlist if os.path.isfile(x)]
    print "Newfiles:"
    print newfiles
    #
    # Link the updated version to staging/data/nustar/fpm:
    #
    #
    dst = caldbstage+'/data/nustar/fpm'
    if os.path.islink(dst):
        os.remove(dst)
    print "Creating symbolic link from {0} to {1}".format(patchinstalldir+"/fpm", dst)
    try:
        os.symlink(patchinstalldir+'/fpm', dst)
    except:
        print "Problem creating symbolic link to "+dst+"; Returning"
        return
    #
    # check that all files in the cif exist in the new caldb
    #
    if ck_file_exists:
        print ("\nChecking existence of files in the caldb.indx file...")
        cif = pyfits.open(caldbstage+'/data/nustar/fpm/caldb.indx')
        missing = ck_file_existence(cif, caldb=caldbstage, quiet=True)
        if len(missing) > 0:
            print "\nThese files are missing from the caldb:"
            for f in missing:
                print f
    #
    # create update file in staging/nustar/to_ingest directory
    #
    update_txt='from: caldbmgr@nasa.gov\n' \
               'Subject: caldb update nustar '+version.strip()+'\n' \
               'Date: '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")+'\n' \
               'To: caldbingest@bigbang.gsfc.nasa.gov\n\n' \
               'instrument=fpm\n' \
               'caldbindexfile=caldb.indx'+version.strip()+'\n\n' \
               'newfile\n'
    for f in newfiles: # get file list from patch tarfile
        try:
            name=f.split('fpm/')[1]
        except IndexError, errmsg:
            print "Could not find a file in tarfile entry {0} ({1})".format(f, errmsg)
        if not 'caldb.indx' in name.strip().lower(): # don't include caldb.indx files in ingest note
            update_txt=update_txt+name+"\n"
    update_txt=update_txt+'endnewfile\n'
    to_ingestdir= os.path.join(caldbstage,'data/nustar/to_ingest')
    if not os.path.exists(to_ingestdir):
        print "Creating {0}".format(to_ingestdir)
        os.makedirs(to_ingestdir)
    ingestfilename = 'caldb_update_nustar_'+version.strip()+'.txt'
    ingestfile = caldbstage+'/data/nustar/to_ingest/'+ingestfilename
    print('Creating {0}'.format(ingestfile))
    fingest = open(ingestfile,'w')
    fingest.write(update_txt)
    fingest.close()
    #print('\nEdit {0} and enter list of new files from the update e-mail from BG'.format(ingestfile))
    #
    # make symlink to caldb_update.txt file
    #
    curdir=os.getcwd()
    to_ingestdir = caldbstage+"/data/nustar/to_ingest"
    updatefile='caldb_update.txt'
    os.chdir(to_ingestdir)
    print "Creating symbolic link from {0} to {1} in {2}".format(ingestfilename, updatefile,to_ingestdir)
    if os.path.islink(updatefile):
        os.remove(updatefile)
    try:
        os.symlink(ingestfile, updatefile)
    except:
        print "Problem creating symbolic link to "+updatefile+"; Returning"
        return
    os.chdir(curdir)
    print "To complete the FPM update, run /FTP/caldb/local/scripts/DATA_DELIVERIES/nustar_caldbupdate.csh"
    return



if __name__ == "__main__":
    #
    # for sshfs mounted directories
    #
    # caldb="/fuse/caldb"
    # stage = 'caldb_stage'
    caldb = "/web_chroot/FTP/caldb"
    stage = "/web_chroot/FTP/caldb/staging"
    patch_nustar_caldb('20170727', caldb=caldb, caldbstage=stage, ck_file_exists=False)
