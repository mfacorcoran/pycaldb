def patch_nustar_caldb(version,
                       nupatchserver = "hassif.caltech.edu",
                       nupatchworkdir = "/pub/nustar/pipeline",
                       caldb = "/FTP/caldb",
                       caldbstage = "/FTP/caldb/staging",
                       ck_file_exists='True'):
    """
    CALLING SEQUENCE:

    patch_nustar_caldb(version,
                       nupatchserver = "hassif.caltech.edu",
                       nupatchworkdir = "/pub/nustar/pipeline",
                       caldb = "/FTP/caldb",
                       caldbstage = "/FTP/caldb/staging")

    The following is from /web_chroot/FTP/.caldb_staging/data/nustar/README_rev1.txt

    20131028
    ========

    installing patch version 20131007 - STEPS:

    a) first created a full version of the existing nustar caldb in the
    /FTP/caldb/data/nustar in the nustar/versions are via rsync:

       [caldbmgr]: pwd
       /web_chroot/FTP/.caldb_staging/data/nustar/versions
       [caldbmgr]: mkdir 20131007
       [caldbmgr]: cd 20131007
       [caldbmgr]: rsync -avz 20130509/CALDB/data/nustar/fpm/bcf .
       [caldbmgr]: rsync -avz 20130509/CALDB/data/nustar/fpm/cpf .
       [caldbmgr]: rsync -avz 20130509/CALDB/data/nustar/fpm/index .

    b) wget the patch (using the location given in Brian Grefenstette's e-mail):

       [caldbmgr]: cd /FTP/caldb/staging/data/nustar/versions/20131007/
       [caldbmgr]: wget ftp://hassif.srl.caltech.edu/pub/nustar/pipeline/20131007/NUSTARCALDB_patch20131007.tgz


    c) Untar the patch:

      [caldbmgr]: pwd
      /web_chroot/FTP/.caldb_staging/data/nustar/versions/20131007

       [caldbmgr]: ls
       fpm  NUSTARCALDB_patch20131007.tgz
       [caldbmgr]: tar -zxvf NUSTARCALDB_patch20131007.tgz
         fpm/caldb.indx
         fpm/index/caldb.indx20131007
         fpm/bcf/arf/nuA20100101v006.arf
         fpm/bcf/arf/nuB20100101v006.arf
         ... stuff deleted ...
         fpm/bcf/instrument/nuAhkrange05_20100101v005.fits
         fpm/bcf/instrument/nuAhkrange04_20100101v005.fits
         fpm/bcf/instrument/nuAhkrange03_20100101v005.fits
         fpm/bcf/instrument/nuAhkrange02_20100101v005.fits

    c) then link the (now full) 20131007 version to staging/data/nustar/fpm:

       [caldbmgr]: pwd
            /web_chroot/FTP/.caldb_staging/data/nustar
       [caldbmgr]: ln -nfs versions/20131007/fpm .
       [caldbmgr]: ls -l
            total 12
            lrwxrwxrwx 1 caldbmgr caldb         21 Oct 28 15:31 fpm -> versions/20131007/fpm
            -rw-r--r-- 1 caldbmgr caldb       1909 Oct 28 15:27 README.txt
            drwxr-xr-x 2 caldbmgr caldb       4096 Oct 28 15:09 to_ingest
            drwxr-xr-x 7 caldbmgr caldb_stage 4096 Oct 28 15:14 versions

    d) from Brian's e-mail list of new files in the patch, created caldb_update_nustar_20131007.txt in the "to_ingest" subdirectory:

         [caldbmgr]: pwd
             /web_chroot/FTP/.caldb_staging/data/nustar/to_ingest

         [caldbmgr]: ls -l
             total 16
             -rw-r--r-- 1 caldbmgr caldb 8839 Aug 14 11:27 caldb_update_nustar_20130509.txt
             -rw-r--r-- 1 caldbmgr caldb 1636 Oct 28 15:09 caldb_update_nustar_20131007.txt
             lrwxrwxrwx 1 caldbmgr caldb   32 Aug 14 11:02 caldb_update.txt -> caldb_update_nustar_20130509.txt

    e) to ingest the data (as in other missions) run /web_chroot/FTP/caldb/local/scripts/DATA_DELIVERIES/nustar_caldbupdate.csh

    :param version: string in 'YYYYMMDD' format, from nustar patch notification (from BG)
    :param nupatchworkdir: URL of directory where patch tar file is located, from BG
    :param caldb: location where the full NuSTAR caldb is stored
    :param caldbstage: location of caldb staging area
    :return:
    """

    import os
    import ftputil
    import tarfile
    import datetime
    import shutil
    from pycaldb import ck_file_existence
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
    newfiles = ptf.getnames()
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
        name=f.split('fpm/')[1]
        if not 'caldb.indx' in name.strip().lower(): # don't include caldb.indx files in ingest note
            update_txt=update_txt+name+"\n"
    update_txt=update_txt+'endnewfile\n'
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

    return



if __name__ == "__main__":
    #
    # for sshfs mounted directories
    #
    caldb="/fuse/caldb"
    stage = '/fuse/caldb_stage'
    patch_nustar_caldb('20160606', caldb=caldb, caldbstage=stage, ck_file_exists=False)
