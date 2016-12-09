def chandra_caldb_update(version, CaldbWorkDir = "/FTP/caldb/staging/cxc", ChandraCaldbDir = "/pub/arcftp/caldb",
    ChandraCaldbHost = "cda.harvard.edu"):
    """
    This function downloads the chandra caldb from the CXC site to the heasarc caldb then untars it and installs 
    the HEASARC version in the HEASARC CALDB
    
    version should be a string like 
        version = "4.7.2"
    """
    
    from astropy.io import fits as pyfits
    from astropy.time import Time
    from astropy.table import Table
    import ftputil
    import time
    import tarfile
    import os
    import shutil
    
    print "Connecting to FTP site\n"
    
    status = 0
    
    host = ftputil.FTPHost(ChandraCaldbHost, "anonymous", "mcorcoran@usra.edu")
    ChandraCaldbFiles = host.listdir(ChandraCaldbDir) # get list of caldb files
    
    ChandraCaldbTarFiles = [f for f in ChandraCaldbFiles if '.tar' in f ] # get list of tar files
    print"Found these Tarfiles:"
    print ChandraCaldbTarFiles
    
    CaldbMainTar = "caldb_"+version.strip()+"_main.tar.gz"
    print CaldbMainTar
 
    try: 
        ChandraCaldbTarFiles.remove(CaldbMainTar) # want caldb main file to be first downloaded and untarred
    except ValueError:
        print "{0} not in CDA FTP directory {1}; returning".format(CaldbMainTar, ChandraCaldbDir)
        status = -99
        host.close()
        return status
    ChandraCaldbTarFiles = [CaldbMainTar]+ChandraCaldbTarFiles
   
    print "ChandraCaldbTarfiles"
    print ChandraCaldbTarFiles
        
    DownLoadDir=CaldbWorkDir+"/chandra-"+version.strip()
    if os.path.isdir(DownLoadDir): # if the directory already exists remove it
        shutil.rmtree(DownLoadDir)
    os.mkdir(DownLoadDir)
 
    print "Changing directory to {0}\n".format(DownLoadDir)
    os.chdir(DownLoadDir)

    print ("\nFile List:")
    print (ChandraCaldbTarFiles)
    for f in ChandraCaldbTarFiles:
        remotefile = ChandraCaldbDir+"/"+f
        localfile = DownLoadDir+"/"+f
        if os.path.exists(f):
            print ("File {0} already exists; deleting prior to download".format(f))
            os.remove(f)
        print "\nGetting file "+f
        getfile = True
        inum = 0 # download counter 
        itrymax = 3 # maximum number of download attempts
        while (getfile and inum<itrymax): # download each time; try a max of itrymax times
            try:
                getfile=False
                #host.download(remotefile,localfile,mode='b')
                host.download(remotefile,localfile)
            except IOError:
                print "IOError on download of file {0}".format(remotefile)
                getfile = True
                inum =+ 1
            if inum >= itrymax: # file did not download so return with error
                status = -99
                print "Could not retrieve file {0} in {1} attempts; returning".format(f, itrymax)
                host.close()
                return status
            else:
                #
                # once file downloaded, then untar it and delete the tar file
                #
                print "Untarring {0}".format(f)
                tarf = tarfile.open(f)
                try:
                    tarf.extractall()
                except:
                    print "ERROR UNTARRING {0}; returning".format(localfile)
                    status = -99
                    return status        
                print "Deleting tar file {0}".format(f)
                os.remove(f)
    
    #
    # don't forget the manifest and readme files
    #
    manfile = [f for f in ChandraCaldbFiles if 'MANIFEST' in f]
    manfile = manfile[0]
    remotefile = ChandraCaldbDir+"/"+manfile
    localfile = DownLoadDir+"/docs/chandra/manifests/"+manfile
    try:
        host.download(remotefile, localfile)
    except IOError:
        print "\nError retrieving {0}; continuing...".format(manfile)
    
    readme = [f for f in ChandraCaldbFiles if 'README' in f]
    readme = readme[0]
    remotefile = ChandraCaldbDir+"/"+readme
    localfile = DownLoadDir+"/docs/chandra/"+readme
    try:
        host.download(remotefile, localfile)
    except IOError:
        print "\nError retrieving {0}; continuing...".format(manfile)

    
    os.chdir(CaldbWorkDir)
    
    host.close()
    
    #
    # move to public FTP side
    #
    staging = "/FTP/caldb/staging/cxc/chandra-"+version.strip()
    public = "/FTP/caldb/data/chandra-"+version.strip()
    if os.path.isdir(public):
        ans = raw_input("\n"+public+' ALREADY EXISTS.  Remove it (Y/n) ?')
        if not ans.lower().strip()[0] == 'n':
            print "Removing {0}".format(public)
            shutil.rmtree(public)
        else: 
            print "Directory not deleted; returning"
            status = -99
            return status

    shutil.move(staging, public)
    #
    # make symlink
    #
    src = "/FTP/caldb/data/chandra-"+version.strip()+"/data/chandra"
    dst = "/FTP/caldb/data/chandra"
    if os.path.islink(dst):
        os.remove(dst) 
    print "Creating symbolic link from {0} to {1}".format(src, dst)
    os.symlink(src, dst)
    
    return status
 