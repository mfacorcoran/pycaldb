def chandra_caldb_update(version, CaldbWorkDir = "/FTP/caldb/staging/cxc",
                         ChandraCaldbDir = "/pub/arcftp/caldb",
                         ChandraCaldbHost = "cda.harvard.edu",
                         user="caldbmgr", host="heasarcdev"):
    """
    This function downloads the chandra caldb from the CXC site to the heasarc caldb then untars it and installs 
    the HEASARC version in the HEASARC CALDB
    
    version should be a string like 
        version = "4.7.2"
    """

    import tarfile
    import os
    import shutil
    import requests
    from bs4 import BeautifulSoup
    import urllib

    errlist=[]
    curuser = os.environ['LOGNAME']
    curhost = os.environ['HOSTNAME']
    if (curuser<>user) or (host not in curhost ):
        error = "This function must be run as {0} from {1}; returning".format(user, host)
        print(error)
        print("Current user is:{0}  Current host is: {1}".format(curuser, curhost))
        errlist.append(error)
        return errlist

    # get the urls of the files to be downloaded from the FULL caldb installation table
    req = requests.get('http://cxc.harvard.edu/ciao/download/caldb.html')
    soup = BeautifulSoup(req.text, 'lxml')
    tab = soup.find_all('table')[0] # first table on this page is the FULL install table
    rows = tab.find_all('tr')
    files_to_download=[]
    for r in rows[1:]: # exclude first row since it's the header
        f = r.find_all('td')[1].a.get('href')
        files_to_download.append(f)
    print('Found these caldb files to download from CXC')
    for f in files_to_download:
        print(f)
    ans = raw_input('Continue (y/N)> ')
    if ans.strip()[0].lower() <> 'y':
        print "User cancelled execution; Returning"
        return

    DownLoadDir=CaldbWorkDir+"/chandra-"+version.strip()
    if not os.path.exists(DownLoadDir):
        os.makedirs(DownLoadDir)
    filename = [f.strip()[f.strip().rfind('/')+1:] for f in files_to_download] # names of all files in FULL download table from CXC

    # retrieve the files to the download directory using wget
    os.chdir(DownLoadDir)
    for f in files_to_download:
        cmd = 'wget {0}'.format(f)
        print cmd
        status = os.system(cmd)
        if status <> 0:
            ans = raw_input('Problem with wget; Continue (y/N)> ')
            if ans.strip()[0].lower() <> 'y':
                print "User cancelled execution; Returning"
                return
    #
    # untar the tar files and delete the tar file after untarring
    #
    os.chdir(DownLoadDir)
    for t in filename:
        if '.tar' in t:
            print "Untarring {0}".format(t)
            tarf = tarfile.open(t)
            try:
                tarf.extractall()
            except:
                print "ERROR UNTARRING {0}; returning".format(t)
                status = -99
                return status
            print "Deleting tar file {0}".format(t)
            os.remove(t)
    
    #
    # don't forget the manifest and readme files
    #
    manfile = 'MANIFEST_{0}_main.txt'.format(version.strip())
    print "Moving {1} into {0}/docs/chandra/manifests{1}".format(DownLoadDir,manfile)
    localfile = DownLoadDir+"/docs/chandra/manifests/"+manfile
    shutil.move("{0}/{1}".format(DownLoadDir,manfile), "{0}/docs/chandra/manifests/{1}".format(DownLoadDir,manfile))

    readme = "README_caldb{0}.txt".format(version.strip())
    print "Moving {1} into {0}/docs/chandra/{1}".format(DownLoadDir,readme)
    shutil.move("{0}/{1}".format(DownLoadDir,readme), "{0}/docs/chandra/{1}".format(DownLoadDir,readme))
    
    os.chdir(CaldbWorkDir)
    
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
 