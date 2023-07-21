"""
See /FTP/caldb/local/software/pycaldb
"""

def update_clockcor(version, file,
                    url='http://www.srl.caltech.edu/NuSTAR_Public/NuSTAROperationSite/Clockfile_v2/',
                    caldb='/FTP/caldb/staging',
                    wwwdir = '/www/htdocs/docs/heasarc/caldb/data/nustar/fpm/index',
                    wwwproddir = '/www.prod/htdocs/docs/heasarc/caldb/data/nustar/fpm/index',
                    templatedir='/Home/home1/corcoran/Python/heasarc/pycaldb/templates'):
    """
    this updates the nustar caldb for a new clock correction file
    @param version: should be of the form YYYYMMDD, i.e. the
    year, month, date of release of the new clock correction file
    which for the clock correction file should be the "valid up to"
    date (for eg, for v43, version = "20150114")
    @param version: date string corresponding to "valid through" date (eg "20150114")
    @param file: name of the new clock correction file (nuCclock20100101v043.fits)
    @param url: url where the clock correction file is located
    @param caldb: location of the local CALDB
    @return: creates a new caldb.indx file with version =  version

    History:
    20200512 MFC: updated url to http://www.srl.caltech.edu/NuSTAR_Public/NuSTAROperationSite/Clockfile_v2/
       to use the fine clock correction implemented at Caltech
    """

    import urllib
    import subprocess
    import os
    from astropy.io import fits as pyfits
    from heasarc.pycaldb.pycaldb import Caldb
    from subprocess import Popen, PIPE
    #import heasoftpy as hsp

    # turn off prompting for batch processing: see
    # https://heasarc.gsfc.nasa.gov/lheasoft/scripting.html
    os.environ['HEADASPROMPT']='/dev/null'

    # file should be of the form nuCclock20100101v043.fits
    nuclockdir=caldb+'/data/nustar/fpm/bcf/clock'
    nucifdir = caldb+'/data/nustar/fpm/index'
    nuclockfile=url+'/'+file
    #nuclockfile_caldb=nuclockdir+'/'+file
    nuclockfile_caldb = os.path.join(caldb,'data/nustar/fpm/bcf/clock',file)
    ###
    # 1a) download the new clock file to the clock directory
    ###
    print ("downloading "+nuclockfile+" to "+nuclockfile_caldb)
    #urllib.urlretrieve(nuclockfile, nuclockfile_caldb)
    try:
        urllib.request.urlretrieve(nuclockfile, filename=nuclockfile_caldb)
    except Exception as e:
        print('Error retrieving {0} ({1}); Stopping'.format(nuclockfile, e))
        return
    ###
    # 1b) ftverify the file
    ###
    print ('FTverifying '+nuclockfile_caldb)
    # returns output of ftverify
    ver = subprocess.check_output('ftverify infile = {0} prstat=no'.format(nuclockfile_caldb), shell=True, universal_newlines=True)
    #print (a[a.find('Verification found'):]) # print summary of Verification results
    #if (a.find('0 warning') and a.find('0 error')) < 0:
    if 'verification OK' not in ver:
        print ('Warning or Error found by ftverify: Stopping')
        return
    print(ver)
    print ('Continuing with caldb ingest')
    newindx='caldb.indx'+version.strip()
    ###
    # 2a) copy the current index file as the new index file
    ###
    print ('Copying '+caldb+'/data/nustar/fpm/caldb.indx to '+nucifdir+'/'+newindx)
    stat = subprocess.call(['cp', caldb+'/data/nustar/fpm/caldb.indx', nucifdir+'/'+newindx])
    if stat != 0:
        print ('Error in copying caldb.indx to '+newindx+': Stopping')
        return
    ###
    # 2b) then run udcif to add the new clock correction file to the new caldb.indx file
    # (remember to set the $CALDB variable as /web_chroot/FTP/caldb if on heasarcdev, but use /FTP on heasarcdev-ol8)
    # (this is a bit kludgy - should have a python version of udcif)
    ###
    print('Updating the cif using udcif')
    os.environ['CALDB']=caldb.strip()
    curdir=os.getcwd()
    os.chdir(nuclockdir)
    print ("Changing directory to "+nuclockdir)
    cmd = 'udcif {file} ../../index/{newindx} quality=0 editc=y '.format(file=file,newindx=newindx)
    print(cmd)
    stat=subprocess.call(['udcif',file,'../../index/'+newindx, 'quality=0', 'editc=y'])
    # if stat != 0:
    #     print('Error running udcif: Stopping')
    #     return
    ###
    # 2d) update the CALDBVER keyword in the caldb.indx file
    ###
    os.chdir(nucifdir) # change directory to cif index directory
    cifhdu = pyfits.open(newindx)
    cifhdu[1].header['CALDBVER'] = version
    print ("Updating CALDBVER keyword to {0} and writing to {1}".format(version, nucifdir+'/'+newindx))
    cifhdu.writeto(nucifdir+'/'+newindx,output_verify='fix', overwrite=True, checksum=True)
    os.chdir(curdir)
    ###
    # 2d) ftverify the  new caldb.indx file
    ###
    print ('\nFTVerifying '+caldb+'/data/nustar/fpm/index/'+newindx)
    a = subprocess.check_output('ftverify '+nucifdir+'/'+newindx, shell=True, universal_newlines=True) # returns output of ftverify as a file object
    print (a[a.find('Verification found'):]) # print summary of Verification results
    if (a.find('0 warning') and a.find('0 error')) < 0:
        print ('Warning or Error found by ftverify: Stopping')
        print ('Stopping')
        return
    ###
    # 3) make a link to the new caldb.indx file
    ###
    print ('Making link from '+nucifdir+'/'+newindx+' to '+caldb+'/data/nustar/fpm/caldb.indx')
    os.chdir(caldb+'/data/nustar/fpm')
    stat=subprocess.call(['ln','-nfs','index/'+newindx,'caldb.indx'])
    if stat != 0:
        print ('Error in making caldb.indx symbolic link: Stopping')
        return
    ###
    # 4) Create the tar files
    ###
    os.chdir(caldb)
    print ('\nChanging directory to '+caldb)
    print ('\nCreating tar files {0}/tmp/goodfiles_nustar_fpm.tar'.format(caldb))
    stat=subprocess.call(['tar','-cvf','tmp/goodfiles_nustar_fpm.tar',
            '--exclude=".*"',
            'data/nustar/fpm/bcf', 
            'data/nustar/fpm/cpf',
            'data/nustar/fpm/index', 
            'data/nustar/fpm/caldb.indx'])
    if stat != 0:
        print ('Error in making tarfile')
    flist=subprocess.check_output(['tar', 'tvf', caldb+'/tmp/goodfiles_nustar_fpm.tar'],universal_newlines=True)
    ascii_file='tmp/goodfiles_nustar_fpm.tar.ascii'
    afile=open(ascii_file,'w')
    afile.write(flist)
    afile.close()
    stat=subprocess.call(['tar','cvf',
            'tmp/goodfiles_nustar_fpm_clockcor.tar',
            '--exclude=".*"',
            'data/nustar/fpm/bcf/clock/{file}'.format(file=file),
            'data/nustar/fpm/index', 
            'data/nustar/fpm/caldb.indx'])
    if stat != 0:
        print ('Error in making clockcor directory tarfile')
    stat = subprocess.call(['gzip',caldb+'/tmp/goodfiles_nustar_fpm.tar'])
    stat = subprocess.call(['gzip', caldb+'/tmp/goodfiles_nustar_fpm_clockcor.tar'])
    stat = subprocess.call(['mv',caldb+'/tmp/goodfiles_nustar_fpm.tar.gz',caldb+'/data/nustar/fpm/goodfiles_nustar_fpm.tar.gz'])
    stat = subprocess.call(['mv',caldb+'/tmp/goodfiles_nustar_fpm.tar.ascii',caldb+'/data/nustar/fpm/goodfiles_nustar_fpm.tar.ascii'])
    stat = subprocess.call(['mv',caldb+'/tmp/goodfiles_nustar_fpm_clockcor.tar.gz',caldb+'/data/nustar/fpm/goodfiles_nustar_fpm_clockcor.tar.gz'])
    ###
    # 5) Create the html versions of the CIF files
    ###
    nucaldb = Caldb(telescope='nustar', instrument='fpm')
    cifvername = "caldb.indx{0}".format(version)
    nucaldb.set_cif(cifvername)
    print("Writing the html CIF (dev version)")
    outdir = wwwdir
    fname = nucaldb.html_summary(missionurl='https://www.nustar.caltech.edu',outdir=outdir, templatedir=templatedir, clobber=True)
    cmd = 'ln -nfs index/{0} {1}/index.html'.format(fname, outdir)
    print(cmd)
    os.system(cmd)
    print("Writing the html CIF (PROD version)")
    outdir = wwwproddir
    fname = nucaldb.html_summary(missionurl='https://www.nustar.caltech.edu', outdir=outdir, templatedir=templatedir,clobber=True)
    cmd = 'ln -nfs index/{0} {1}/index.html'.format(fname, outdir)
    print(cmd)
    os.system(cmd)
    ###
    # 6) rename the tar files with the current version
    ###
    old_tf = Popen('ls /FTP/caldb/data/nustar/fpm/goodfiles_nustar_fpm_2*.tar*', shell=True, stdout=PIPE,  universal_newlines=True).communicate()[0].strip()
    cmd = 'mv {0} /FTP/caldb/data/nustar/fpm/goodfiles_nustar_fpm_{1}.tar.gz'.format(old_tf, version)
    os.system(cmd)
    old_cf = Popen('ls /FTP/caldb/data/nustar/fpm/goodfiles_nustar_fpm_clockcor_2*.tar*', shell=True, stdout=PIPE,  universal_newlines=True).communicate()[0].strip()
    cmd = 'mv {0} /FTP/caldb/data/nustar/fpm/goodfiles_nustar_fpm_clockcor_{1}.tar.gz'.format(old_cf, version)
    os.system(cmd)

    os.chdir(curdir)
    print ("\nDo the following:\n")
    #print "create /docs/heasarc/caldb/data/nustar/fpm/index.html from current caldb.indx file using the Caldb html_summary() method on heasarcdev"
    print ("Create /www/htdocs/docs/heasarc/caldb/nustar/docs/release*.txt file for current release")
    print ("update /www/htdocs/docs/heasarc/caldb/nustar/docs/nustar_caldbhistory.html"          )
    print ("update http://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/nustar/nustar_caldb.html"     )
    print ("edit the supported missions and the what's new page, and the rss feed"               )
    print ("add rss item to /www/htdocs/docs/nustar/news/nustar.rss"                             )
    print ("update astroupdate include file (use RESPONSIBLE_PARTY = NuSTAR for astroupdate.rss)")
    print ("run /www.prod/htdocs/docs/rss/nustar_rss2.pl (and /www/ version)"                    )
    print ("run /www.prod/htdocs/docs/rss/software_rss2.pl and /www/ version"                    )
    return


def get_clockcor(nucchtml='http://www.srl.caltech.edu/NuSTAR_Public/NuSTAROperationSite/clockfile.php'):
    """
    Gets list of all the nustar clock corrections listed on nucchtml,
    returns a dictionary ordered by most recent [0] to oldest [-1] clock correction file
    :param nucchtml: URL where nustar clock corrections are listed
    :return:
    """
    import requests
    from bs4 import BeautifulSoup

    cchtml = requests.get(nucchtml)
    ccsoup = BeautifulSoup(cchtml.text,'lxml')
    li = ccsoup.find_all('li')
    nuccfile = []
    nuccvalidity = []
    for l in li:
        if 'nuCclock' in l.text:
            nuccfile.append(l.text.split('(')[0].strip())
            nuccvalidity.append(l.text.split('(valid up to')[1].strip().replace(')','').replace('-',''))
    nucc = {'clock_correction_file':nuccfile, 'valid_through':nuccvalidity}
    nucdf = pd.DataFrame(nucc)
    nucdf.sort_values(by='clock correction file', inplace=True)
    nucdf.reset_index(inplace=True)
    return nucdf


def check_clockcor_version(ClockCorFile):
    """
    for the specified clockcorversion and clockcorfile check to see that the the clockcorfile is in the CALDB,
    and in the latest caldb.indx file
    :param ClockCorVersion:
    :param ClockCorFile:
    :return:
    """
    # TODO - complete check_clockcor_version
    pass


def sync_web(version, DoCopy=False, user='mcorcora', host="heasarcdev"):
    """
    :param version: version of nustar caldb (YYYYMMDD)
    :param DoCopy: if False just print the cmds that will be executed when DoCopy = True
    :return: 
    """
    import os

    errlist=[]
    curuser = os.environ['LOGNAME']
    curhost = os.environ['HOSTNAME']
    if (curuser != user) or (host not in curhost ):
        error = "This function must be run as {0} from {1}; returning".format(user, host)
        print(error)
        print("Current user is:{0}  Current host is: {1}".format(curuser, curhost))
        errlist.append(error)
        return errlist
    files_to_sync = [
        ("caldb_wotsnew.html","/www/htdocs/docs/heasarc/caldb","/www.prod/htdocs/docs/heasarc/caldb"),
        ("caldb_supported_missions.html", "/www/htdocs/docs/heasarc/caldb", "/www.prod/htdocs/docs/heasarc/caldb"),
        ("caldb.rss", "/www/htdocs/docs/heasarc/caldb", "/www.prod/htdocs/docs/heasarc/caldb"),
        ("release_{0}.txt".format(version), "/www/htdocs/docs/heasarc/caldb/nustar/docs/","/www.prod/htdocs/docs/heasarc/caldb/nustar/docs/"),
        ("nustar_caldbhistory.html","/www/htdocs/docs/heasarc/caldb/nustar/docs/","/www.prod/htdocs/docs/heasarc/caldb/nustar/docs/"),
        ("nustar_caldb.html","/www/htdocs/docs/heasarc/caldb/nustar","/www.prod/htdocs/docs/heasarc/caldb/nustar"),
        ("cif_nustar_fpm_{0}.html".format(version),"/FTP/caldb/local/software/pro/DATA_DELIVERIES/pro/mission_summaries/work","/www/htdocs/docs/heasarc/caldb/data/nustar/fpm/index"),
        ("cif_nustar_fpm_{0}.html".format(version), "/FTP/caldb/local/software/pro/DATA_DELIVERIES/pro/mission_summaries/work","/www.prod/htdocs/docs/heasarc/caldb/data/nustar/fpm/index"),
        ("nustar.rss","/www/htdocs/docs/nustar/news/","/www.prod/htdocs/docs/nustar/news/"),
        ("associated_data.html","/www/htdocs/docs/heasarc/astro-update/inc","/www.prod/htdocs/docs/heasarc/astro-update/inc"),
        ("astro-update.rss","/www/htdocs/docs/heasarc/astro-update/","/www.prod/htdocs/docs/heasarc/astro-update/")
    ]
    print ("Syncing website:" )
    for f in files_to_sync:
        print ("copying {0} from {1} to {2}".format(f[0],f[1],f[1]) )
        cmd = "cp {1}/{0} {2}/{0}".format(f[0],f[1],f[2])
        print ("{0}\n".format(cmd) )
        if DoCopy:
            status = os.system(cmd)
            if status != 0:
                errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    print ("Creating nustar and heasarc newsfeeds (dev version)" )
    cmd = "/www/htdocs/docs/rss/nustar_rss2.pl"
    print("{0}\n".format(cmd))
    if DoCopy:
        status = os.system(cmd)
        if status  !=  0:
            errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    print ("Creating nustar and heasarc newsfeeds (dev version)")
    cmd = "/www/htdocs/docs/rss/nustar_rss2.pl"
    print("{0}\n".format(cmd))
    if DoCopy:
        status = os.system(cmd)
        if status  !=  0:
            errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    return errlist

if __name__ == "__main__":
    import os
    #caldb = '/web_chroot/FTP/caldb' # appropriate for running as caldbmgr on heasarcdev
    # nucc = nustar_get_clockcor()
    # nucckeys = nucc.keys()
    # for i in range(len(nucc[nucckeys[0]])):
    #     print nucc[nucckeys[0]][i], nucc[nucckeys[1]][i]

    # os.system('rm /Users/mcorcora/software/caldb/data/nustar/fpm/index/caldb.indx20200428')
    caldb =  '/Users/mcorcora/software/caldb'
    ccurl = 'http://www.srl.caltech.edu/NuSTAR_Public/NuSTAROperationSite/Clockfile_v2/nuCclock20100101v101.fits.gz'
    ncurl = os.path.split(ccurl)[0]
    cfile = os.path.split(ccurl)[1]
    templatedir = '/Users/mcorcora/software/github/heasarc/pycaldb/templates'
    update_clockcor('20200510', cfile, url=ncurl, caldb=caldb,
                       templatedir=templatedir)

    #update_clockcor('20180126', 'nuCclock20100101v078.fits', caldb=caldb)
    #update_clockcor('20171204', 'nuCclock20100101v076.fits', caldb=caldb)
    #update_clockcor('20171002', 'nuCclock20100101v075.fits', caldb=caldb)
    #update_clockcor('20170817', 'nuCclock20100101v074.fits', caldb=caldb)
    #update_clockcor('20170720', 'nuCclock20100101v073.fits', caldb=caldb)
    #update_clockcor('20170614', 'nuCclock20100101v072.fits', caldb=caldb)
    #update_clockcor('20170503', 'nuCclock20100101v071.fits', caldb=caldb)
    #update_clockcor('20170419','nuCclock20100101v070.fits',caldb=caldb)
    #update_clockcor('20170222','nuCclock20100101v069.fits',caldb=caldb)
    #update_clockcor('2017-01-20','nuCclock20100101v068.fits',caldb=caldb)
    #update_clockcor('20161207','nuCclock20100101v066.fits', caldb=caldb)
    #update_clockcor('20160922','nuCclock20100101v062.fits', caldb=caldb)
    #update_clockcor('20160731','nuCclock20100101v060.fits', caldb=caldb)
    #update_clockcor('20151008', 'nuCclock20100101v052.fits', caldb=caldb)
    #update_clockcor('20150904', 'nuCclock20100101v051.fits', caldb=caldb)
    #update_clockcor('20150114', 'nuCclock20100101v043.fits', caldb='/fuse/caldb_staging/data/nustar/versions/20150114')
    #update_clockcor('20150114','nuCclock20100101v043.fits',caldb='/Volumes/USRA16/caldb')
