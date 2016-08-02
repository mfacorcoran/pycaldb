"""
See /FTP/caldb/local/software/pycaldb
"""

def update_clockcor(version, file,
                    url='http://www.srl.caltech.edu/NuSTAR_Public/NuSTAROperationSite/Clockfile/',
                    caldb='/FTP/caldb'):
    """
    this updates the nustar caldb for a new clock correction file
    @param version: should be of the form YYYYMMDD, i.e. the
    year, month, date of release of the new clock correction file
    which for the clock correction file should be the "valid up to"
    date (for eg, for v43, version = "20150114")
    @param file: name of the new clock correction file (nuCclock20100101v043.fits)
    @param url: url where the clock correction file is located
    @param caldb: location of the local CALDB
    @return:
    """
    import urllib
    import subprocess
    import os
    from astropy.io import fits as pyfits
    # file should be of the form nuCclock20100101v043.fits
    nuclockdir=caldb+'/data/nustar/fpm/bcf/clock'
    nucifdir = caldb+'/data/nustar/fpm/index'
    nuclockfile=url+'/'+file
    nuclockfile_caldb=nuclockdir+'/'+file
    """
    1a) download the new clock file to the clock directory
    """
    print "downloading "+nuclockfile+" to "+nuclockfile_caldb
    urllib.urlretrieve(nuclockfile, nuclockfile_caldb)
    """
    1b) ftverify the file
    """
    print 'FTverifying '+nuclockfile_caldb
    a = subprocess.check_output('ftverify '+nuclockfile_caldb, shell=True) # returns output of ftverify as a file object
    print a[a.find('Verification found'):] # print summary of Verification results
    if (a.find('0 warning') and a.find('0 error')) < 0:
        print 'Warning or Error found by ftverify: Stopping'
        return
    print 'Continuing with caldb ingest'
    newindx='caldb.indx'+version.strip()
    """
    2a) copy the current index file as the new index file
    """
    print 'Copying '+caldb+'/data/nustar/fpm/caldb.indx to '+nucifdir+'/'+newindx
    stat = subprocess.call(['cp', caldb+'/data/nustar/fpm/caldb.indx', nucifdir+'/'+newindx])
    if stat<>0:
        print 'Error in copying caldb.indx to '+newindx+': Stopping'
        return
    """
    2b) then run udcif to add the new clock correction file to the new caldb.indx file
    (remember to set the $CALDB variable as /web_chroot/FTP/caldb if on heasarcdev)
    (this is a bit kludgy - should have a python version of udcif)
    """
    os.environ['CALDB']=caldb.strip()
    curdir=os.getcwd()
    os.chdir(nuclockdir)
    print "Changing directory to "+nuclockdir
    stat=subprocess.call(['udcif',file,'../../index/'+newindx, 'quality=0', 'editc=y'])
    """
    2d) update the CALDBVER keyword in the caldb.indx file
    """
    os.chdir(nucifdir) # change directory to cif index directory
    cifhdu = pyfits.open(newindx)
    cifhdu[1].header['CALDBVER'] = version
    print "Updating CALDBVER keyword to {0} and writing to {1}".format(version, nucifdir+'/'+newindx)
    cifhdu.writeto(nucifdir+'/'+newindx,output_verify='fix', clobber=True, checksum=True)
    os.chdir(curdir)
    """
    2d) ftverify the  new caldb.indx file
    """
    print '\nFTVerifying '+caldb+'/data/nustar/fpm/index/'+newindx
    a = subprocess.check_output('ftverify '+nucifdir+'/'+newindx, shell=True) # returns output of ftverify as a file object
    print a[a.find('Verification found'):] # print summary of Verification results
    if (a.find('0 warning') and a.find('0 error')) < 0:
        print 'Warning or Error found by ftverify: Stopping'
        print 'Stopping'
        return
    """
    3) make a link to the new caldb.indx file
    """
    print 'Making link from '+nucifdir+'/'+newindx+' to '+caldb+'/data/nustar/fpm/caldb.indx'
    os.chdir(caldb+'/data/nustar/fpm')
    stat=subprocess.call(['ln','-nfs','index/'+newindx,'caldb.indx'])
    if stat<>0:
        print 'Error in making caldb.indx symbolic link: Stopping'
        return
    """
    4) Create the tar files
    """
    os.chdir(caldb)
    print '\nChanging directory to '+caldb
    print '\nCreating tar files'
    stat=subprocess.call(['tar','-cvf','tmp/goodfiles_nustar_fpm.tar', 'data/nustar/fpm/bcf', 'data/nustar/fpm/cpf',
                          'data/nustar/fpm/index', 'data/nustar/fpm/caldb.indx'])
    if stat<>0:
        print 'Error in making tarfile'
    flist=subprocess.check_output(['tar', 'tvf', caldb+'/tmp/goodfiles_nustar_fpm.tar'])
    ascii_file='tmp/goodfiles_nustar_fpm.tar.ascii'
    afile=open(ascii_file,'w')
    afile.write(flist)
    afile.close()
    stat=subprocess.call(['tar','cvf', 'tmp/goodfiles_nustar_fpm_clockcor.tar','data/nustar/fpm/bcf/clock',
                          'data/nustar/fpm/index', 'data/nustar/fpm/caldb.indx'])
    if stat<>0:
        print 'Error in making clockcor directory tarfile'
    stat = subprocess.call(['gzip',caldb+'/tmp/goodfiles_nustar_fpm.tar'])
    stat = subprocess.call(['gzip', caldb+'/tmp/goodfiles_nustar_fpm_clockcor.tar'])
    stat=subprocess.call(['mv',caldb+'/tmp/goodfiles_nustar_fpm.tar.gz',caldb+'/data/nustar/fpm/goodfiles_nustar_fpm.tar.gz'])
    stat=subprocess.call(['mv',caldb+'/tmp/goodfiles_nustar_fpm.tar.ascii',caldb+'/data/nustar/fpm/goodfiles_nustar_fpm.tar.ascii'])
    stat=subprocess.call(['mv',caldb+'/tmp/goodfiles_nustar_fpm_clockcor.tar.gz',caldb+'/data/nustar/fpm/goodfiles_nustar_fpm_clockcor.tar.gz'])

    os.chdir(curdir)

if __name__ == "__main__":
    caldb = '/web_chroot/FTP/caldb' # appropriate for running as caldbmgr on heasarcdev
    update_clockcor('20151008', 'nuCclock20100101v052.fits', caldb=caldb)
    #update_clockcor('20150904', 'nuCclock20100101v051.fits', caldb=caldb)
    #update_clockcor('20150114', 'nuCclock20100101v043.fits', caldb='/fuse/caldb_staging/data/nustar/versions/20150114')
    #update_clockcor('20150114','nuCclock20100101v043.fits',caldb='/Volumes/USRA16/caldb')
