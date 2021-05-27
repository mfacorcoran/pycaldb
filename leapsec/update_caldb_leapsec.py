def update_caldb_leapsec(NewLeapSecondDate, NewLeapSecond, updater="MFC", outdir=".", clobber=True, leapsecfile=None):
    """
    Given the date of a new leapsecond (for example 2017-01-01T00:00:00)
    creates a new leapsecond file for transfer to the caldb

    NewLeapSecDate = date of new leap second in ISOT format YYYY-MM-DDTHH:MM:SS
    NewLeapSecond = amount of new leap second (usually 1.0)

    writes a new FITS file to the current working directory by default
    """
    from astropy.io import fits as pyfits
    from astropy.time import Time
    from astropy.table import Table
    import ftputil
    import time
    if leaspsecfile is not None:
        #
        # this block retrieves the lastest leapsecond file from /FTP/caldb/data/gen/bcf
        # based on the leapsecond file naming convention, "leapsec_mmddyy.fits"
        #
        LSdir = "FTP/caldb/data/gen/bcf/"
        host = ftputil.FTPHost('heasarc.gsfc.nasa.gov', "anonymous", "mcorcoran@usra.edu")
        genbcf = host.listdir(LSdir)  # get directory listing
        host.close()

        LeapsecFileList = [f for f in genbcf if 'leapsec' in f]
        LeapsecFileYear = [y.split("_")[1].split(".fits")[0][4:6] for y in LeapsecFileList]
        LeapsecFileYear = [('19' + y if int(y) > 50 else '20' + y) for y in LeapsecFileYear]
        LeapsecFileMonth = [m.split("_")[1].split(".fits")[0][2:4] for m in LeapsecFileList]
        LeapsecFileDay = [d.split("_")[1].split(".fits")[0][0:2] for d in LeapsecFileList]

        maxjd = 0.0

        for i in range(len(LeapsecFileList)):
            tiso = LeapsecFileYear[i] + "-" + LeapsecFileMonth[i] + "-" + LeapsecFileDay[i]
            fjd = Time(tiso).jd
            if fjd > maxjd:
                maxjd = fjd
                LatestLSF = LeapsecFileList[i]
        print("Latest Leapsecond File = {0}".format(LatestLSF))
    else:
        LSdir = os.path.split(leapsecfile)[0]
        LatestLSF = os.path.split(leapsecfile)[-1]

    hdu = pyfits.open("http://heasarc.gsfc.nasa.gov/" + LSdir + "/" + LatestLSF) # open the file

    orig_header = hdu[1].header
    mjdref = orig_header['MJDREF']

    UpdateDate = time.strftime('%Y-%m-%d %H:%M:%S')
    outfile = 'leapsec_' + NewLeapSecondDate[8:10] + NewLeapSecondDate[5:7] + NewLeapSecondDate[2:4] + '.fits'


    #
    # create new row for new leapsecond information
    #
    newdate = NewLeapSecondDate.split('T')[0]
    newtime = NewLeapSecondDate.split('T')[1]
    newmjd = Time(NewLeapSecondDate, format='isot').mjd
    newsecs = (newmjd - mjdref) * 86400
    # add the leapseconds to newsecs since each leapsecond day is 86400+N*leapseconds
    # where N is the row number starting at 0 for the first row
    nPreviousLeapSecs = len(hdu[1].data)
    newsecs = newsecs + nPreviousLeapSecs
    newLS = NewLeapSecond
    NewLeapSecondRow = [newdate, newtime, newmjd, newsecs, newLS]
    #
    # append new leapsecond to table and write to output file
    #
    tbdata = hdu[1].data
    t = Table(tbdata)  # convert hdu data to a python Table to add new row
    t.add_row(NewLeapSecondRow)  # add row of data
    hdunew = pyfits.table_to_hdu(t)  # convert table back to hdu (with minimal header)

    hdunew.columns.change_unit('SECONDS', 's') # table_to_hdu doesn't seem to preserve the Unit
    hdunew.columns.change_unit('LEAPSECS', 's') # table_to_hdu doesn't seem to preserve the Unit

    hdunew.header = orig_header  # use header from original file
    hdunew.header['COMMENT'] = UpdateDate+": "+updater+" ADDED "+NewLeapSecondDate+" LEAP SECOND"
    hdunew.header['HISTORY'] = "File modified by user "+updater+" on "+UpdateDate
    pyfits.writeto(outdir + "/" + outfile, hdunew.data, hdunew.header, clobber=clobber, checksum=True)
    return outfile

def mk_leapsectab(leaplist = 'https://hpiers.obspm.fr/iers/bul/bulc/ntp/leap-seconds.list',
                  mjdref=40587.0,
                  outname=None,
                  clobber=False):
    """
    creates the leap second table from the list at leaplist, and optionally writes it to a fits file
    :param leaplist: URL of the location of the list of leapseconds
    :param mjdref: reference time for output SECONDS (default corresponds to 1970-01-01 00:00:00.000)
    :param outname: output name of FITS file (including directory path)
    :return: leapsectab
    """
    __version__='1.0'
    __name__ = 'mk_leapsectab'
    __author__="M. F. Corcoran"
    from astropy.time import Time
    from astropy.table import Table
    import requests
    import numpy as np

    # get the list of leapseconds
    if 'https://' in leaplist:
        req = requests.get(leaplist)
        lsfile=req.text
        # parse the leapsecond list
        lines=lsfile.split("\n")
    else:
        with open(leaplist,'r') as f:
            lines = f.readlines()
            lines = [l.strip('\n') for l in lines]

    leapsecs={'DATE':[],'TIME':[],'MJD':[],'SECONDS':[],'LEAPSECS':[],'NTPTIME':[],'DTAI':[]}
    for l in lines:
        #print(l[0])
        if len(l)>0:
            if "#" != l[0]:
                leapsecs['NTPTIME'].append(float(l[:15].strip()))
                leapsecs['DTAI'].append(float(l[15:20].strip()))
                #leapsecs['DAYMONYEAR'].append(l[25:].strip())
    # make the leapsecs dictionary
    # mjd 15020 = 1900.0
    leapsecs['MJD'] = np.asarray(leapsecs['NTPTIME'])/86400+15020
    leapsecs['SECONDS'] = (np.asarray(leapsecs['MJD'])-mjdref)*86400+np.asarray(leapsecs['DTAI'])-10
    leapsecs['DATE'] = [x.split()[0] for x in Time(leapsecs['MJD'],format='mjd').iso]
    leapsecs['TIME'] = [x.split()[1] for x in Time(leapsecs['MJD'],format='mjd').iso]
    cumtime = np.insert(np.arange(len(leapsecs['MJD'])-1),0,0)
    leapsecs['LEAPSECS'] = np.asarray(leapsecs['DTAI'])-10 - cumtime
    # create table
    leapsectab = Table(leapsecs)
    leapsectab['SECONDS'].unit='s'
    leapsectab['LEAPSECS'].unit = 's'
    leapsectab['NTPTIME'].unit = 's'
    leapsectab['DTAI'].unit = 's'
    # update metadata
    leapsectab.meta['EXTNAME'] = "LEAPSECONDS"
    leapsectab.meta['TELESCOP'] = "GEN"
    leapsectab.meta['INSTRUME'] = "INS"
    leapsectab.meta['CREATOR'] = __name__
    leapsectab.meta['DATE'] = Time.now().fits
    leapsectab.meta['CCLS0001']= 'BCF'
    leapsectab.meta['CCNM0001']= 'LEAPSECS'
    leapsectab.meta['CDTP0001']= 'DATA'
    leapsectab.meta['CVSD0001']= '1970-01-01'
    leapsectab.meta['CVST0001']= '00:00:00'
    for i in np.arange(1,10):
        leapsectab.meta['CBD{i}0000'.format(i=i)] = 'NONE'
    leapsectab.meta['CDES0001']= 'Table of times at which Leap seconds occurred'
    leapsectab.meta['TIMESYS']='1970-01-01T00:00:00'
    leapsectab.meta['MJDREF']= mjdref
    leapsectab.meta['REFYEAR']= 1970.0
    trunclines = [x[0:79].replace('\t','').replace('#',' ') for x in lines]
    trunclines.append('File {0}'.format(leaplist))
    coms = [' ','From {0}'.format(leaplist),' ']
    coms.extend(trunclines)
    leapsectab.meta['comments'] =coms
    leapsectab.meta['history'] = 'Last update: {0} (UT)'.format(Time.now().iso)
    if outname is not None:
        # write fits file to outname
        leapsectab.write(outname, format='fits', overwrite=clobber)
    return leapsectab

def cmp_leaplist(current='https://hpiers.obspm.fr/iers/bul/bulc/ntp/leap-seconds.list',
                 archived = 'https://heasarc.gsfc.nasa.gov/FTP/caldb/data/gen/bcf/leapsec.fits'):
    """
    compares the IERS and CALDB leapsecond tables
    :param current: IERS current file
    :param archived: CALDB leapsecond file
    :return: True if the tables are the same, False if not
    """
    from astropy.table import Table
    curtab = mk_leapsectab(leaplist=current)
    #print('reading {0}'.format(archived))
    arctab = Table.read(archived, hdu='LEAPSECONDS')
    # test is a logical array = all True if all the elements are the same
    # minimum of a logical array containing all Trues is True; False otherwise
    testcols = ['LEAPSECS', 'MJD', 'NTPTIME']
    test = curtab[testcols] == arctab[testcols]
    #print('Are the files the same? {0}'.format(test.min()))
    return test.min()

def check_leaplist(iers='https://hpiers.obspm.fr/iers/bul/bulc/ntp/leap-seconds.list',
                 caldb = 'https://heasarc.gsfc.nasa.gov/FTP/caldb/data/gen/bcf/leapsec.fits',
                   verbose = False, testing=False):
    """
    checks that the IERS and CALDB leapseconds file have the same entries, send e-mail notification
    to caldbmgr@bigbang.gsfc.nasa.gov if not
    :param iers: IERS leapsecond file
    :param caldb: CALDB leapsecond fits file
    :return:
    """
    from astropy.table import Table
    import os
    same = cmp_leaplist(current=iers, archived=caldb)
    if testing:
        same=False
    if not same:
        if verbose:
            print('IERS and CALDB file have DIFFERENT LEAPSECOND entries')
        arctab = Table.read(caldb, hdu='LEAPSECS')
        curtab = mk_leapsectab(leaplist=iers)
        with open('/tmp/leapsec_notification.txt', 'w') as f:
            f.write('CALDB Last Leapsecond on {0:.5f}\n'.format(arctab['MJD'].max()))
            f.write('IERS  Last Leapsecond on {0:.5f}\n'.format(curtab['MJD'].max()))
        cmd = 'mail -s "IERS LEAPSECONDS TABLE UPDATED" caldbmgr@bigbang.gsfc.nasa.gov < /tmp/leapsec_notification.txt'
        if verbose:
            print(cmd)
        stat = os.system(cmd)
    else:
        stat=0
        if verbose:
            print('IERS and CALDB file have the same LEAPSECOND entries')
    return stat


if __name__ == "__main__":
    file = update_caldb_leapsec('2017-01-01T00:00:00', 1.0)
    print(file)