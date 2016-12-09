def update_caldb_leapsec(NewLeapSecondDate, NewLeapSecond, updater="MFC", outdir=".", clobber=True):
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

if __name__ == "__main__":
    file = update_caldb_leapsec('2017-01-01T00:00:00', 1.0)
    print file