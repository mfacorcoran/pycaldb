def cif_summary(cif,max_cal_qual=4):
    """
    reads a cif and returns some interesting info
    
    KEYWORDS
    max_cal_qual - maximum value of calibration data quality.  default is to only consider good quality data (cal_qual <= 4);
    
    """
    from astropy.io import fits as pyfits
    import numpy as np
    
    hdu=pyfits.open(cif)
    hdu.info()
        
    telescop=hdu[1].data['TELESCOP']
    
    #
    # check for uniqueness of telescop value - stop if not unique
    #
    tel_un=list(set(telescop))
    if not len(tel_un) == 1:
#        telescop=tel_un[0]
#    else:
        print 'TELESCOP values not unique; returning'
        print 'TELESCOP =',tel_un
        cifdict = -1
        return cifdict

        
    numrows = len(telescop)
    
    try:
        detnam=hdu[1].data['DETNAM']
    except KeyError:
        detnam = np.asarray(['INDEF']*numrows)
    try:
        filter=hdu[1].data['FILTER']
    except KeyError:
        filter = np.asarray(['INDEF']*numrows)
    try:
        instrume = hdu[1].data['INSTRUME']
    except KeyError:
        instrume = np.asarray(['INDEF']*numrows)
    files = hdu[1].data['CAL_FILE']
    qual = hdu[1].data['CAL_QUAL']
    cnam = hdu[1].data['CAL_CNAM']
    cdir = hdu[1].data['CAL_DIR']
    cdate = hdu[1].data['CAL_DATE']
    cdesc = hdu[1].data['CAL_DESC']
    cdev = hdu[1].data['CAL_DEV'] ; cdev=cdev.strip().upper()
   
    print 'Number of entries in CIF =',numrows
    
    #
    # get ONLINE files
    #
    ind = np.where(np.asarray(cdev) == 'ONLINE')[0]
    print 'Number of entries with CAL_DEV = ONLINE is {0}'.format(len(ind))
    files = files[ind]
    qual = qual[ind]
    cnam = cnam[ind]
    cdir = cdir[ind]
    cdate = cdate[ind]
    cdesc = cdesc[ind]
    instrume = instrume[ind]
    detnam = detnam[ind]
    filter = filter[ind]
    #
    # screen for quality
    #
    ind=np.where(qual <= max_cal_qual)[0]
    print 'Number of entries with CAL_QUAL <= {0} is {1}'.format(max_cal_qual,len(ind))
    files = files[ind]
    qual = qual[ind]
    cnam = cnam[ind]
    cdir = cdir[ind]
    cdate = cdate[ind]
    cdesc = cdesc[ind]
    instrume = instrume[ind]
    detnam = detnam[ind]
    filter = filter[ind]
    #
    #find uniq filenames
    #
    f=[]
    for i in range(len(cdir)):
        f.append(cdir[i].strip()+'/'+files[i].strip())
    funiq = list(set(f))
    print 'Number of unique files in caldb.indx = {0}'.format(len(funiq))
    cifdict =  {'TELESCOP':telescop,'INSTRUME':instrume, 'DETNAM':detnam, 'FILTER':filter,
    'CAL_DIR':cdir,'CAL_FILES':files,'CNAM':cnam,'CAL_QUAL':qual,'CDATE':cdate,'CDESC':cdesc}
    
    return cifdict

if __name__ == "__main__":
    acif = "http://heasarc.gsfc.nasa.gov/FTP/caldb/data/chandra/acis/index/caldbN0439.indx"
    adict = cif_summary(acif)
    print "\nkeys:"
    print adict.keys()


