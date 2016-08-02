from astropy.io import fits as pyfits
import os
import datetime
import pdb
import numpy as np
from astropy.time import Time
import glob

def get_cif(telescope,instrument, cif=''):
    """
    reads the caldb.indx file for a given telescope and instrument, for either a local caldb installation
    or a remote installation access via http or ftp
    : param telescope = telescope specification for cif
    : param instrument = instrument specification for telescope
    : param cif =  allows user to request a specific cif
            (cif = "http://heasarc.gsfc.nasa.gov/FTP/caldb/data/nustar/fpm/index/caldb.indx20160502")
    : return
    """
    try:
        caldb=os.environ['CALDB']
    except KeyError:
        print "Your $CALDB environment variable does not seem to be set"
        return
    if not cif: # if cif not specified get latest cif from telescope and instrument
        cif=caldb+"/data/"+telescope.strip().lower()+"/"+instrument.strip().lower()+"/caldb.indx"
    try:
        hdulist = pyfits.open(cif)
    except:
        hdulist=0
    return hdulist

def quizcif(telescope, instrument, cal_cnam, detnam='',cal_cbd=['','','','','','','','',''],
filter="none", cal_vsd="today", cal_vst="now", cal_qual=0, caldb='ftp://heasarc.gsfc.nasa.gov/caldb', caldbver=''):
    """
    quizcif is a python replacement for quzcif
    Returns files from the caldb.indx file for a given telescope, instrument, which matches the cal_cnam, filter
    calibration boundary etc
    """
    hdulist=get_cif(telescope,instrument, caldb=caldb)
    if hdulist==0:
        print "CIF not found for telescope %s instrument %s" % (telescope, instrument)
        print "Returning"
        return 0
    cols = hdulist[1].columns
    cifdata = hdulist[1].data
    #
    # get entries that match cal_cnam
    #
    icnam=np.where(cifdata.field("cal_cnam")==cal_cnam)
    #print cifdata[icnam].field("cal_cnam")
    cifdata=cifdata[icnam]
    print cifdata.field("cal_cnam")
    #
    # get entries that match cal_qual
    #
    iqual=np.where(cifdata.field("cal_qual")==cal_qual)
    cifdate=cifdata[iqual]
    print cifdata.field("cal_qual")
    #
    # vsd & vst
    #
    # to come...plan is to convert vsd+vst into jd, along with the date/time specified by the user and compare
    # 
    now = datetime.datetime.now()
    if cal_vsd=="today":
        vsd=now.date()
    else:
        vsd=cal_vsd
    if cal_vst=="now":
        vst=now.time()
    else:
        vst=cal_vst
    vsdt=str(vsd)+' '+str(vst)
    vsdt=Time(vsdt,format='iso',scale='utc')
    cif_vsdt=cifdata.field("cal_vsd")+' '+cifdata.field("cal_vst")
    cif_vsdt= Time(cif_vsdt, format='iso', scale='utc')

    #print vsdt.jd
    #print cif_vsdt.jd

    dift= cif_vsdt.jd-vsdt.jd
    #dift=vsdt.jd-cif_vsdt.jd
    print "DIFT =" 
    print dift
    ivsd=np.where(dift<0)
    print ivsd[0]
    if ivsd[0].size==0: # if no elements match then ivsd[0] is an array of size 0
        print "Time criterion not matched by CALDB"
        return
    cifdata=cifdata[ivsd] # filter on acceptable vsds.  Any cif_vsdt<vsdt is acceptable
    #
    # cbds - have to extract all the boundary value entries, then find the ones that match
    # the specified boundary array entries
    # for example we might specify as input: 
    # cal_cbd=['theta=50', 'phi=30'] (only equalities of input parameters are supported)
    #
    # and have cbds in the cif like cbd='theta(40-60)degrees           phi(10-50)degrees'
    #   which would match the input, while
    #   cbd='theta(40-60)degrees           phi(10-25)degrees'
    #   would not match the input (note the implicit "AND"-ing of the boundaries, not 'OR'-ing
    #
    bnd= cifdata.field("cal_cbd") # get a string array of all the boundary values
    # find all rows of bnd that contain the first input boundary parameter
    #  a. find the name of the first input boundary parameter by splitting on = equality, or < or > as necessary
    param=cal_cbd[0].split('=')
    param=param[0]
    print param
    
def get_cbds(telescope, instrument, caldb='http://heasarc.gsfc.nasa.gov/FTP/caldb'):
    """
    access the cal_cbd string (9 values of 72 characters each), extract the parameter name, then return a python dictionary
    with each parameter name as the key and the minimum and maximum values allowable as the key values
    (if an equality, then minval=maxval as expected)
    cbd={'parname1':[min,max], 'parname2':[min,max],...}
    to test for the existence of a parname:
    cbd.has_key('parname2') (evaluates to true if parname2 is a key in the dictionary)
    """
    #print cal_cbd
    #urllib.urlcleanup() # this apparently removes the file 
    hdulist=get_cif(telescope, instrument, caldb=caldb)
    tbdata = hdulist[1].data
    cbds=tbdata.field('CAL_CBD')
    nel=len(cbds)
    # return array of boundaries
    cbds=cbds.split() # split the cbd values on whitespace
    print "Splitting cbds"
    return (cbds, tbdata)

def cmp_cbd(cbd, pname, pval, punit="", chatter=0):
    """ 
    for specfied param name (pname), parameter value (pval), an optional unit string, and a caldb boundary
    list from a cif of form PARAM(VAL-RANGE)UNIT, returns true if pname = cbd parameter name AND pval is matches the 
    cbd parameter value or range, AND (optionally) the unit strings match
    
    Currently, it simply checks for a match in the parameter name for one parameter - need to add pval check
    and need to add capability to include more than one pname (i.e. for AND or OR conditions)
    """
    a=parse_cbd(cbd)
    if pname.upper()==a[0].upper(): # param names match
        if a[1]<=pval<=a[2]:
            if chatter>0:
                print pname.upper(), pval, punit, " Found"
            if punit<>a[3]: 
                print "WARNING: UNIT MISMATCH"
                print "Specified Unit: "+punit
                print "Boundary Value Unit: "+a[3]
            check=bool(1)
        else: 
            #print "No Match"
            check=bool(0)
    else:
        check=bool(0) # no param names match
    return check

def parse_cbd(bound, chatter=0):
    """
    for a boundary string of the form PARAM(VALUE)UNIT, returns the parameter name, min and max allowed value,
    (and unit string if any) as a python list (param, minval, maxval, unit)
    """
    #print "parse_cbd: Bound = %s" % bound
    test=np.asarray(bound.split("(")) # split the string at the "(" into a list; use numpy.asarray to convert to an array
    param=test[0]
    cbd=("NONE", "NONE", "NONE", "")
    if param != "NONE": 
        bval=test[1].split(")") 
        #print "\n"
        if bval[0].find("-")>=0: # .find returns a -1 if the substring is not found, otherwise the position of the string
            ptype="Range" 
            vals=bval[0].split("-")
            minval=vals[0]
            maxval=vals[1]
        else:
            ptype="Equality"
            minval=bval[0]
            maxval=bval[0]
        if chatter>0:
            print "Parameter = "+param
            print "Parameter Boundary Value = "+bval[0]
            print "Type = "+ptype
        if ptype =='Range':
            if chatter>0:
                print "Found Parameter Range"
                print("Minimum = %f, maximum = %f") % (float(minval), float(maxval))
        unit=""
        if len(bval[1]) > 0:
            unit=bval[1]
        if chatter>0:
            print "Unit = "+unit
        if minval.isdigit(): 
            minval=float(minval)
        if maxval.isdigit(): 
            maxval=float(maxval)
        cbd=(param, minval, maxval, unit)
    return cbd #returns python list (string, number or string, number or string, string)

def ck_file_existence(cif,caldb="/FTP/caldb", quiet=True):
    """
    for a cif hdulist structure (as returned, for example, by get_cif(), or from cif=pyfits.open(PATH_TO_CIF))
    retrieves the list of files from the cif and their corresponding directories
    then checks the specified caldb for them
    """
    cifdata=cif[1].data
    calfile=cifdata.field('cal_file')
    caldir=cifdata.field('cal_dir')
    files=[]
    for icf in range(len(caldir)):
        files.append(caldir[icf].strip() + '/' + calfile[icf].strip())
    ufiles = list(set(files)) # get all unique file names
    missing=[]
    for f in ufiles:
        #print "File = %s" % f
        file_path=caldb+'/'+f
        if os.path.exists(file_path):
            if not quiet:
                print "%s exists in %s" % (f,caldb)
            hdu=pyfits.open(file_path, checksum=True)
            hdu.verify('exception')
        else:
            print "\nFILE %s in caldb.indx file DOES NOT EXIST in %s" % (f,caldb)
            status = -99
            missing.append(f)
    print "Checking Complete"
    return missing

def ck_calfile_presence(filename,telescope,instrument, caldb = 'http://heasarc.gsfc.nasa.gov/FTP/caldb',
                        cif=''):
    """
    for a specified filename, searches the specified cif for the presence of the file, then
    determines if the file exists in the caldb
    :param filename:
    :param cif:
    :return:
    """
    if not cif:
        cif = caldb + '/data/' + telescope + '/' + instrument
        cifhdu = get_cif(cif=cif)
    else:
        ciuhdy = get_cif(cif=cif)
    cfiles = cifhdu[1].data['CAL_FILE']
    cal_xno = cifhdu[1].data['CAL_XNO']
    cal_cnam = cifhdu[1].data['CAL_CNAM']
    status = -1
    for cf in range(len(cfiles)):
        if filename in cfiles[cf]:
            print "{0} found in CIF {1}".format(filename, cif)
            status = 0
            #
            # verify that file exists in the caldb, at the correct location,
            # with the correct numbers of FITS extensions
            #
            # TODO: GET CALDIR, XNO, CNAM, from CIF and check that the file exists, with the correct extension and CNAM
            pass
        pass
    return




def get_calqual(cif,file):
    """
    for a cif hdulist structure (as returned, for example, by get_cif(), or from cif=pyfits.open(PATH_TO_CIF))
    retrieves the CAL_QUAL for the extensions in the file
    """
    data=cif[1].data
    calfiles = data['CAL_FILE']
    caldir = data['CAL_DIR']
    calqual = data['CAL_QUAL']
    calxno = data['CAL_XNO']
    cf=[]
    for c in calfiles:
        cf.append(c.strip()) # remove any whitespaces in the cal file name
    ind = np.where(np.asarray(cf)==file.strip())[0]
    if len(ind) == 0:
        print "File {0} not found in CIF; returning".format(file)
        xno=0
        cqual=5
        return xno,cqual
    else:
        print "File = {0}".format(file)
        xno = calxno[ind]
        cqual = calqual[ind]
        cdir = caldir[ind]
        for i in range(len(ind)):
            print "{0}/{1}[{2}] CAL_QUAL = {3} (Row = {4})".format(cdir[i].strip(),file,calxno[i],cqual[i],ind[i])
    return xno, cqual


def get_calpars (toolname, package='heasoft'):
    """
    for the given toolname (like xrtmkarf), this returns a dictionary of the parameters
    which except the string CALDB
    :param toolname:
    :return:
    """
    envname = 'LHEASOFT'
    calpars = {} # define calpars dictionary
    if package.lower().strip() == 'heasoft':
        try:
            pfiledir = os.getenv(envname)+'/syspfiles'
        except:
            print "Package = {0} but problem accessing environment variable ${1}".format(package, envname)
            print "Package = {0} but problem accessing envirionment variable ${1}".format(package, envname)
            return 0
        parfile = pfiledir+'/'+toolname+'.par'
        with open (parfile,'r') as pf:
            pars=pf.readlines()
        parnames={}
        calpars={}
        for p in pars:
            if 'CALDB' in p:
                pname = p.split(',')[0]
                pdesc = p.split(',')[6]
                parnames.update({pname:pdesc})
        calpars[toolname]=parnames
    return calpars


def test_pycaldb(dummy, caldb='ftp://heasarc.gsfc.nasa.gov/caldb'):
    import numpy as np
    print "\n TESTING cmp_cbd\n"
    cbd="PHI(10-20)arcmin"
    a= cmp_cbd(cbd,'phi',20)
    print a
    print "\n TESTING get_cif\n"
    cifhdu=get_cif('fermi','lat',cifout='/tmp/cif4',caldb=caldb)
    cols=cifhdu[1].columns
    tbdata=cifhdu[1].data
    cbd=tbdata.field('CAL_CBD')
    n=1
    cbd[len(cbd)-1][0+n*70:70+n*70] # can divide the boundary string into 9 string of length 70 characters
    bnd=['theta=100','phi=200']
    quizcif('fermi','lat','EDISP',cal_cbd=bnd, caldb=caldb)
    print "\n TESTING get_cbds and cmp_cbd\n"
    (cbds, tbdata)=get_cbds('fermi','lat', caldb=caldb)
    #print "Printing boundary values"
    cal_cnam=tbdata.field('CAL_CNAM')
    cal_qual=tbdata.field('CAL_QUAL')
    cal_file=tbdata.field('CAL_FILE')
    #Now match CNAM and CQUAL
    cnam="EDISP"
    cnamematch=np.asarray(np.where(cal_cnam.upper()==cnam.upper()))[0]
    select_index=cnamematch
    cal_qual[select_index[0]]=5 # set first quality to 5 just as a test
    cqualmatch=np.asarray(np.where(cal_qual[select_index]==0))[0]
    select_index=select_index[cqualmatch]
    #Now match boundary parameters
    for row in range(len(cbds[select_index])):
        for item in range(len(cbds[select_index[row]])):
            cbdtest=parse_cbd(item)
            if cbdtest[0]=='CLASS':
                print "Current CBD is %s, CNAM is %s (%i, %s)" % (cbds[select_index[row]][item], cal_cnam[select_index[row]], row, item)
            if cmp_cbd(cbds[select_index[row]][item], "CLASS", "P6_v1_diff_front"):
                print "\nFound parameter %s, row = %i, boundary value in row = %i" % (cbds[row][item], row, item)
                print "CBD for row %i is %s" % (select_index[row], cbds[select_index[row]])
                print "CAL_CNAM is %s" % cal_cnam[select_index[row]]
                print "CAL_FILE is %s" % cal_file[select_index[row]]             
                # if cmp_cbd
    a="Done"
    return a


if __name__ == "__main__":
    test_pycaldb(0)
