import os
import datetime
import numpy as np
from astropy.time import Time
from astropy.io import fits as pyfits
from astropy.table import Table
import sys
import urllib2
import glob
import jinja2
from heasarc.utils import utils as hu

def read_update_notice(telescope,updatedate,notice_dir = '/Users/corcoran/Desktop/caldb_updates'):
    """reads the update notice
    
    reads the update notification e-mail text file and returns a dictionary of telescope, instrument, 
    and list of files for the update.  See
    
     "Automated Delivery of Calibration Data to the CALDB" 
     (https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_2003_001/cal_gen_2003_001.html)
     for specification on the format of the update notification e-mail
    
    :param telesope: caldb telescope name
    :param updatedate: date of update in YYYYMMDD format as specified in the subject of the e-mail; 
    not to be confused with the update VERSION, also in YYYYMMDD format.  The update VERSION gives the version number
    of the individual caldb.indx file describing the index, for example the string 20150304 in caldb.indx20150304
    :param notice_dir: directory where the notification e-mail is located
    :return: update_dict, a dictionary of instruments and files for the specified update
    """
    update_email = '{0}/caldb_update_{1}_{2}.txt'.format(notice_dir,telescope, updatedate)
    update_dict = dict()
    update_dict[telescope] = dict()
    update_dict[telescope][updatedate]=dict()
    with open(update_email, mode='r') as f:
        updatelist = f.readlines()
    # get indices of instrument lines
    instindex=[]
    instrument=[]
    # identify instrument block
    for i,l in enumerate(updatelist):
        if 'instrument=' in l.lower():
            instindex.append(i)
            inst = l.strip().lower().split('=')[1]
            instrument.append(inst)
    # number of instrument blocks found
    numinst = len(instindex)
    if numinst < 1:
        print "Missing instrument block; returning"
        update_dict = -1
        return update_dict
    # at least 1 instrument block exists
    instindex.append(len(updatelist))
    for i in range(numinst):
        instblock = updatelist[instindex[i]:instindex[i+1]]
        cif = [x.strip().lower().split('caldbindexfile=')[1] for x in instblock if 'caldbindexfile=' in x]
        cif = cif[0].strip()
        if len(cif) == 0:
            print 'Calibration index file not found for instrument {0}; returning'.format(instrument[i])
            update_dict = -1
            return update_dict
        # find start and end of newfile block
        fileindx = [x for x, y in enumerate(instblock) if 'newfile' in y]
        filelist = [x.strip() for x in instblock[fileindx[0]+1:fileindx[1]]]
        # prepend name of cif to file list
        filelist.insert(0,cif)
        try:
            update_dict[telescope][updatedate][instrument[i]]=filelist
        except Exception, errmsg:
            print errmsg
            print instrument, filelist
    return update_dict

def mk_cifname(telescope, instrument, version="", caldb="http://heasarc.gsfc.nasa.gov/FTP/caldb"):
    """Creates fully specified calibration index file name
    
    Creates fully specified calibration index file name given telescope, instrument, and
    (optionally) caldb index file version (caldb.indx20170405) and caldb top level directory path or url
    
    :param telescope:  CALDB name of telescope 
    :param instrument: CALDB name of the instrument on telescope
    :param version: version of the index file to use (optional; for example caldb.indx20170405)
    :param caldb: caldb top level directory path or url
    :return: name of cif (with full directory path or full URL)
    """
    caldb      = caldb.strip()
    telescope  = telescope.strip().lower()
    instrument = instrument.strip().lower()
    version    = version.strip()

    if not version:
        cifname = "{0}/data/{1}/{2}/caldb.indx".format(caldb,telescope, instrument)
    else:
        cifname = "{0}/data/{1}/{2}/index/{3}".format(caldb,telescope, instrument, version)
    return cifname

def get_cif(telescope,instrument, version='', caldb='http://heasarc.gsfc.nasa.gov/FTP/caldb'):
    """ convert caldb.indx file to FITS HDU
    
    reads the caldb.indx file for a given telescope and instrument, for either a local caldb installation
    or a remote installation access via http or ftp
    
    :param telescope: telescope specification for cif
    :param instrument: instrument specification for telescope
    :param version: version of cif to use (for example caldb.indx20160502)
    :return: returns a FITS hdu representation of the caldb.indx file
    """
    cif = mk_cifname(telescope, instrument, version=version, caldb=caldb)
    try:
        hdulist = pyfits.open(cif)
    except Exception, errmsg:
        sys.exit('Could not open file {0} ({1})'.format(cif,errmsg))
    return hdulist

def cif_to_df(telescope, instrument, version='', caldb='http://heasarc.gsfc.nasa.gov/FTP/caldb'):
    """ convert caldb.indx file to a dataframe
    
    Creates a pandas dataframe from a calibration index file (cif) for a given mission/instrument
    
    :param telescope: caldb standard telescope name
    :param instrument: caldb standard instrument name for that telescope
    :param cif: (optional) specific cif to retrieve
    :return: pandas dataframe containing the data from the 1st extension in the CIF
    """
    cif = get_cif(telescope, instrument, version=version, caldb=caldb)
    try:
        cifdf = Table(cif[1].data).to_pandas()
    except TypeError, errmsg:
        sys.exit("Error accessing CIF {0}; exiting".format(cif))
    return cifdf

def cif_diff(telescope, instrument, newcif, oldcif,caldb='http://heasarc.gsfc.nasa.gov/FTP/caldb'):
    """get differences between new and old versions of a caldb.indx file
    
    cif_diff compares a new caldb.indx file to an older one and returns 
        - names of files in newcif which are not in the old cif
        - 
    
    :param telescope: caldb telescope name
    :param instrument: caldb instrument name
    :param newcif: version of caldb.indx file (for eg, caldb.indx20170405)
    :param oldcif: version of caldb.indx file to compare against (for eg, caldb.indx20170203)
    :param caldb: top-level caldb directory (or url)
    :return: 
    """
    # TODO - complete (do we also want to compare the files in the new/old tar files?)
    try:
        cifdf = cif_to_df(telescope,instrument,version=newcif,caldb=caldb)
    except:
        print("Error opening {0}; returning".format(newcif))
        return
    try:
        cifdf_old = cif_to_df(telescope,instrument,version=oldcif,caldb=caldb)
    except:
        print "Error opening {0}; returning".format(newcif)
        return
    files    = set(cifdf.CAL_DIR.str.strip() + '/' + cifdf.CAL_FILE.str.strip())
    oldfiles = set(cifdf_old.CAL_DIR.str.strip() + '/' + cifdf_old.CAL_FILE.str.strip())
    missing = [x for x in files if x not in oldfiles]
    deleted = [x for x in oldfiles if x not in files]
    diff_dict=dict()
    diff_dict['MISSING']=missing
    diff_dict['DELETED']=deleted
    diff_dict['CIF'] = mk_cifname(telescope, instrument, version=newcif, caldb=caldb)
    diff_dict['OLD_CIF'] = mk_cifname(telescope, instrument, version=oldcif, caldb=caldb)
    return diff_dict

def cifstats(telescope, instrument, version='', caldb='http://heasarc.gsfc.nasa.gov/FTP/caldb'):
    """ get summary statistics for a caldb.indx file
    
    returns some interesting statistics for the specified CIF:
    number of Unique filenames
    number of Unique files (a file with the same name in 2 separate directories is counted as 2 different files; 
    hopefully this doean't happen)
    number of files with quality = 0
    number of files with quality = 5
    number of files with other quality values
    number of uniques CNAM values
    
    :param telescope: caldb TELESCOPE value
    :param instrument: caldb Instrument value
    :param version: caldb.indx version (eg, caldb.indx20170405)
    :param caldb: caldb top-level directory (or url)
    :return: dictionary of summary statistics for the caldb.indx file
    """
    cifdf = cif_to_df(telescope, instrument, version=version, caldb=caldb)
    cifstat= dict()
    cifstat['Unique_filenames'] = len(set(cifdf.CAL_FILE))
    calfiles = []
    for cdir, cfile in zip(cifdf.CAL_DIR,cifdf.CAL_FILE):
        calfiles.append("{0}/{1}".format(cdir, cfile))
    cifstat['Unique_files']=len(set(calfiles))
    cifstat['Number_files_qual_0'] = len(set(cifdf[cifdf.CAL_QUAL == 0].CAL_FILE))
    cifstat['Number_files_qual_5'] = len(set(cifdf[cifdf.CAL_QUAL == 5].CAL_FILE))
    cifstat['Number_of_files_other_quality'] = cifstat['Unique_filenames'] - (cifstat['Number_files_qual_0']+cifstat['Number_files_qual_5'])
    cifstat['Number_unique_CNAMS']=len(set(cifdf.CAL_CNAM))
    return cifstat

def make_missioncif_page(telescope,instrument,version, missionurl,
                         outdir="/software/github/heasarc/pycaldb/html/cifs",
                         templatedir = "/software/github/heasarc/pycaldb/templates",cal_qual = 5,clobber=False):
    """ Make an html-formatted caldb.indx file
    
    Creates an html table from a calibration index file (CIF) the specified mission, instrument & version
    for all the "good" (CAL_QUAL = 0) entries
    
    :param telescope: caldb standard telescope name
    :param mission:  caldb standard mission name
    :param version: version of the cif (should be in YYYYMMDD format except for chandra in form of 4.7.3)
    :param outdir: output directory to hold the html version of the cif
    :param templatedir: directory where the jinja templates are located
    :param clobber: True to overwrite an existing html cif file
    :return: status (0 if no errors)
    """
    status = 0
    try:
        caldb = os.environ['CALDB']
    except KeyError:
        print "CALDB environment variable not defined; returning"
        return
    # if cifname not specified use the current cif; otherwise use the specified cif
    #
    # chandra has an overall caldb version - 4.7.3 for eg
    # but each instrument has an individual caldb name caldbN0440.indx for acis for 4.7.3
    # so we need to account for that
    if telescope.strip().lower() == 'chandra':
        ciffile = raw_input("Enter name of CIF for {0} {1} for version {2} > ".format(telescope, instrument, version))
        cifname = "{0}/data/{1}/{2}/index/{3}".format(caldb,telescope.strip().lower(), instrument.strip().lower(), ciffile.strip())
    else:
        cifname = "{0}/data/{1}/{2}/index/caldb.indx{3}".format(caldb,telescope.strip().lower(), instrument.strip().lower(),version)
    cifdf = cif_to_df(telescope.strip().lower(),instrument.strip().lower(), cif=cifname)
    if not 'FILTER' in cifdf.columns:
        print ''
        ans = raw_input('CIF missing FILTER values; set FILTER to "INDEF" ([y]/n) > ')
        if ans.strip().lower() <> 'n':
            cifdf['FILTER'] = "INDEF"
    goodcifdf = cifdf[cifdf['CAL_QUAL'] <> cal_qual]
    goodcifdf.reset_index(inplace=True) # resets the index to be 0 -> len(goodcifdf)
    # Render html file
    try:
        templateLoader = jinja2.FileSystemLoader(searchpath=templatedir)
    except Exception, errmsg:
        print ('Could not load templates in {0}; Exiting ({1})'.format(templatedir, errmsg))
        status = -1
        return status
    templateEnv = jinja2.Environment(loader=templateLoader, trim_blocks=True,)
    # define filter for DETNAM since DETNAM not included in chandra cifs
    # see http://jinja.pocoo.org/docs/2.9/api/#writing-filters
    def set_detnam (row):
        try:
            detnam = row.DETNAM
        except AttributeError:
            detnam = 'INDEF'
        return detnam
    templateEnv.filters['set_detnam']  = set_detnam
    template = templateEnv.get_template('missioncif_table_template.html')
    output_html = template.render(cifdf = goodcifdf, telescope=telescope,
                                  instrument=instrument, version=version,
                                  cifname=cifname, missionurl=missionurl)
    outname = "{0}_{1}_{2}_{3}.html".format("cif",telescope.strip().lower(),instrument.strip().lower(),version.strip().lower())
    fname = '{0}/{1}'.format(outdir, outname)
    if os.path.isfile(fname):
        if not clobber:
            sys.exit("File {0} exists, and clobber = False; not overwritten".format(fname))
    try:
        with open(fname, mode='w') as fout:
            fout.write(output_html)
    except IOError, errmsg:
        print('Could not write {0}; Exiting ({1})'.format(fname, errmsg))
        status = -1
        return status
    print "Wrote {0}".format(fname)
    return status


    print "Wrote {0}/{1}".format(outdir, outname)
    return


def quizcif(telescope, instrument, cal_cnam, detnam='',cal_cbd=['','','','','','','','',''],
            filter="none", cal_vsd="today", cal_vst="now", cal_qual=0,
            caldb='ftp://heasarc.gsfc.nasa.gov/caldb', caldbver=''):

    """ Function to retrieve calibration files matching specified constraints (UNDER DEVELOPMENT)
    
    quizcif is a python replacement for quzcif
    Returns files from the caldb.indx file for a given telescope, instrument, which matches the cal_cnam, filter
    calibration boundary etc

    :param telescope: name of telescope/mission
    :param instrument: name of instrument on telescope
    :param cal_cnam: caldb codename of files
    :param detnam: detector name (if needed)
    :param cal_cbd: calibration file boundary values
    :param filter: Name of filter for instrument (if needed)
    :param cal_vsd: Validity start date
    :param cal_vst: Validity start time
    :param cal_qual: Quality flag for files 
    :param caldb: Main caldb directory
    :param caldbver: Version of calibration file (usually of form YYYYMMDD)
    :return: 
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
    """ Get calibration boundaries
    
    Accesses the calibration boundary values for a specified telescope/instrument
    
    :param telescope: Name of Telescope/Mission (nustar)
    :param instrument: Name of Instrument (fpm)
    :param caldb: location of caldb data
    :return: calibration boundary values and table of CIF data
    """
    #print cal_cbd
    #urllib.urlcleanup() # this apparently removes the file 
    hdulist=get_cif(telescope, instrument, caldb=caldb)
    tbdata = hdulist[1].data
    cbds=tbdata.field('CAL_CBD')
    #nel=len(cbds)
    # return array of boundaries
    cbds=cbds.split() # split the cbd values on whitespace
    print "Splitting cbds"
    return (cbds, tbdata)

def cmp_cbd(cbd, pname, pval, punit="", chatter=0):
    """ compare calibration boundary values
    
    for specfied param name (pname), parameter value (pval), an optional unit string, and a caldb boundary
    list from a cif of form PARAM(VAL-RANGE)UNIT, returns true if pname = cbd parameter name AND pval is matches the 
    cbd parameter value or range, AND (optionally) the unit strings match
    
    Currently, it simply checks for a match in the parameter name for one parameter - need to add pval check
    and need to add capability to include more than one pname (i.e. for AND or OR conditions)

    :param cbd: list of calibration boundaries (from get_cbds, for example)
    :param pname: name of parameter to test
    :param pval: value of named parameter to test
    :param punit: physical unit (optional)
    :param chatter: level of verbosity
    :return: check value (True if parameter value falls within one of the tabulated boundaries)
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
    """ parse a boundary value
    
    for a boundary string of the form PARAM(VALUE)UNIT, returns the parameter name, min and max allowed value,
    (and unit string if any) as a python list (param, minval, maxval, unit)

    :param bound: caldb boundary value string
    :param chatter: level of verbosity
    :return: list of (param, minval, maxval, unit)
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
    """Check that file in caldb.indx exists in CALDB
    
    for a cif hdulist structure (as returned, for example, by get_cif(), or from cif=pyfits.open(PATH_TO_CIF))
    retrieves the list of files from the cif and their corresponding directories
    then checks the specified caldb for them

    :param cif: HDU representation of the CALDB index file
    :param caldb: location of CALDB data
    :param quiet: if True suppress status messages
    :return: status flag and list of missing files
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
    return status, missing

def ck_calfile_incif(caldir, filename,telescope,instrument, caldb = 'http://heasarc.gsfc.nasa.gov/FTP/caldb',
                        version=''):
    """ verifies that a specified file is actually in the caldb.indx file
    
    for a specified filename, checks that this file exists in the specified calibration index file (cif),
    then checks that the CALDB keyword values are consistent between the cif and the actual file 
    determines if the file exists in the caldb
    
    usage:
    ck_calfile_incif('data/hitomi/sxi/cpf/background','ah_sxi_nxbpnt_20140101v001.pha',
    caldb = 'http://heasarc.gsfc.nasa.gov/FTP/caldb',cif='caldb.indx20170405')
    
    :param caldir: path to file in the caldb relative to caldb top level directory/url
    :param filename: specified file 
    :param telescope: Name of TELESCOP
    :param instrument: Name of Instrumentn
    :param caldb: caldb top level directory/url
    :param version: version of cif to use (for example caldb.indx20170405); if blank use caldb.indx
    :return: 
    """
    # get telescope and instrument from the caldir
    telescope = caldir.split('/')[1]
    instrument = caldir.split('/')[2]

    if not version:
        cifhdu = get_cif(telescope, instrument)
    else:
        cif = get_cif(telescope, instrument, version=version)
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
    """ Get Quality Flag for caldb file
    
    For a cif hdulist structure (as returned, for example, by get_cif(), or from cif=pyfits.open(PATH_TO_CIF))
    and for a specified file name, retrieves the CAL_QUAL for all the extensions in the file

    :param cif: HDU representation of the CALDB index 
    :param file: name of file for which quality flag is returned
    :return: list of extensions and caldb quality values for those extensions
    """
    data=cif[1].data
    calfiles = data['CAL_FILE']
    caldir = data['CAL_DIR']
    calqual = data['CAL_QUAL']
    calxno = data['CAL_XNO']
    #cf=[]
    #for c in calfiles:
    #    cf.append(c.strip()) # remove any whitespaces in the cal file name
    cf = [x.strip() for x in calfiles]
    ind = np.where(np.asarray(cf)==file.strip())[0]
    if len(ind) == 0:
        print "File {0} not found in CIF; returning".format(file)
        xno=0
        cqual= -99
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
    """ parameters that accept the CALDB string as input
    
    For the given toolname (like xrtmkarf), this returns a dictionary of the parameters
    which except the string CALDB
    
    :param toolname: name of HEASoft tool/program
    :param package: Name of software package
    :return: dictionary of tools and the parameter names that accept CALDB as input
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
            pars = pf.readlines()
        parnames = {}
        calpars = {}
        for p in pars:
            if 'CALDB' in p:
                pname = p.split(',')[0]
                pdesc = p.split(',')[6]
                parnames.update({pname:pdesc})
        calpars[toolname]=parnames
    return calpars

def create_caldb_tar(telescop, instrume, version, tarName = "",
                     tardir = "", exclude_calqual = 5,
                     caldbdirs=['bcf','cpf'],
                     caldb = "http://heasarc.gsfc.nasa.gov/FTP/caldb"):
    """ Create CALDB tar file
    
    creates a tar file of "good" files for given mission/instrument
    from the HEASARC caldb.

    This should provide a means to generate the tar file so a) caldbmgr can generate
    new tar files for new caldb updates and b) allow users the ability to
    generate a tar file of previous caldb versions to install locally

    so to do this we want to
        a) access the caldb.indx file for a specified version
        b) read the file, get the list of all the good files for that version
        c) create the tar file for the goodfiles (handle case of remote or local caldb)

    :param mission: name of the mission in the caldb (str) ex: nustar
    :param instrume: name of the mission's instrument (str) ex: fpma
    :param version: version of caldb.indx file to use to create tar file (for example caldb.indx20170504)
    :param tarName: name of tarfile to be created - will be constructed if not specified
    :param tardir: name of output directory for tarfile; if not specified use CWD
    :param exclude_calqual: quality value to exclude 
    (some missions like Chandra use cal_qual values from 0 to 4 as acceptable values)
    :param caldbdirs: subdirectories of $CALDB/data/<telescop>/<instrume> to tar 
        (the index directory is always retrieved)
    :param caldb: top-level caldb directory or url
    :return: status (0 if no errors)
    """
    import tarfile
    from ftplib import FTP
    import shutil
    # TODO: allow selection of one or more subdirectories to tar; useful for clockcor
    #cif = get_cif(telescop, instrume)
    tel=telescop.strip().lower()
    instr = instrume.strip().lower()
    ver = version.strip()
    # ver should be of the form 20170503
    # ver = version.split('caldb.indx')[1].strip()
    cwd = os.getcwd()
    status = 0
    exclude_calqual = int(exclude_calqual) # make sure this is an integer
    if not tardir:
        tardir = os.getcwd()
    else:
        try:
            os.chdir(tardir)
        except OSError:
            print "Could not change directory to {0}".format(tardir)
            print "returning"
            status = -1
            return
        print 'Current directory = {0}'.format(os.getcwd())

    cif = "{0}/data/{1}/{2}/index/caldb.indx{3}".format(caldb,tel,instr,ver)
    try:
        hdu = pyfits.open(cif)
    except Exception, errmsg:
        print("Could not open {0}; returning".format(cif))
        status = -1
        return status
    if not tarName: # if tarname not specified, use "standard" naming
        tarName = "{3}/goodfiles_{0}_{1}_{2}.tar.gz".format(tel, instr, ver, tardir)
    tar = tarfile.open(tarName, "w:gz")
    if (("http://" in caldb) or ("ftp://" in caldb)): # if remote, download the index file directory
         tindname = "{0}/index.tar".format(tardir)
         tindex = open(tindname, 'w')
         ftp = FTP("heasarc.gsfc.nasa.gov")
         ftp.login(user="anonymous", passwd ="")
         # heasarc ftp allows retrieval of an entire directory if .tar appended to directory name
         ftp.retrbinary("RETR /FTP/caldb/data/{0}/{1}/index.tar".format(tel, instr), tindex.write)
         tindex.close()
         t = tarfile.open(tindname,'r')
         try:
             t.extractall()
         except Exception, errmsg:
             print "Problem extracting files from index.tar:{0}".format(errmsg)
             status = -1
         cifs = glob.glob('index/caldb*')
         arcname = ['index/' + c[c.rfind('/') + 1:] for c in cifs]
         for c, a in zip(cifs, arcname):
             tar.add(c, "data/{0}/{1}/{2}".format(tel, instr, a))
         if os.path.islink("caldb.indx"):
             print "Removing {0}/caldb.indx".format(tardir)
             os.remove("{0}/caldb.indx".format(tardir))
         try:
             os.symlink("index/caldb.indx{0}".format(ver), "caldb.indx")
         except OSError, emsg:
             print "Could not create symbolic line to caldb.indx"
             print "Error message was: {0}".format(emsg)
             print "returning"
             status = -1
             return status
         tar.add("caldb.indx", "data/{0}/{1}/caldb.indx".format(tel, instr))
         shutil.rmtree(tardir+'/index') # remove local index directory
         os.remove(tardir+'/index.tar')
    else:
        cifdir = "{0}/data/{1}/{2}/index".format(caldb, tel, instr)
        cifs = glob.glob("{0}/caldb*".format(cifdir))
        arcname = ['index/'+c[c.rfind('/')+1:]  for c in cifs]
        for c, a in zip(cifs, arcname):
            tar.add(c, "data/{0}/{1}/{2}".format(tel, instr, a))
        if os.path.isfile("{0}/caldb.indx".format(tardir)):
            print "Removing {0}/caldb.indx".format(tardir)
            os.remove("{0}/caldb.indx".format(tardir))
        try:
            os.symlink("index/caldb.indx{0}".format(ver), "caldb.indx")
        except OSError, errmsg:
            lcif = "{0}/caldb.indx".format(tardir)
            print errmsg
            print "Does file {0} exist? {1}".format(lcif, os.path.isfile(lcif))
            status = -1
        tar.add("caldb.indx","data/{0}/{1}/caldb.indx".format(tel, instr))
    os.remove("{0}/caldb.indx".format(tardir))
    # get calfiles
    ciftab = hdu[1].data
    cqual = ciftab['CAL_QUAL']
    cfiles = ciftab['CAL_FILE']
    cdir = ciftab['CAL_DIR']
    # get cfile_to_tar, list of unique caldb files to tar based on specified cal quality
    # and create goodfiles array
    cfile_to_tar = []
    cdir_to_tar = []
    goodfiles = []
    for cq, cf, cd in zip(cqual, cfiles, cdir):
        testfile = "{0}/{1}".format(cd,cf)
        if (cq <> exclude_calqual) and (testfile not in cfile_to_tar):
            cfile_to_tar.append(testfile)
            cdir_to_tar.append(cd)
            goodfiles.append("{0}".format(testfile))
    # now create list of selected goodfiles based on the specified caldbdirs list
    goodfiles_selected = []
    for c in caldbdirs:
        tmp = [x for x in goodfiles if c in x]
        goodfiles_selected.extend(tmp)
    # goodfiles_selected now contains the list of qualtity selected caldb files in the specified subdirectories caldbdirs
    for goodfile in goodfiles_selected:
        if ("http://" in caldb) or ("ftp://" in caldb): # if remote access, download the file
            # download the goodfile from the heasarc
            response = urllib2.urlopen("{0}/{1}".format(caldb,goodfile))
            data = response.read()
            localgoodfile = goodfile[goodfile.rfind("/")+1:] # strip off the directory path for the local file
            tmp_cfile = "{0}/{1}".format(tardir, localgoodfile)
            out = open(tmp_cfile,'w')
            out.write(data)
            out.close()
            print "Adding {0} as {1}".format(tmp_cfile, goodfile)
            try:
                tar.add(tmp_cfile, arcname=goodfile)
            except Exception, errmsg:
                print "Error adding {0}".format(goodfile)
                print errmsg
                status = -1
            os.remove(tmp_cfile)
        else: # using locally-mounted caldb
            #tmp_cfile = "{0}/{2}".format(caldb, cd, cf)
            cfile_to_add = "{0}/{1}".format(caldb,goodfile)
            print "Adding {0} as {1}\n".format(cfile_to_add, goodfile)
            try:
                tar.add(cfile_to_add, arcname=goodfile)
            except Exception, errmsg:
                print "Error adding {0}".format(goodfile)
                print errmsg
                status = -1

    tar.close()
    os.chdir(cwd)
    return status

def create_update_tarfile(telescope, update, notice_dir_root='/FTP/caldb/staging/data/'):
    """ creates tarfile of latest update files 
    
    Creates tarfile of the files in the latest update files.  This file can be untarred in $CALDB to update 
    the previous version of a local caldb to the new version
    
    :param telescope: CALDB mission designation (example: 'swift')
    :param update: CALDB update designation (example: '20170501')
    :param notice_dir_root: directory root where update notices are stored 
    :return: 
    """
    telescope = telescope.lower().strip()
    update    = update.lower().strip()
    notice_dir_root = notice_dir_root.lower().strip()
    upd = read_update_notice(telescope, update, notice_dir="{0}/{1}/to_ingest".format(notice_dir_root,
                                                                                      telescope))
    date = upd[telescope].keys()[0]
    inst = upd[telescope][date].keys()[0]
    tarname = 'update_{0}_{1}_{2}.tar'.format(telescope,inst, date)
    inflist = ["/FTP/caldb/data/{1}/{2}/{0}".format(x, telescope, inst) for x in upd[telescope][date][inst][1:]]
    outlist = [x.replace('/FTP/caldb/data','data') for x in inflist]
    hu.create_tarfile(inflist, tarname, tarlist=outlist)
    return

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

def sync_caldbweb(DoCopy=False, user='mcorcora', host="heasarcdev", mission='', version=''):
    """ Syncs the caldb website
    
    After editing/updating appropriate html files on the development side (caldbwww), this routine moves 
    the files to the public side (caldbwwwprod) and also generates the cadlb and heasarc rss feeds
    
    The html files that need to be updated are::
    
        caldb_wotsnew.html
        caldb.rss
        caldb_supported_missions.html
        the astro-update associated data file (inc/associated_data.html)
        
    For nustar updates, you also need to edit these files::
    
        release_<VERSION>.txt
        nustar_caldbhistory.html
        nustar_caldb.html
        nustar.rss

    
    :param mission: mission name (aka TELESCOP)
    :param version: version of caldb (usually of form YYYYMMDD)
    :param DoCopy: if False just print the cmds that will be executed when DoCopy = True
    :return: 
    """
    import os
    mission = mission.lower().strip()
    errlist=[]
    curuser = os.environ['LOGNAME']
    curhost = os.environ['HOSTNAME']
    if (curuser<>user) or (host not in curhost ):
        error = "This function must be run as {0} from {1}; returning".format(user, host)
        print(error)
        print("Current user is:{0}  Current host is: {1}".format(curuser, curhost))
        errlist.append(error)
        return errlist
    # generate the caldb news item
    print "# Generating CALDB news item (develop version)"
    cmd = "~mcorcora/bin/caldbdev_rss2.pl"
    print cmd
    if DoCopy:
        status = os.system(cmd)
        if status <> 0:
            errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    print "# Generating CALDB news item (public version)"
    cmd = "~mcorcora/bin/caldbprod_rss2.pl"
    print cmd
    if DoCopy:
        status = os.system(cmd)
        if status <> 0:
            errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    # each element of files_to_sync has 3 components:
    # name of file, start directory, end directory
    files_to_sync = [
        ("caldb_wotsnew.html","/www/htdocs/docs/heasarc/caldb","/www.prod/htdocs/docs/heasarc/caldb"),
        ("caldb_supported_missions.html", "/www/htdocs/docs/heasarc/caldb", "/www.prod/htdocs/docs/heasarc/caldb"),
        ("caldb.rss", "/www/htdocs/docs/heasarc/caldb", "/www.prod/htdocs/docs/heasarc/caldb"),
        ("associated_data.html","/www/htdocs/docs/heasarc/astro-update/inc","/www.prod/htdocs/docs/heasarc/astro-update/inc"),
        ("astro-update.rss","/www/htdocs/docs/heasarc/astro-update","/www.prod/htdocs/docs/heasarc/astro-update")
    ]
    if mission == 'nustar':
        files_to_sync.append(
            [
                ("release_{0}.txt".format(version), "/www/htdocs/docs/heasarc/caldb/nustar/docs/","/www.prod/htdocs/docs/heasarc/caldb/nustar/docs/"),
                ("nustar_caldbhistory.html","/www/htdocs/docs/heasarc/caldb/nustar/docs/","/www.prod/htdocs/docs/heasarc/caldb/nustar/docs/"),
                ("nustar_caldb.html","/www/htdocs/docs/heasarc/caldb/nustar","/www.prod/htdocs/docs/heasarc/caldb/nustar"),
                ("cif_nustar_fpm_{0}.html".format(version),"/FTP/caldb/local/software/pro/DATA_DELIVERIES/pro/mission_summaries/work","/www/htdocs/docs/heasarc/caldb/data/nustar/fpm/index"),
                ("cif_nustar_fpm_{0}.html".format(version), "/FTP/caldb/local/software/pro/DATA_DELIVERIES/pro/mission_summaries/work","/www.prod/htdocs/docs/heasarc/caldb/data/nustar/fpm/index"),
                ("nustar.rss","/www/htdocs/docs/nustar/news/","/www.prod/htdocs/docs/nustar/news/")
            ]
        )
    print "\n # Syncing website:"
    for f in files_to_sync:
        print "# copying {0} from {1} to {2}".format(f[0],f[1],f[1])
        cmd = "cp {1}/{0} {2}/{0}".format(f[0],f[1],f[2])
        print "{0}\n".format(cmd)
        if DoCopy:
            status = os.system(cmd)
            if status<>0:
                errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    if mission == "nustar":
        print "# Creating {0} and heasarc newsfeeds (dev version)".format(mission)
        cmd = "/www/htdocs/docs/rss/nustar_rss2.pl"
        print("{0}\n".format(cmd))
        if DoCopy:
            status = os.system(cmd)
            if status <> 0:
                errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
        print "# Creating nustar and heasarc newsfeeds (public version)"
        cmd = "/www.prod/htdocs/docs/rss/nustar_rss2.pl"
        print("{0}\n".format(cmd))
        if DoCopy:
            status = os.system(cmd)
            if status <> 0:
                errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    # create heasarc news feed from astroupdate
    print "\n # Creating heasarc newsfeed (dev version)".format(mission)
    cmd = "/www/htdocs/docs/rss/software_rss2.pl"
    print("{0}\n".format(cmd))
    if DoCopy:
        status = os.system(cmd)
        if status <> 0:
            errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    print "\n #Creating heasarc newsfeed (public version)"
    cmd = "/www.prod/htdocs/docs/rss/software_rss2.pl"
    print("{0}\n".format(cmd))
    if DoCopy:
        status = os.system(cmd)
        if status <> 0:
            errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    return errlist

def test_makemissioncif(telescope='Swift',instrument='SC', version='20170505',missionurl='http://swift.gsfc.nasa.gov'):
    make_missioncif_page(telescope, instrument,version, missionurl=missionurl, clobber=True)
    return

if __name__ == "__main__":
    #create_caldb_tar('nustar','fpm','20161021', tardir='/Users/corcoran/Desktop/tmp/caltartest',
    #                 tarName='goodfiles_nustar_fpm_94nov9.tar.gz',
    #                 caldb= '/caldb')
    # test
    #create_caldb_tar('swift', 'xrt', '20160609 ', tardir='/FTP/caldb/staging/tmp', caldb='/FTP/caldb')
    #test_pycaldb(0)
    # dummy  command
    #test_makemissioncif()
    #cifdf = cif_to_df('nustar','fpm')
    # cstats = cifstats('nustar','fpm')
    # cstatsold = cifstats('nustar','fpm', version='caldb.indx20160606')
    # for k in cstats.keys():
    #     print "{0:30s} latest = {1} 20160606 = {2}".format(k, cstats[k],cstatsold[k])
    ud = read_update_notice('swift', '20081203')
    for k in ud['swift'].keys():
        for i in ud['swift'][k].keys():
            print k, i, ud['swift'][k][i]
