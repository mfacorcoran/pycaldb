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
from heasarc import utils as hu
# from heasarc.utils.webutils import webutils as hw
import time
import ftplib

class Caldb(object):
    """
    Defines a CALDB object; optionally specify the telescope (& instrument)
    """
    def __init__(self,caldb=None,
                 caldbconfig = None,
                 caldbalias = None,
                 telescope = None, instrument = None):
        """ Returns a CALDB object
        A CALDB object defines the structure of a calibration database.  The structure is defined as:
           * the location of the caldb (i.e. the value of caldb environment variable)
           * the path+name of the caldb.config file
           * the path+name of the caldb alias file
           * the CALDB data directories (list of the directories directly under $CALDB/data)
        :param telescope: CALDB telescope name (TELESCOP keyword)
        :param instrument: CALDB instrument name (INSTRUME keyword)
        :param version:
        :param caldb: name of caldb top level directory; if None get from the CALDB environment
        :param caldbconfig: name of caldb config file; if None get from the CALDB environment
        :param caldbalias: name of caldb alias file; if None get from the CALDB environment
        """
        if not caldb:
            try:
                caldb = os.environ['CALDB']
            except KeyError, errmsg:
                print ("{0} missing ".format(errmsg))
                caldb = ''
        if not caldbconfig:
            try:
                caldbconfig = os.environ['CALDBCONFIG']
            except KeyError, errmsg:
                print ("{0} missing ".format(errmsg))
                caldbconfig = ''
        if not caldbalias:
            try:
                caldbalias = os.environ['CALDBALIAS']
            except KeyError, errmsg:
                print ("{0} missing ".format(errmsg))
                caldbalias = ''
        self.caldb = caldb.strip()
        self.caldbconfig = caldbconfig
        self.caldbalias = caldbalias
        self.telescope = telescope
        self.instrument = instrument
        self.is_installed = None
        self.columns = ['TELESCOP', 'INSTRUME', 'DETNAM', 'FILTER', 'CAL_DEV', 'CAL_DIR',
                        'CAL_FILE', 'CAL_CLAS', 'CAL_DTYP', 'CAL_CNAM', 'CAL_CBD',
                        'CAL_XNO', 'CAL_VSD', 'CAL_VST', 'REF_TIME', 'CAL_QUAL',
                        'CAL_DATE', 'CAL_DESC']
        self.keywords =['TELESCOP', 'INSTRUME',
                         'DETNAM','FILTER','CCLS*',
                         'CDTP*','CCNM*','CBD1*',
                         'CVSD*','CVST*','CDES*']
        if self.telescope:
            # if telescope defined check that the current caldb is configured for it
            if not self.check_config():
                #print "Configuration issue with Telescope {0} or instrument {1}".format(self.telescope, self.instrument)
                print "Use set_telescope(), set_instrument() methods to set telescope and/or instrument attributes"
                self.telescope = None
                self.instrument = None
        return


    def get_config(self):
        """
        returns a dictionary of the caldb configuration as documented in caldb.config
        :return:
        """
        cfig = get_caldbconfig(self.caldbconfig)
        if self.telescope:
            cfig=cfig[self.telescope]
        return cfig


    def check_config(self):
        """
        returns True if the CALDB is configured for the specified telescope and/or instrument, False otherwise
        :return:
        """
        # read the caldb config file
        cfig = get_caldbconfig(self.caldbconfig)
        if self.telescope:
            if self.telescope in cfig.keys():
                configured = True
            else:
                configured = False
            if configured:
                print "caldb.config configured for TELESCOP = {0}".format(self.telescope)
                # if the caldb is configured for telescope, check if its also configured
                # for the telescope's instrument, if the instrument is specified
                if self.instrument:
                    if self.instrument in cfig[self.telescope].keys():
                        print(" ... and for INSTRUME = {inst}".format(inst=self.instrument))
                        configured = True
                    else:
                        configured = False
                        print("  ... but not configured for instrument {inst} of telescope {tel}".format(inst=self.instrument,
                                                                                                      tel = self.telescope))
                        print "Available INSTRUME values for {0}:".format(self.telescope)
                        print(" ".join(cfig[self.telescope].keys()))

        return configured

    def get_telescopes(self):
        """
        returns the names of the missions defined in the caldb.config file
        :return:
        """
        return self.get_config().keys()

    def set_telescope(self, mission, verbose=True):
        """
        sets the caldb to the caldb for a specific mission;
        prints the available instruments if verbose = True
        :return:
        """
        cfig= get_caldbconfig(self.caldbconfig)
        try:
            missions = cfig.keys()
        except KeyError:
            print "Could not get missions from caldbconfig file '{0}'".format(self.caldbconfig)
            return
        if mission in missions:
            self.telescope = mission
            self.instrument = ''
            print "Setting telescope to {0}; Available instruments are".format(mission)
            print "  ".join(cfig[mission].keys())
        else:
            print "Mission {0} not found in {1}".format(mission, self.caldb)
            print "Available missions are: "
            print " ".join(missions)
            return
        return


    def get_instruments(self, verbose=True):
        """
        if self.telescope defined, returns the list of defined instruments from the caldb.config file
        :return:
        """
        #cfig = get_caldbconfig(self.caldbconfig)
        if self.telescope:
            cfig = self.get_config()
            instruments = cfig.keys()
            if verbose:
                print "Available instruments for {0} are".format(self.telescope)
                print ' '.join(instruments)
        else:
            print("Telescope not set for Caldb object.  Use set_telescop() method to set it")
            instruments=''
        return instruments

    def set_instrument(self, instrument):
        """
        for one of the defined instruments
        gets the versions of the caldb.indx file available
        :return:
        """
        if not self.telescope:
            print "Telescope not set in Caldb object; use set_telescope() to set it"
            return
        try:
            instruments = get_caldbconfig(self.caldbconfig)[self.telescope]
        except:
            print "Problem getting instruments for {0}; returning".format(self.telescope)
            return
        if instrument in instruments:
            self.instrument = instrument
        else:
            print "Instrument {0} is not defined for Telescope {1}".format(instrument, self.telescope)
            print "Available instruments for {0} are".format(self.telescope)
            for i in instruments:
                print i
        return

    def get_insdir(self):
        """
        gets the directory in the caldb that corresponds to the specified instrument
        :return: path to instrument directory
        """
        if not self.telescope:
            print "Telescope not defined for Caldb object; use set_telescope() to set it"
            return
        if not self.instrument:
            print "Instrument not defined for Caldb object; use set_instrument() to set it"
            return
        instrument  = self.instrument
        insdir = os.path.join(self.caldb,self.get_config()[instrument][1])
        return insdir

    def get_versions(self):
        """Get the versions of the caldb.indx files available for the telescope and instrument

        Retrieves the available versions of the caldb.indx files as listed in the index subdirectory

        :return: returns the available caldb.indx* files from the index directory as a list
        """

        # TODO - allow user to return the version that's currently being used (may not be the latest version)

        # get versions of available caldb.indx files
        try:
            insdir = self.get_insdir()
        except AttributeError:
            print("Error getting instrument directory for {tel} & {inst}; returning".format(tel=self.telescope, inst=self.instrument))
            return []
        try:
            indxdir = os.path.join(insdir, 'index')
        except AttributeError:
            print("Error getting instrument directory for {tel} & {inst}; returning".format(tel=self.telescope, inst=self.instrument))
            return []
        if "heasarc.gsfc.nasa.gov" in self.caldb:
            # using remote access so need to get caldb.indx file list using ftplib
            ftp = ftplib.FTP("heasarc.gsfc.nasa.gov")
            ftp.login("anonymous", "anonymous")
            try:
                wdir = indxdir.split('heasarc.gsfc.nasa.gov/')[1] # get heasarc directory path
                ftp.cwd(wdir)
                versions = ftp.nlst()
                ftp.close()
                return versions
            except Exception, errmsg:
                print "problem accessing {0} on heasarc.gsfc.nasa.gov via FTP {1}".format(indxdir,errmsg)
        else:
            indxdir = os.path.join(insdir,'index')
            if os.path.exists(indxdir):
                versions = glob.glob("{0}/*".format(indxdir))
                try:
                    versions = [os.path.split(x)[1] for x in versions]
                except:
                    print "problem accessing {0}".format(indxdir)
                return versions
            else:
                #local indxdir does not exist
                print "{0} does not exist".format(indxdir)
                return

    def set_cif(self, version=None):
        """Sets the cif attribute of the Caldb object to the specified version; if none, use current caldb.indx file
        :param version: optional, version of the caldbindx file to use (eg: 'caldb.indx20170727')
        :return: sets the cif attribute to the specified caldb version
        """
        if version:
            self.cif = os.path.join(self.get_insdir(),'index',version)
        else:
            self.cif = os.path.join(self.get_insdir(),'caldb.indx')
        self.version = version
        return

    def get_cif(self):
        """
        get the calibration index file as a Pandas Datafram

        :return: Pandas dataframe of the cif
        """
        cifhdu = self.get_cifhdu()
        try:
            cifdf = Table(cifhdu[1].data).to_pandas()
        except TypeError, errmsg:
            print "Error converting CIF {0} to Table; exiting".format(self.cif)
            return
        if 'INSTRUME' not in cifdf.columns:
            cifdf['INSTRUME'] = 'INDEF'
        # strip whitespace from string columns
        stripcols = cifdf.columns
        for s in stripcols:
            try:
                cifdf[s] = [x.strip() for x in cifdf[s]]
            except:
                pass
        # make the CBDs a list of lists rather than a string
        cbds = [x.split() for x in cifdf.CAL_CBD]
        cifdf.CAL_CBD = cbds
        return cifdf

    def get_cifhdu(self, cif = None):
        if not cif:
            try:
                print "Retrieving {0}".format(self.cif)
            except AttributeError:
                print "cif needs to be set using the set_cif() method for the CALDB instance"
                return
        try:
            cifhdu = pyfits.open(self.cif)
        except Exception, errmsg:
            print ('Could not open file {0} ({1})'.format(self.cif, errmsg))
            cifhdu = None
        return cifhdu


    def mkcif(self, cifout, src_cif=None, outdir=None, New=False, clobber=False):
        """
        creates a copy of cif named cifout in the caldb for a given telescope and instrument; if
        cifdir is not specified, use the index directory for the specified telescope and instrument

        Uses:
            1) create a new cif for a new telescope/instrument: in this case the telescope & instrument
            values are defined, a caldb directory structure with files has been set up
            2) create a copy of an existing cif which may then be updated with new caldb files

        :param telescope: name of the telescope/mission (value of the TELESCOP keyword)
        :param instrument: name of the instrument (value of the INSTRUME keyword)
        :param src_cif: name of cif to be copied
        :param outdir: location to write the cif; if not specified, define from the caldbconfig file
        :param New: if True, create a blank cif
        :return:
        """

        # if cifdir not defined, create cif director from caldbconig file

        if not outdir:
            outdir= os.path.join(self.get_insdir(),'index')
        if not self.check_config():
            print ("Returning")
            return

        if New:
            c1 = pyfits.Column(name='TELESCOP', format='10A     ')
            c2 = pyfits.Column(name='INSTRUME', format='10A     ')
            c3 = pyfits.Column(name='DETNAM  ', format='20A     ')
            c4 = pyfits.Column(name='FILTER  ', format='10A     ')
            c5 = pyfits.Column(name='CAL_DEV ', format='20A     ')
            c6 = pyfits.Column(name='CAL_DIR ', format='70A     ')
            c7 = pyfits.Column(name='CAL_FILE', format='40A     ')
            c8 = pyfits.Column(name='CAL_CLAS', format='3A      ')
            c9 = pyfits.Column(name='CAL_DTYP', format='4A      ')
            c10 = pyfits.Column(name='CAL_CNAM', format='20A     ')
            c11 = pyfits.Column(name='CAL_CBD ', format='630A70  ')
            c12 = pyfits.Column(name='CAL_XNO ', format='I       ')
            c13 = pyfits.Column(name='CAL_VSD ', format='10A     ')
            c14 = pyfits.Column(name='CAL_VST ', format='8A      ')
            c15 = pyfits.Column(name='REF_TIME', format='D       ')
            c16 = pyfits.Column(name='CAL_QUAL', format='I       ')
            c17 = pyfits.Column(name='CAL_DATE', format='10A     ')
            c18 = pyfits.Column(name='CAL_DESC', format='70A     ')

            cols = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9,
                                 c10, c11, c12, c13, c14, c15, c16, c17, c18])

            hdu0 = pyfits.PrimaryHDU()
            hdu1 = pyfits.BinTableHDU(from_columns(cols))

            hdu1.header['TELESCOP'] = self.telescope
            hdu1.header['EXTNAME'] = 'CIF'
            hdu1.header['CREATOR'] = 'Pycaldb CRCIF'
            hdu1.header['CIFVERSN'] ='1.1'
            hdu1.header['TELESCOP'] = telescope
            hdu1.header['INSTRUME'] = instrument

            cif = pyfits.HDUList([hdu0, hdu1])
            cif.writeto(os.path.join(outdir,cifout))

        else:

            dst = os.path.join(outdir, cifout)
            if os.path.isfile(dst):
                print("File {dst} exists and clobber = {clob}".format(dst=dst, clob=clobber))
                if clobber:
                    print("File will be overwritten")
                else:
                    print("File not overwritten; returning")
                    return
            if not self.version:
                print "Available CIFs are:\n "
                print " ".join(self.versions())
                cifsrc =  raw_input("Enter name of CIF to copy > ")
            else:
                cifsrc = self.version
            cifsrc = os.path.join(self.get_insdir(),'index',cifsrc)
            cifhdu = self.get_cifhdu(cif = cifsrc)
            cifhdu.writeto(dst, overwrite = clobber)
            print("Copied {cifsrc} to {dst}").format(cifsrc=cifsrc, dst=dst)
            return

    def cifstats(self):
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
        try:
            cifdf = self.get_cif()
        except AttributeError:
            print "Error in retrieving named caldb index file; returning"
            return
        cifstat = dict()
        cifstat['Number_Unique_filenames'] = len(set(cifdf.CAL_FILE))
        cifstat['Unique_Instrument_names'] = [x.strip() for x in list(set(cifdf.INSTRUME))]
        calfiles = []
        for cdir, cfile in zip(cifdf.CAL_DIR, cifdf.CAL_FILE):
            calfiles.append("{0}/{1}".format(cdir, cfile))
        cifstat['Number_Unique_files'] = len(set(calfiles))
        cifstat['Number_files_qual_0'] = len(set(cifdf[cifdf.CAL_QUAL == 0].CAL_FILE))
        cifstat['Number_files_qual_5'] = len(set(cifdf[cifdf.CAL_QUAL == 5].CAL_FILE))
        cifstat['Number_of_files_other_quality'] = cifstat['Number_Unique_filenames'] - (
        cifstat['Number_files_qual_0'] + cifstat['Number_files_qual_5'])
        cifstat['Number_unique_CNAMS'] = len(set(cifdf.CAL_CNAM))
        return cifstat

    # def get_stats(self):
    #     print self.caldb
    #     return cifstats(self.telescope, self.instrument,
    #                     version = self.version, caldb=self.caldb)
    
    def get_cnames(self):
        """returns a list of all the unique, defined CAL_CNAME values from the caldb.indx file

        :param caldb: optionally specify the caldb path (URLs allowed)
        :param cifname: optionally specify the name of the caldb.indx file (uses current file if not specified)
        :return: list of all the unique defined CAL_CNAME values from the caldb.indx file
        """
        cdf = self.get_cif()
        try:
            cnames=list(set(cdf['CAL_CNAM']))
        except TypeError:
            print "Could not retrieve CIF; returning"
            return
        return cnames

    def find_calfiles(self,cname, cqual=0, return_files=True):
        """
        returns a list of files from a caldb.indx object having the specified
        CAL_CNAM and CAL_QUAL

        :param cname: CAL_CNAM value of interest (eg: "hkrange") (UPPER CASE)
        :param cqual: calibration quality value CAL_QUAL (integer; default 0)
        :param return_files: if False, returns the entire dataframe matching the specified cname, cqual;
        if True, returns only the unique filenames with fully-qualified subdirectory path
        :return:

        """
        cname = cname.upper().strip()
        cdf = self.get_cif()
        # convert to upper case for comparison
        cdfc = cdf[cdf['CAL_CNAM'] == cname]
        if len(cdfc)==0:
            print "Could not find any files with CAL_CNAM  = {0} in caldbindx".format(cname)
            return -1
        cdfcq = cdfc[cdfc['CAL_QUAL'] == cqual]
        if len(cdfcq) == 0:
            print "Could not find any files with CAL_CNAM = {0} and CAL_QUAL ={1}".format(cname, cqual)
            return -1
        if return_files:
            files = []
            for i in range(len(cdfcq)):
                fname = "{0}/{1}".format(cdfcq.CAL_DIR.iloc[i],cdfcq.CAL_FILE.iloc[i])
                files.append(fname)
            files = list(set(files)) # return unique file names
            return files
        else:
            return cdfcq

    def get_calfile_info(self, calfilename):
        """
        returns info from the caldbindx file for the specified calfilename
        :param calfilename: name of cal file of interest, with fully specified directory path relative to $CALDB
        :return:
        """
        cdf = self.get_cif()
        cfsplit = os.path.split(calfilename)
        cdir = cfsplit[0]
        cfile = cfsplit[1]
        try:
            fileinfodf = cdf[(cdf.CAL_FILE == cfile) & (cdf.CAL_DIR == cdir)]
        except AttributeError:
            print "Cif data frame has no CAL_FILE attribute; returning"
            return
        if not len(fileinfodf):
            print "Could not find {0} in {1}".format(cfile, self.cif)
            return
        # try:
        #     fileinfodf = filedf[filedf.CAL_DIR == cdir]
        # except:
        #     print "Could not find {0} in {1}".format(cdir, self.cif)
        #     return
        print "For {0}:".format(calfilename)
        print "\nInformation from {0}".format(self.cif)
        print "  Number of extensions indexed in CIF = {0}".format(max(fileinfodf['CAL_XNO']))
        fname = os.path.join(self.caldb, fileinfodf.CAL_DIR.iloc[0], fileinfodf.CAL_FILE.iloc[0])
        try:
            hdu = pyfits.open(fname)
        except Exception, errmsg:
            print "Could not retrieve {0} ({1}".format(fname, errmsg)
            return
        print "  Number of extension in File = {0}".format(len(hdu)-1)
        print "  CNAMEs for each extension from CIF: "
        for i in range(len(hdu)-1):
            fxno = fileinfodf[fileinfodf.CAL_XNO == i+1]
            cname = fxno.CAL_CNAM.iloc[0]
            cbd = fxno.CAL_CBD.iloc[0]
            print "    Extension #{0} CNAME = {1} CBDs = {2}".format(i+1, cname, ','.join(cbd))
        return

    def cif_diff(self, caldb2cmp, print_report=True):
        """
        compares the current Caldbindx cif (as given by the caldb.cif attribute) to the specified version and returns a dictionary showing:
        - MISSING_FROM_COMPARISON:names of files & extension numbers in current
              caldb.cif which are not in the specified version
        - MISSING_FROM_CURRENT: names of files & extension numbers in specified version of the cif which
            are not in the current caldb.cif
        TODO
        - CHANGED_QUALITY: names of files & extension numbers in the current Caldbindx object which have changed their quality flags compared
        to the comparison Caldbindx object TODO

        :param caldb2cmp: Caldbindx object to be compared to the current Caldbindx instance
        :param print_report: if True, prints a nicely formatted report of the differences to the screen
        :return:

        """
        #
        # TODO - add CHANGED_QUALITY
        #  how to do this:
        #     files2cmp = list(set(cif2cmp['CAL_FILE'] # get all unique file names
        #     for each file in the CIF_to_CMP, find all the files in the CIF
        #     for f in files2cmp:
        #            testcif2cmp = cif2cmp[['CAL_FILE','CAL_XNO', 'CAL_QUAL'][cif2cmp['CAL_FILE']==f]
        #            testcif = cif2cmp[['CAL_FILE','CAL_XNO', 'CAL_QUAL'][cif2cmp['CAL_FILE']==f]
        #     then see if testcif2cmp is the same as testcif.
        #            if not testcif2cmp.equals(testcif)
        #
        # If it is, move on, if not find the differences
        diff_dict = dict()
        cifdf = self.get_cif()
        diff_dict['CIF'] = self.cif


        if not self.cif:
            print "CIF not set for main caldb object; use set_cif() method to set it; returning"
            return
        if not caldb2cmp.cif:
            print "CIF not set for comparison caldb object; use set_cif() method to set it; returning"
            return

        cifdf2cmp = caldb2cmp.get_cif()
        diff_dict['ComparisonCIF'] = caldb2cmp.cif

        # get unique filenames from cif
        files = set(cifdf.CAL_DIR.str.strip() + '/' + cifdf.CAL_FILE.str.strip())

        # get unique filenames from comparison cif
        cifdf2cmpfiles = set(cifdf2cmp.CAL_DIR.str.strip() + '/' + cifdf2cmp.CAL_FILE.str.strip())

        # find files in cif which are not in comparison cif; then find those files in comparison not in cif
        missing_from_comparison = [x for x in files if x not in cifdf2cmpfiles]
        missing_from_current = [x for x in cifdf2cmpfiles if x not in files]

        # find changes for files that appear in both cif and cif2cmp
        diff_dict['Changes'] = dict()
        if missing_from_comparison:
          diff_dict['MISSING_FROM_Comparison_CIF'] = missing_from_comparison
        else:
            diff_dict['MISSING_FROM_Comparison_CIF'] = ['All files in {0} found in {1}'.format(self.cif, caldb2cmp.cif)]
        if missing_from_current:
            diff_dict['MISSING_FROM_CIF'] = missing_from_current
        else:
            diff_dict['MISSING_FROM_CIF'] = ['No files in {1} missing from {0}'.format(self.cif, caldb2cmp.cif)]

        # check changes in number of extensions, cal_qual
        # for each file in the CIF_to_CMP, find all the entries for that file in the CIF

        # first get list of files which are in both cifs by doing an inner merge

        filescmp = list(set(cifdf2cmp['CAL_FILE']))  # get a list of all unique file names from comparison CIF

        # for files that are in both the CIF and the comparison CIF, see if any values have changed
        for f in filescmp:
            # get rows from cif data frame where file = filename being checked
            testcif = cifdf[['CAL_FILE','CAL_XNO', 'CAL_QUAL']][cifdf['CAL_FILE']==f]
            # get rows from comparison cif data frame where file = filename being checked
            testcif2cmp = cifdf2cmp[['CAL_FILE','CAL_XNO', 'CAL_QUAL']][cifdf2cmp['CAL_FILE']==f]
            testcif.reset_index(inplace=True)
            testcif2cmp.reset_index(inplace=True)
            # then see if each entry in testcif2cmp is the same as the corresponding entry in testcif
            try:
                cols2cmp = testcif2cmp.columns
            except:
                print f, len(testcif), len(testcif2cmp), testcif2cmp.columns
            for i in range(len(testcif2cmp)):
                try:
                    ind = np.where(testcif.iloc[i]<>testcif2cmp.iloc[i])
                except:
                    print "For file {2}: Number of rows in testcif = {0}; number of rows in testcif2cmp ={1}".format(len(testcif), len(testcif2cmp), f)
                try:
                    val_diff = ind[0][0]
                except:
                    val_diff = None
                if val_diff:
                    # print "File {0} changed values:".format(f)
                    # print "Is"
                    cols=cols2cmp.values[val_diff]
                    # print testcif[cols]
                    # print "Was"
                    # print testcif2cmp[cols]
                    changedname = "{0}_{1}".format(testcif2cmp['CAL_FILE'].iloc[i], testcif2cmp['CAL_XNO'].iloc[i])
                    if type(cols)==str:
                        changedstring = "{0} Was {1} Now {2} ".format(cols, testcif2cmp[cols].iloc[i], testcif[cols].iloc[i])
                    else:
                        changedstring = ["{0} Was {1} Now {2} ".format(c, testcif2cmp[c].iloc[i], testcif[c].iloc[i]) for c in cols]
                    diff_dict['Changes'][changedname] = changedstring
        if print_report:
            print "CIF = {0}".format(diff_dict['CIF'])
            print "Comparison CIF ={0}".format(diff_dict['ComparisonCIF'])
            print "Files in CIF missing from comparison CIF:"
            for m in diff_dict['MISSING_FROM_Comparison_CIF']:
                print "   {0}".format(m)
            print "Files in Comparison CIF not in CIF:"
            for m in diff_dict['MISSING_FROM_CIF']:
                print "   {0}".format(m)
            print "Changes from Comparison CIF to CIF:"
            for c in diff_dict['Changes'].keys():
                print "   {0}".format(diff_dict['Changes'][c])
        return diff_dict

    def html_summary(self, missionurl='', outdir=".",
                     templatedir="/software/github/heasarc/pycaldb/templates",
                     cal_qual=5, clobber=False):
        """ write an html summary of the specified version of the Caldb instance

        NOTE: the instrument attribute of the CALDB instance should be set to the subdirectory
        of $CALDB/data/<telescope> in which the
        caldb index files are located
        (i.e., for NuSTAR, self.instrument = 'fpm' rather than 'fpma')

        :param missionurl: URL of the mission/telescope home page
        :param outdir: location of output html version of cif
        :param templatedir: location of FLASK template
        :param cal_qual: entries in cif with CAL_QUAL < this value will be displayed in the html file
        :param clobber: if True, overwrite existing file
        :return:

        CHANGES
        20170614 MFC - fixed error in call to cif_to_df; fixed templatedir and outdir on heasarcdev
        20170630 MFC - fixed error in definition of templatedir and outdir when running on heasarcdev

        """
        # status = make_missioncif_page(self.telescope, self.instrument, self.version, missionurl=missionurl,
        #                      outdir=outdir,
        #                      templatedir=templatedir,
        #                      cal_qual=cal_qual, clobber=clobber)

        status = 0
        try:
            version = self.version
        except AttributeError:
            print "Version needs to be explicitly set via the caldb instance's set_cif(version=version) method; returning"
            status = -1
            return status
        if not version:
            print "Version needs to be explicitly set via the caldb instance's set_cif(version=version) method; returning"
            status = -1
            return status
        if 'heasarcdev' in os.environ['HOST']:
            #
            # set directories appropriately if running on heasarcdev
            #
            templatedir = "/Home/home1/corcoran/Python/heasarc/pycaldb/templates"
            if not outdir:
                outdir = "/www/htdocs/docs/heasarc/caldb/data/{0}/{1}/index".format(self.telescope,self.instrument)
        try:
            ins_subdir = self.get_insdir()
        except:
            print "Error: Could not get subdirectory for instrument; returning"
            return
        instrument = os.path.split(ins_subdir)[1]
        cifdf = self.get_cif()
        if not 'FILTER' in cifdf.columns:
            print ''
            ans = raw_input('CIF missing FILTER values; set FILTER to "INDEF" ([y]/n) > ')
            if ans.strip().lower() <> 'n':
                cifdf['FILTER'] = "INDEF"
        goodcifdf = cifdf[cifdf['CAL_QUAL'] <> cal_qual]
        goodcifdf.reset_index(inplace=True)  # resets the index to be 0 -> len(goodcifdf)
        # Render html file
        try:
            templateLoader = jinja2.FileSystemLoader(searchpath=templatedir)
        except Exception, errmsg:
            print ('Could not load templates in {0}; Returning ({1})'.format(templatedir, errmsg))
            status = -1
            return status
        templateEnv = jinja2.Environment(loader=templateLoader, trim_blocks=True, )

        # define filter for DETNAM since DETNAM not included in chandra cifs
        # see http://jinja.pocoo.org/docs/2.9/api/#writing-filters
        def set_detnam(row):
            try:
                detnam = row.DETNAM
            except AttributeError:
                detnam = 'INDEF'
            return detnam

        templateEnv.filters['set_detnam'] = set_detnam
        template = templateEnv.get_template('missioncif_table_template.html')
        output_html = template.render(cifdf=goodcifdf, telescope=self.telescope.upper(),
                                      instrument=instrument, version=self.version,
                                      cifname=self.cif, missionurl=missionurl)
        outname = "cif_{0}_{1}_{2}.html".format(self.telescope, instrument,self.version.replace('caldb.indx',''))
        outname=outname.replace('.indx','') # needed for chandra

        fname = '{0}/{1}'.format(outdir, outname)
        if os.path.isfile(fname):
            if not clobber:
                print("File {0} exists, and clobber = False; not overwritten".format(fname))
                return
        try:
            with open(fname, mode='w') as fout:
                fout.write(output_html)
        except IOError, errmsg:
            print('Could not write {0}; Returning ({1})'.format(fname, errmsg))
            status = -1
            return status
        print "Wrote {0}".format(fname)
        return fname

    def create_caldb_tar(tarName="",
                         tardir=".", calQual=0,
                         caldbdirs=['bcf', 'cpf'],
                         caldb="http://heasarc.gsfc.nasa.gov/FTP/caldb"):
        """ Create CALDB tar file

        creates a tar file of "good" files for given mission/instrument
        from the specified caldb.

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
        :param caldbdirs: subdirectories of $CALDB/data/<telescop>/<instrume> to tar
            (the index directory is always retrieved)
        :param caldb: top-level caldb directory or url
        :return: status (0 if no errors)
        """
        import tarfile
        from ftplib import FTP
        import shutil
        tel = self.telescope
        instr = self.instrume
        # ver should be of the form 20170503
        ver = version.split('caldb.indx')[1].strip()
        cwd = os.getcwd()
        status = 0
        calQual = int(calQual)  # make sure this is an integer

        cif = self.cif
        try:
            hdu = pyfits.open(cif)
        except Exception, errmsg:
            print("Could not open {0}; returning".format(cif))
            status = -1
            return status
        if not tarName:  # if tarname not specified, use "standard" naming
            tarName = "{3}/goodfiles_{0}_{1}_{2}.tar.gz".format(tel, instr, ver, tardir)
        tar = tarfile.open(tarName, "w:gz")
        if (("http://" in caldb) or ("ftp://" in caldb)):  # if remote, download the index file directory
            tindname = "{0}/index.tar".format(tardir)
            tindex = open(tindname, 'w')
            ftp = FTP("heasarc.gsfc.nasa.gov")
            ftp.login(user="anonymous", passwd="")
            # heasarc ftp allows retrieval of an entire directory if .tar appended to directory name
            ftp.retrbinary("RETR /FTP/caldb/data/{0}/{1}/index.tar".format(tel, instr), tindex.write)
            tindex.close()
            t = tarfile.open(tindname, 'r')
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
            shutil.rmtree(tardir + '/index')  # remove local index directory
            os.remove(tardir + '/index.tar')
        else:
            cifdir = "{0}/data/{1}/{2}/index".format(caldb, tel, instr)
            cifs = glob.glob("{0}/caldb*".format(cifdir))
        arcname = ['index/' + c[c.rfind('/') + 1:] for c in cifs]
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
        tar.add("caldb.indx", "data/{0}/{1}/caldb.indx".format(tel, instr))
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
            testfile = "{0}/{1}".format(cd, cf)
            if (cq == calQual) and (testfile not in cfile_to_tar):
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
            if ("http://" in caldb) or ("ftp://" in caldb):  # if remote access, download the file
                # download the goodfile from the heasarc
                response = urllib2.urlopen("{0}/{1}".format(caldb, goodfile))
                data = response.read()
                localgoodfile = goodfile[goodfile.rfind("/") + 1:]  # strip off the directory path for the local file
                tmp_cfile = "{0}/{1}".format(tardir, localgoodfile)
                out = open(tmp_cfile, 'w')
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
            else:  # using locally-mounted caldb
                # tmp_cfile = "{0}/{2}".format(caldb, cd, cf)
                cfile_to_add = "{0}/{1}".format(caldb, goodfile)
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


class CaldbFile(object):
    """
    A caldb file object has these attributes:
        - a file name
        - a caldb path relative to $CALDB (note that the actual subdirectory holding the file is not specified by the file information)
        - a detector to which it belongs (optional)
        - an instrument to which the detector belongs
        - a mission to which the instrument belongs
        - a number of elements (primary HDU + extensions)
        - a file size in bytes
        - a time of creation
        - a time of last modification
        - a caldb class (CCLSmmmmm) (bcf, cpf etc)
        - each element (primary hdu + extensions) of the file has:
        - a caldb contents descriptor label (CCNMmmmm)
        - A caldb validity start date (CVSDmmmmm keyword),
        - a caldb validity start time on the start date (CVSTmmmmm)
        - a caldb data type (CDTPmmmmm)
        - zero or more calibration boundary limits (CBDnmmmm)
        - a description of the element contents (CDESmmmm)

    """
    def __init__(self, filename):
        """
        :param filename: the name of the file (with directory path if the file is not in the local directory)
        """
        self.file_size = os.path.getsize(filename)
        self.mod_time = time.ctime(os.path.getmtime(filename))
        self.create_time = time.ctime(os.path.getctime(filename))
        cdu = pyfits.open(filename)
        self.num_elements=(len(cdu)) # don't count primary array as an extension
        telescope = cdu[0].header['TELESCOP']
        instrument = cdu[0].header['INSTRUME']
        self.telescop = telescope
        self.instrume = instrument
        calinfo = [] # calinfo is an array of dictionariers which summarizes the caldb file elements
        self.filename = filename
        self.path= os.path.split(filename)[0]
        cdict = dict()
        for i, c in enumerate(cdu):
            cdict = get_calkeys(filename, i)
            calinfo.append(cdict)
        self.calinfo = calinfo
        try:
            calclass =  calinfo[1]['CCLS0001']
            self.path = os.path.join('data', telescope, instrument, calclass)
        except Exception, errmsg:
            print "Required CALDB Keyword CCLS0001 missing from first extension; cannot create CALDB path"
            self.path = None
        return

    def summarize(self):
        """
        prints out a summary view of the information stored in a calfile
        :param self:
        :return:
        """
        for i, cc in enumerate(self.calinfo):
            numbds = len(cc['CBD'])
            print("Ext #{0} CCNM: {1} CVSD: {2} # CBDs: {3}".format(i, cc['CCNM0001'], cc['CVSD0001'], numbds))
            if numbds:
                print "    Boundaries = ",
                for cb in cc['CBD'].values():
                    print " {0}".format(cb),
                print "\n",
        return

def file_exists(self):
    """
    returns true if file exists in the caldb
    :return:
    """
    if os.sep in self.filename:
        file = os.path.join(self.caldb,self.filename)
    else:
        # TODO: get path to file from the caldbindx file
        pass
    return os.path.isfile(file)

def get_cbds(telescope, instrument, caldb=None):
    """ Get calibration boundaries from a caldb file

    Accesses the calibration boundary values for a specified telescope/instrument

    :return: calibration boundary values and table of CIF data
    """
    # print cal_cbd
    # urllib.urlcleanup() # this apparently removes the file
    if not caldb:
        caldb = os.environ['CALDB']
    hdulist = get_cif(telescope, instrument, caldb=caldb)
    tbdata = hdulist[1].data
    cbds = tbdata.field('CAL_CBD')
    # nel=len(cbds)
    # return array of boundaries
    cbds = cbds.split()  # split the cbd values on whitespace
    print "Splitting cbds"
    return (cbds, tbdata)


def get_caldbconfig(caldbconfig):
    """
    Parses the caldb.config file, returns a dictionary of:
    caldbconfig={'telescope':{'instrument1':['caldb_variable', 'caldb_index_directory','cifname','store_dev','store_dir']}}
    From the caldb.config file:
    # Each non-comment line has seven tokens.  The first two tokens are the
    # mission and instrument (these tokens must be entirely in uppercase);
    # the remaining tokens describe the location of the index file and the
    # calibration storage directory in a system independent manner.  Tokens
    # 3, 4, and 5 are the device (can also be an environment variable),
    # directory (subdirectories are separated by '/' characters), and file
    # name of the Calibration index file.  Tokens 6 and 7 are the device and
    # directory of the calibration storage directory, and are not used;
    # they are available for historical reasons only
    #
    #
    # For example, a typical line in this file might look like:
    #
    # HOSR HOSEHD caldb data/hosr/hosehd caldb.indx caldb data/hosr/hosehd
    #
    # If this file were read on a UNIX machine, then this line would
    # indicate that data for the instrument HOSEHD, of the mission HOSR will
    # be indexed in the file
    #
    # /caldb/data/hosr/hosehd/caldb.indx
    #
    # and the calibration files for this instrument will be located beneath
    # the directory
    #
    # /caldb/data/hosr/hosehd

    :param caldbconfig: path to/name of the caldb config file
    :return: caldbconfig dictionary
    """
    try:
        with open(caldbconfig, 'r') as f:
            ccon = f.readlines()
    except IOError, errmsg:
        print "Problem opening {0} ({1}); returning".format(caldbconfig, errmsg)
        return -1
    ccon = [x.strip() for x in ccon if '#' not in x] #remove comments and strip newlines
    ccon = [x.split() for x in ccon]
    cfig=dict()
    # Create an dictionary of all the telescopes
    for c in ccon:
        try:
            cfig[c[0].lower().strip()]=dict()
        except IndexError, errmsg:
            print "Could not parse TELESCOP from caldb config line = {0}".format(c)
            pass
    # for each telescope dictionary create subdictionary of instrument info
    for c in ccon:
        try:
            cfig[c[0].lower().strip()][c[1].lower().strip()] = c[2:]
        except IndexError, errmsg:
            print "Could not parse instrument from caldb config line = {0}".format(c)
    return cfig

def read_update_notice(telescope,updatedate,notice_dir = '/Users/corcoran/Desktop/caldb_updates',
                       chatter = 0):
    """
    reads the update notice

    reads the update notification e-mail text file and returns a dictionary of telescope, instrument,
    and list of files for the update.  See
    
    "Automated Delivery of Calibration Data to the CALDB"
    (https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_2003_001/cal_gen_2003_001.html)
    for specification on the format of the update notification e-mail
    
    :param telesope: caldb telescope name
    :param updatedate: date of update in YYYYMMDD format as specified in the subject of the e-mail; not to be confused with the update VERSION, also in YYYYMMDD format.  The update VERSION gives the version number of the individual caldb.indx file describing the index, for example the string 20150304 in caldb.indx20150304
    :param notice_dir: directory where the notification e-mail is located
    :return: update_dict, a dictionary of instruments and files for the specified update

    VERSIONS
    0.1 MFC 20170614 Initial Version

    """
    telescope = telescope.strip().lower()
    if 'heasarcdev'in os.environ['HOST']:
        notice_dir = '/FTP/caldb/staging/data/{0}/to_ingest'.format(telescope)
    update_email = '{0}/caldb_update_{1}_{2}.txt'.format(notice_dir,telescope, updatedate)
    update_dict = dict()
    update_dict[telescope] = dict()
    update_dict[telescope][updatedate]=dict()
    if chatter > 0:
        print("Opening update e-mail {0}".format(update_email))
    try:
        with open(update_email, mode='r') as f:
            updatelist = f.readlines()
            if chatter > 0:
                print updatelist
    except IOError, errmsg:
        print ("Problem opening {0} ({1}; returning".format(update_email, errmsg))
        return
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

def read_cif(telescope=None,instrument=None, version=None, cifname = None, caldb=None):
    """ convert caldb.indx file to FITS HDU
    
    reads the caldb.indx file for a given telescope and instrument, for either a local caldb installation
    or a remote installation access via http or ftp
    
    :param telescope: telescope specification for cif
    :param instrument: instrument specification for telescope
    :param version: version of cif to use (for example caldb.indx20160502)
    :return: returns a FITS hdu representation of the caldb.indx file
    """
    if not cifname:
        try:
            cifname = mk_cifname(telescope, instrument, version=version, caldb=caldb)
        except:
            print "Error in running mk_cifname in get_cif; returning"
            return -1
    if not caldb:
        try:
            caldb = os.environ['CALDB']
        except:
            print 'CALDB not specfied and $CALDB environment variable not set; returning'
            return -1
    try:
        hdulist = pyfits.open(cifname)
    except Exception, errmsg:
        print ('Could not open file {0} ({1})'.format(cifname,errmsg))
        return -1
    return hdulist

def cif_to_df(telescope, instrument, version=None, cifname = None):
    """ convert caldb.indx file to a dataframe
    
    Creates a pandas dataframe from a calibration index file (cif) for a given mission/instrument
    
    :param telescope: caldb standard telescope name
    :param instrument: caldb standard instrument name for that telescope
    :param cif: (optional) specific cif to retrieve
    :return: pandas dataframe containing the data from the 1st extension in the CIF
    """
    if not cifname:
        cifname = mk_cifname(telescope, instrument, version)
    print "Retrieving CIFNAME = {0}".format(cifname)
    cif = get_cif(cifname=cifname)
    try:
        cifdf = Table(cif[1].data).to_pandas()
    except TypeError, errmsg:
        sys.exit("Error accessing CIF {0}; exiting".format(cif))
    if 'INSTRUME' not in cifdf.columns:
        cifdf['INSTRUME'] = 'INDEF'
    # strip whitespace from string columns
    stripcols = cifdf.columns
    for s in stripcols:
        try:
            cifdf[s] =  [x.strip() for x in cifdf[s]]
        except:
            pass
    return cifdf

def cifstats(telescope='', instrument='', version='',
             caldb='http://heasarc.gsfc.nasa.gov/FTP/caldb', cifname = None):
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
    cifdf = cif_to_df(telescope, instrument, version=version, caldb=caldb, cifname=cifname)
    cifstat= dict()
    cifstat['Number_Unique_filenames'] = len(set(cifdf.CAL_FILE))
    cifstat['Unique_Instrument_names']= [x.strip() for x in list(set(cifdf.INSTRUME))]
    calfiles = []
    for cdir, cfile in zip(cifdf.CAL_DIR,cifdf.CAL_FILE):
        calfiles.append("{0}/{1}".format(cdir, cfile))
    cifstat['Number_Unique_files']=len(set(calfiles))
    cifstat['Number_files_qual_0'] = len(set(cifdf[cifdf.CAL_QUAL == 0].CAL_FILE))
    cifstat['Number_files_qual_5'] = len(set(cifdf[cifdf.CAL_QUAL == 5].CAL_FILE))
    cifstat['Number_of_files_other_quality'] = cifstat['Number_Unique_filenames'] - (cifstat['Number_files_qual_0']+cifstat['Number_files_qual_5'])
    cifstat['Number_unique_CNAMS']=len(set(cifdf.CAL_CNAM))
    return cifstat


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

def get_calkeys(filename, extno):
    """
    Get the caldb keywords from the given extension number in a caldb file
    :param filename: name of the caldb file
    :param extno: integer extension number to search (0 for the primary HDU)
    :return: dictionary of the caldb keyword values
    """
    calkeys = ['CCLS0001',
               'CDTP0001',
               'CCNM0001',
               'CVSD0001',
               'CVST0001',
               'CDES0001']
    caldict=dict()
    cdu = pyfits.open(filename)
    for ck in calkeys:
        try:
            caldict[ck]=cdu[extno].header[ck]
        except KeyError, errmsg:
            print ("{0} File = {1} Ext. = {2}".format(errmsg, filename, extno))
            caldict[ck] =''
    # get the set of calibration boundary keywords as a pyfits header object
    caldict['CBD'] = cdu[extno].header['CBD*']
    return caldict

def get_cbds(telescope, instrument, caldb='http://heasarc.gsfc.nasa.gov/FTP/caldb'):
    """ Get calibration boundaries
    
    Accesses the calibration boundary values in the caldb.indx file
    for a specified telescope/instrument
    
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


def ck_calfile_incif(calfile, telescope, instrument,
                     caldb='http://heasarc.gsfc.nasa.gov/FTP/caldb',
                     version='', verbose=False):
    """ verifies that a specified file is actually in the caldb.indx file

    for a specified filename, checks that this file exists in the specified calibration index file (cif),
    then checks that the CALDB keyword values are consistent between the cif and the actual file
    determines if the file exists in the caldb

    usage:
    ck_calfile_incif('data/hitomi/sxi/cpf/background','ah_sxi_nxbpnt_20140101v001.pha',
    caldb = 'http://heasarc.gsfc.nasa.gov/FTP/caldb',cif='caldb.indx20170405')

    :param calfile: caldb file relative to $CALDB
    :param telescope: Name of TELESCOP
    :param instrument: Name of Instrument
    :param caldb: caldb top level directory/url
    :param version: version of cif to use (for example caldb.indx20170405); if blank use caldb.indx
    :return:
    """
    # get telescope and instrument from the caldir
    if not version:
        cifhdu = get_cif(telescope, instrument, caldb=caldb)
    else:
        cif = get_cif(telescope, instrument, version=version, caldb=caldb)
    cfiles = cifhdu[1].data['CAL_FILE']
    cal_dir = cifhdu[1].data['CAL_DIR']
    # create list of files in cif with their path relative to $CALDB
    ciffiles = ["{0}/{1}".format(x[0], x[1]) for x in zip(cal_dir, cfiles)]
    if calfile in ciffiles:
        if verbose:
            print "Found {0} in caldb.indx".format(calfile)
        status = True
    else:
        if verbose:
            print "{0} not in caldb.indx".format(calfile)
        status = False
    return status


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


def sync_caldbweb(mission, instrument, version, DoCopy=False, user='mcorcora', host="heasarcdev",
                  htmlcaldbdev='/www/htdocs/docs/heasarc/caldb/data',
                  htmlcaldbprod= '/www.prod/htdocs/docs/heasarc/caldb/data', prompt=True):
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
    :param instrument: CALDB-standard instrument name
    :param version: version of caldb (usually of form YYYYMMDD, except for Chandra)
    :param DoCopy: if False just print the cmds that will be executed when DoCopy = True
    :param Prompt: if True will prompt the user to continue if DoCopy
    :return: returns a list of error messages in case of partially successful copies, or -1 if user quits

    CHANGES
    20170614 added Prompt keyword if DoCopy chose; also fixed file list for nustar, added instrument parameter
    """
    import os
    mission = mission.lower().strip()
    errlist = []
    curuser = os.environ['LOGNAME']
    curhost = os.environ['HOSTNAME']
    if (curuser <> user) or (host not in curhost):
        error = "This function must be run as {0} from {1}; returning".format(user, host)
        print(error)
        print("Current user is:{0}  Current host is: {1}".format(curuser, curhost))
        errlist.append(error)
        return errlist
    # generate the caldb news item
    print "# Generating CALDB news item (develop version)"
    cmd = "~{0}/bin/caldbdev_rss2.pl".format(user)
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
        ("caldb_wotsnew.html", "/www/htdocs/docs/heasarc/caldb", "/www.prod/htdocs/docs/heasarc/caldb"),
        ("caldb_supported_missions.html", "/www/htdocs/docs/heasarc/caldb", "/www.prod/htdocs/docs/heasarc/caldb"),
        ("caldb.rss", "/www/htdocs/docs/heasarc/caldb", "/www.prod/htdocs/docs/heasarc/caldb"),
        ("associated_data.html", "/www/htdocs/docs/heasarc/astro-update/inc",
         "/www.prod/htdocs/docs/heasarc/astro-update/inc"),
        ("astro-update.rss", "/www/htdocs/docs/heasarc/astro-update", "/www.prod/htdocs/docs/heasarc/astro-update")
    ]
    cifhtmldev = "{0}/{1}/{2}/index".format(htmlcaldbdev, mission.strip().lower(), instrument.strip().lower())
    cifhtmlprod = "{0}/{1}/{2}/index".format(htmlcaldbprod, mission.strip().lower(), instrument.strip().lower())
    files_to_sync.append(("cif_nustar_fpm_{0}.html".format(version),
                          cifhtmldev, cifhtmlprod))
    if mission == 'nustar':
        files_to_sync.append(("release_{0}.txt".format(version), "/www/htdocs/docs/heasarc/caldb/nustar/docs/",
                             "/www.prod/htdocs/docs/heasarc/caldb/nustar/docs/"))
        files_to_sync.append(("nustar_caldbhistory.html", "/www/htdocs/docs/heasarc/caldb/nustar/docs",
                             "/www.prod/htdocs/docs/heasarc/caldb/nustar/docs"))
        files_to_sync.append(("nustar_caldb.html", "/www/htdocs/docs/heasarc/caldb/nustar",
                             "/www.prod/htdocs/docs/heasarc/caldb/nustar"))
        files_to_sync.append(("nustar.rss", "/www/htdocs/docs/nustar/news/",
                             "/www.prod/htdocs/docs/nustar/news/"))
    print "\n # Syncing website:"
    for f in files_to_sync:
        print "# copying {0} from {1} to {2}".format(f[0], f[1], f[2])
        cmd = "cp {1}/{0} {2}/{0}".format(f[0], f[1], f[2])
        print "{0}\n".format(cmd)
        if DoCopy:
            status = os.system(cmd)
            if prompt:
                ans = hu.yesno('Continue?')
                if ans:
                    pass
                else:
                    print "Stopping at user request"
                    status = -1
                    return status
            if status <> 0:
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
            if prompt:
                ans = hu.yesno('Continue?')
                if ans:
                    pass
                else:
                    print "Stopping at user request"
                    status = -1
                    return status
            if status <> 0:
                errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    # create heasarc news feed from astroupdate
    print "\n# Creating heasarc newsfeed (dev version)".format(mission)
    cmd = "/www/htdocs/docs/rss/software_rss2.pl"
    print("{0}\n".format(cmd))
    if DoCopy:
        status = os.system(cmd)
        if status <> 0:
            errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    print "\n# Creating heasarc newsfeed (public version)"
    cmd = "/www.prod/htdocs/docs/rss/software_rss2.pl"
    print("{0}\n".format(cmd))
    if DoCopy:
        status = os.system(cmd)
        if prompt:
            ans = hu.yesno('Continue?')
            if ans:
                pass
            else:
                print "Stopping at user request"
                status = -1
                return status
        if status <> 0:
            errlist.append('Problem executing {0}; status={1}'.format(cmd, status))
    return errlist

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
                     tardir = "", calQual = 0,
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
    tel=telescop.strip()
    instr = instrume.strip()
    #ver = version.strip()
    # ver should be of the form 20170503
    ver = version.split('caldb.indx')[1].strip()
    cwd = os.getcwd()
    status = 0
    calQual = int(calQual) # make sure this is an integer
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
        if (cq == calQual) and (testfile not in cfile_to_tar):
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


def test_makemissioncif(telescope='Swift',instrument='SC', version='20170629',missionurl='http://swift.gsfc.nasa.gov'):
    outdir="/software/github/heasarc/pycaldb/html/cifs"
    make_missioncif_page(telescope, instrument,version, missionurl=missionurl, clobber=False, outdir=outdir)
    return

def check_cif_diff():
    fermiv10 = Caldb(telescope='fermi', instrument='lat')
    fermiv10.set_cif('caldb_v10r0p5_20150509.indx')
    fermiv11 = Caldb(telescope='fermi', instrument='lat', caldb='/Volumes/SXDC/caldb_test/fermi')
    fermiv11.set_cif('caldb_v11r5p3_20180215.indx')
    diff_dict = fermiv11.cif_diff(fermiv10)
    return

    
if __name__ == "__main__":# create_caldb_tar('nustar','fpm','20161021', tardir='/Users/corcoran/Desktop/tmp/caltartest',
    #                 tarName='goodfiles_nustar_fpm_94nov9.tar.gz',
    #                 caldb= '/caldb')
    # test
    # create_caldb_tar('swift', 'xrt', '20160609 ', tardir='/FTP/caldb/staging/tmp', caldb='/FTP/caldb')
    # test_pycaldb(0)
    #
    # test_makemissioncif()
    #
    # cifdf = cif_to_df('nustar','fpm')
    # cstats = cifstats('nustar','fpm')
    # cstatsold = cifstats('nustar','fpm', version='caldb.indx20160606')
    # for k in cstats.keys():
    #     print "{0:30s} latest = {1} 20160606 = {2}".format(k, cstats[k],cstatsold[k])
    # ud = read_update_notice('swift', '20081203')
    # for k in ud['swift'].keys():
    #    for i in ud['swift'][k].keys():
    #        print k, i, ud['swift'][k][i]
    # ud = read_update_notice('Swift', '20170629')
    #
    # cfg = get_caldbconfig(cfig='caldb.config', caldb='/caldb')
    # for k in cfg.keys():
    #     print k, cfg[k].keys()

    #f = '/caldb/data/swift/xrt/bcf/instrument/swxhkranges0_20010101v008.fits'
    #cd = get_calkeys(f,1)
    #print cd
    #
    # debug cif_diff
    # cifold = Caldbindx('swift', 'xrt', cifname='/caldb/data/swift/xrt/index/caldb.indx20160121')
    # cifnew = Caldbindx('swift', 'xrt', cifname='/caldb/data/swift/xrt/index/caldb.indx20160609')
    # cdiff = cifnew.cif_diff(cifold)
    # for k in cdiff.keys():
    #     print "{0} => {1}".format(k, cdiff[k])

    # check Caldbindx object creation
    # os.environ['CALDBCONFIG'] = '/software/caldb/caldb.config'
    # os.environ['CALDBALIAS'] = '/software/caldb/alias_config.fits'
    # os.environ['CALDB'] = 'ftp://heasarc.gsfc.nasa.gov/caldb'
    # nucif = Caldb()
    # nucif.set_telescope('nustar')
    # nucif.set_instrument('fpma')
    # nucif.set_cif()
    # diff_dict = nucif.cif_diff(version='caldb.indx20170727')
    # for k in diff_dict.keys():
    #     print k, diff_dict[k]

    # check find_calfiles method
    # caldb = Caldb(telescope='asca', instrument='sis0')
    # caldb.set_cif()
    # caldb.find_calfiles('FULL',)

    # uvcaldb=Caldb(telescope='swift',instrument='uvota')
    # uvcaldb.set_cif(version='caldb.indx20170922')
    # outdir='/software/github/heasarc/pycaldb/html/cifs'
    # uvcaldb.html_summary(missionurl='https://swift.gsfc.nasa.gov', clobber=True, outdir=outdir)

    #check_cif_diff()

    testfile = CaldbFile('/Users/corcoran/program/HEASARC/missions/NICER_HEASARC/calibration/caldbdev/20180404/nixtisoyuz20170601v001.fits')



