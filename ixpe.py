def ixpe_caldbupdate(tarName, verbose=True, 
                        scratchdir = '/Users/mcorcora/software/caldbtest/data/ixpe/tars',
                        caldb = '/Users/mcorcora/software/caldbtest/',
                        mission='ixpe'):
    """
    updates the ixpe caldb from an update tar file (usually from Kurt Dietz)
       * checks that all the files in an ixpe caldbupdate are unique
       * checks that all files have checksums
       * ftverifies the files

       if those tests pass then untars the update tar file into the ixpe caldb and creates
       the goodfiles tar file and also an update tar file (smaller) for the gpd and xrt
    :param tarName: name of tarfile containing the files in the IXPE caldb update (eg ixpe_caldb_update_20230526_01.tgz)

    Note that for ixpe, the file in the tar file are relative to the ixpe/ directory
    i.e. ixpe/xrt/bcf/reef/ixpe_m1_20210103_reef_03.fits
    """
    import tarfile
    import heasoftpy as hsp
    import os
    from heasarc.pycaldb import pycaldb as pc

    status = 0

    try:
        tf = tarfile.open(tarName)
    except Exception as e:
        print(f'Problem opening {tarName}\n ({e});\n Returning')
        return -1
    # get filenames (without including directory names or 'caldb.indx' files)
    files = [m for m in tf.getnames() if (os.path.split(m)[-1] != 'caldb.indx') and ('.' in m)]
    # extract the files from the update tar file to scratchdir
 
    tf.extractall(path=scratchdir)
    dupfiles = []
    for f in files:
        cfile = os.path.join(caldb,'data',f)
        if os.path.exists(cfile):
            if verbose:
                print(f'\nDUPLICATE: {cfile}', end='\n\n')
            dupfiles.append(f)
    if len(dupfiles) > 0:
        return {'DUPLICATES':dupfiles}
    
    # no duplicates found so check that all files/extensions have checksums
    print('No duplicate files found')

    missingcksum = []
    for f in files:
        cksumout = hsp.ftchecksum(infile = f)
        if 'missing' in cksumout.stdout:
            if verbose:
                print(f'{f} MISSING CHECKSUM')
            missingcksum.append(f)
        if len(missingcksum)>0:
            return {'MISSING_CKSUM':missingcksum}

    print('All extensions have checksums')

    # all files/extensions have checksums, so ftverify the files which have been untarred in scratchdir
    noverif = []

    print('FTVERIFYING')   
    for f in files:
        cfile = os.path.join(scratchdir,f)
        out = hsp.ftverify(infile = cfile, prstat='no', numerrs=0, numwrns=0)
        if out.stdout.split(':')[0] != 'verification OK':
            if verbose:
                print(f'{cfile} FAILED Verification')
                print(out, end='\n\n')
            noverif.append(f)
        if len(noverif) > 0:
            return {'VERIFICATION_FAIL':noverif}
    
    # all files verify, so untar into the ixpe caldb
    print ('All files have Valid checksums and pass FITS verification; UNTARRING', end='\n\n')

    tf.extractall(path=os.path.join(caldb,'data'))
    
    #check that new files exist in caldb
    missingfromcaldb = []             
    
    for f in files:
        if not os.path.exists(os.path.join(caldb,'data',f)):
            if verbose:
                print(f'{f} MISSING from {caldb+"/data"}; Returning')
                missingfromcaldb.append(f)
    if len(missingfromcaldb) > 0:
        return {'MISSING':missingfromcaldb}
    
    print('All extracted tarfiles in update now in CALDB')
    # create update tarfile by cd'ing to $CALDB, then for a given instrument
    # find all the files for that instrument in the update tar file
    # then create a tarfile with only those files in it

    # get instruments in tarfile

    os.chdir(caldb)
    inst = list(set([f.split('/')[1] for f in files]))

    # create update tarfile for all instruments
    for i in inst:
        # get current version of caldb for instrument
        indx = [x for x in files if f'{i}/index/' in x]
        try:
            version = indx[0].replace(f'ixpe/{i}/index/caldb.indx','')
        except Exception as e:
            print(f'Could not get version for {i} ({e}); returning')
            return 1
        
        #create name of output tarfile
        taroutfile = os.path.join(caldb,'data','ixpe',i,f'goodfiles_{mission}_{i}_{version}_update.tar.gz')
        print(f'Creating {taroutfile}')

       # create tarfile object
        tfout = tarfile.open(taroutfile,'w:gz')
        ifiles = [f for f in files if i in f]

        # add files to tarfile
        for ifi in ifiles:
            tfout.add(os.path.join('data',ifi))

        # add caldb.indx file
        tfout.add(os.path.join('data',mission,i,'caldb.indx'))

        #write tarfile
        tfout.close()

    # create tarfile of all good files for each instrument

    for i in inst:
        print('Instrument = ',i,end='\n\n')
        Caldb = pc.Caldb(mission, i)
        Caldb.set_cif(f'caldb.indx{version}')
        goodcal = Caldb.get_cif(return_table=True,calqual=0)
        taroutfile = os.path.join(caldb,'data','ixpe',i,f'goodfiles_{mission}_{i}_{version}.tar')
        print(f'Creating {taroutfile}')
        cmd = f'tar -cvf {taroutfile} data/{mission}/{i}/index/caldb.indx{version}'
        print(cmd)        
        os.system(cmd)
        cmd = f'tar -rf {taroutfile} data/{mission}/{i}/caldb.indx'
        print(cmd)        
        os.system(cmd)
        j = -1
        #create list of calfiles
        calfiles = []
        for x,y in zip(goodcal['CAL_DIR'],goodcal['CAL_FILE']):
            j += 1
            cfile = f'{x.strip()}/{y.strip()}'
            if not os.path.exists(cfile):
                print('{calfile} MISSING from CALDB AFTER UPDATE')
            else:
                calfiles.append(cfile)
        # get unique list of calfiles
        calfiles = list(set(calfiles))
        # add unique list of calfiles to tar file
        for calfile in calfiles:
            cmd = f'tar -rf {taroutfile} {calfile}'
            print(cmd)
            os.system(cmd)
        # create 
    return status


if __name__ == "__main__":
    import os
    print('Starting test')
    os.environ['CALDB'] = "/Users/mcorcora/software/caldbtest"
    os.chdir('/Users/mcorcora/software/caldbtest')
    cmd = 'source test_ixpecaldb.csh'
    print (cmd)
    os.system(cmd)
    tarName = '/Users/mcorcora/software/caldbtest/tars/ixpe_caldb_update_20230526_01.tgz'
    ixpe_caldbupdate(tarName, verbose=True, 
                        scratchdir = '/Users/mcorcora/software/caldbtest/scratch',
                        caldb = '/Users/mcorcora/software/caldbtest/',
                        mission='ixpe')

