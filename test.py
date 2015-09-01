ddir='/Users/corcoran/research/ETA_CAR/Swift/2014/work'
cdir=os.getcwd()
os.chdir(ddir)
#seq=['00091911078','00091911079']
seq=glob.glob('000*')
#seq=[0:2]
for i in seq:
    fname=i+'/xspec/ec_srcbin10.pha'
    try:
        hdr=pyfits.getheader(fname,ext=1)
        dat=pyfits.getdata(fname,ext=1)
        totcounts=' {0:.0f}'.format(sum(dat['COUNTS']))
        dobs=hdr['DATE-OBS']
        t=Time(dobs,form='iso',scale='utc')
        dend=hdr['DATE-END']
        te=Time(dend,form='iso',scale='utc')
        dur=(te.jd-t.jd)
        d=' {0:.3f}'.format(dur)
        tmid=dur/2.0+t.jd - 2400000.0
        tm=' {0:.3f}'.format(tmid)
        m=Time(tmid+2400000.0,format='jd',scale='utc')
        endt=Time(dend,form='iso',scale='utc')
        phi=' {0:.3f}'.format((tmid-(epoch-2400000.0))/period)
        expo=' {0:.0f}'.format(hdr['EXPOSURE'])
        print i+' '+m.iso+tm+d+phi+expo+totcounts+'\n'