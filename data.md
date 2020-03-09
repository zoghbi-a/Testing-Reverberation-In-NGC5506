Using:
- `heasoft 6.26`
- `xspec 12.10.1f`
- `xmmsas_20180620_1732-17.0.0`
- `ccf` update `28-05-2019`

# XMM Data

## Download the data
- Use xamin to create a list of all obsids


```python
import numpy as np
import os
import subprocess as subp
import glob
import time
import re
import aztools as az
from ftplib import FTP
from astropy.io import fits as pyfits
import astropy.time as atime
import sys

base_dir = '/home/abzoghbi/data/ngc5506/spec_timing_paper'
sys.path.append(base_dir)
from timing_helpers import *
%load_ext autoreload
%autoreload 2
```

    The autoreload extension is already loaded. To reload it, use:
      %reload_ext autoreload



```python
data_dir = 'data/xmm'
os.system('mkdir -p %s'%data_dir)
obsids = ['0013140101', '0201830201', '0201830401', '0554170101', '0761220101',
          '0013140201', '0201830301', '0201830501', '0554170201']
obsids = np.array(obsids)
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
```

    There are 9 observations
    0013140101, 0201830201, 0201830401, 0554170101, 0761220101, 0013140201, 0201830301, 0201830501, 0554170201


- We use `ftplib` to get the data from heasarc (may take some time)


```python
os.chdir(base_dir)
ftp = FTP('legacy.gsfc.nasa.gov', 'anonymous', 'anonymous@gmail.com')
ftp.cwd('xmm/data/rev0')
failed = []
for o in obsids:
    tar_file = '%s/%s.tar'%(data_dir, o)
    # download file only if not already downloaded
    if not os.path.exists(tar_file):
        try:
            ftp.retrbinary('RETR %s.tar'%o ,open(tar_file, 'wb').write)
        except:
            print('failed downloading %s'%o)
            os.system('rm %s >/dev/null 2>&1'%tar_file)
            failed.append(o)

```


```python
for f in failed:
    obsids = np.delete(obsids, np.argwhere(obsids==f)[0,0])
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
```

    There are 9 observations
    0013140101, 0201830201, 0201830401, 0554170101, 0761220101, 0013140201, 0201830301, 0201830501, 0554170201


## Process the PN data
We use our shell script `xmm_process`. Split it into two parts so we can run things in parallel across observations. The first creates `ccf.cif`, and the second creates the event files


```python

os.chdir('%s/%s'%(base_dir, data_dir))
os.system('mkdir -p log')
procs = []
for o in obsids:
    #if os.path.exists(o): continue
    #os.system('tar -xf %s.tar'%o)
    os.chdir(o)
    os.system('rm -r 3XMM om_mosaic PPS >/dev/null 2>&1')
    os.system('mv ODF odf')
    os.chdir('odf')
    if not os.path.exists('ccf.cif'):
        os.system('gzip -d *gz')
        log_file = '../../log/%s_process.log'%o
        proc = subp.Popen(['/bin/bash', '-i', '-c', 'sasinit; xmm_process > %s 2>&1'%log_file])
        procs.append(proc)
    os.chdir('../..')

# wait for the tasks to end
for p in procs: p.wait()
```


```python
os.chdir('%s/%s'%(base_dir, data_dir))
procs = []
for o in obsids:
    os.chdir(o)
    os.system('mkdir -p pn')
    os.chdir('pn')
    if len(glob.glob('*EVL*')) == 0 and len(glob.glob('pn.fits')) == 0:
        log_file = '../../log/%s_process_pn.log'%o
        p = subp.Popen(['/bin/bash', '-i', '-c', 'sasinit; xmm_process pn > %s 2>&1'%log_file])
        procs.append(p)
    os.chdir('../..')

# wait for the tasks to end
for p in procs: p.wait()
```

### Extra Processing


```python
# For obs: 0761220101; most of the data is in exposure 2; 
os.chdir('%s/%s'%(base_dir, data_dir))
procs = []
for o in ['0761220101']:
    os.chdir('%s/pn'%o)
    if len(glob.glob('pn.fits')) != 0: continue
    cmd = 'export SAS_ODF="$PWD/../odf/"; export SAS_CCF="$SAS_ODF/ccf.cif";'
    cmd += 'sasinit; epchain exposure=2'
    log_file = '../../log/%s_process_pn.log'%o
    p = subp.Popen(['/bin/bash', '-i', '-c', '%s > %s 2>&1'%(cmd, log_file)])
    procs.append(p)
    os.chdir('../..')
# wait for the tasks to end
for p in procs: p.wait()
```

## Spectral Extraction
### Standard Filtering & Region
- `xmm_filter.py` does standard background filtering and opens `ds9` and requrest a region file called `ds9.reg`, which contrains the source and background regions in **Physical** coordinate units, so `xmm_spec.py` can understand it.
- Here, we use an annular region for the source, so later we can check for pileup; at this stage we set the inner radius to 0 and outer radius to 50 arcsec. The background is circular with radius of 50 arcsec.


```python
os.chdir('%s/%s'%(base_dir, data_dir))
exists = os.path.exists
for o in obsids:
    print('-- obs %s --'%o)
    os.chdir('%s/pn'%o)
    if not os.path.exists('pn.fits'):
        # combine multiple exposures if needed #
        evts = glob.glob('*PN*EVL*')
        if len(evts) == 1:
            os.system('ln -s %s pn.fits'%evts[0])
        else:
            cmd = 'sasinit; merge set1=%s set2=%s outset=pn.fits'%(evts[0], evts[1])
            subp.call(['/bin/bash', '-i', '-c', cmd])

    
    os.system('mkdir -p spec')
    os.chdir('spec')
    if not exists('pn_filtered.fits') or not exists('ds9.reg'):
        # check if we have a saved region file, or a temporary region file
        # for faster loading
        saved_reg = '../../../log/%s_ds9.reg'%o
        if exists(saved_reg):
            os.system('cp %s ds9.reg'%saved_reg)
            region = '--region'
        else:
            region = '--region'
            tmp_reg = '../../../log/ds9.reg'
            if exists(tmp_reg):
                os.system('cp %s tmp.reg'%tmp_reg)

        subp.call(['/bin/bash', '-i', '-c', 
                'sasinit; xmm_filter.py ../pn.fits pn --std %s'%region])
        if not exists(saved_reg):
            os.system('cp ds9.reg %s'%tmp_reg) 
    os.chdir('../../..')
    
```

    -- obs 0013140101 --
    -- obs 0201830201 --
    -- obs 0201830401 --
    -- obs 0554170101 --
    -- obs 0761220101 --
    -- obs 0013140201 --
    -- obs 0201830301 --
    -- obs 0201830501 --
    -- obs 0554170201 --


- Any obsids that we don't want to include?
    



```python
no_data = []
obsids = np.sort([o for o in obsids if not o in no_data])
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
```

    There are 9 observations
    0013140101, 0013140201, 0201830201, 0201830301, 0201830401, 0201830501, 0554170101, 0554170201, 0761220101


### Check for pileup
We use `epatplot` to check for pileup starting with an annulus source region with inner radius of 0. A pileup is present if the fraction of expected to predicted single or doubles events deviates from 1 by more than 3 sigma. If there is pileup, we increase the inner radius of the annulus to 3 arcsec and repeat the increase in steps of 0.5 arcsec, until the fractions are consistent with 1.


```python
os.chdir('%s/%s'%(base_dir, data_dir))
print('obs_num | obsid | fractions | radius | pileup?')
for iobs, o in enumerate(obsids):
    if os.path.exists('%s/pn/pileup'%o): continue
    os.chdir('%s/pn'%o)
    os.system('mkdir -p pileup')
    os.chdir('pileup')
    
    pileup = True
    radius = 0.0
    reg_lines = open('../spec/ds9.reg').readlines()
    if 'background' in reg_lines[-1]:
        reg_lines = reg_lines[:-2] + [reg_lines[-1], reg_lines[-2]]
    g = re.match('\(.*\)\n', reg_lines[-1])
    dum_rad = reg_lines[-1].split(',')[2]
    reg_lines[-1] = reg_lines[-1].replace(',' + dum_rad + ',', ',%g,')
    reg_text = ''.join(reg_lines)
    

    while pileup:
        with open('ds9.reg', 'w') as fp: fp.write(reg_text%radius)
        subp.call(['/bin/bash', '-i', '-c', 
           'sasinit; xmm_spec.py ../pn.fits ds9.reg --check_pileup > pileup.log 2>&1'])
        line = [l for l in open('pileup.log').readlines() if '+/-' in l][0]
        frac = np.array(np.array(line.split())[[2,4,6,8]], np.double)
        pileup = (frac[0] > 1+3*frac[1]) or (frac[2] < 1-3*frac[3])
        # 1 arcsec = 20 pixels
        text = '{:4d} | {} | {} | {} | {:d}'.format(iobs+1, o, line[10:-1], radius/20., pileup)
        print(text)
        if pileup: radius = np.max([3.0*20, radius+10])
    os.chdir('../spec')
    os.system('cp ds9.reg ds9.reg.0')
    os.system('cp ../pileup/ds9.reg ds9.reg')
    os.chdir('../../../')
            
```

    obs_num | obsid | fractions | radius | pileup?
       1 | 0013140101 |  s: 1.010 +/- 0.011   d: 0.971 +/- 0.015 | 0.0 | 0
       2 | 0013140201 |  s: 1.020 +/- 0.011   d: 0.958 +/- 0.014 | 0.0 | 0
       3 | 0201830201 |  s: 1.000 +/- 0.011   d: 0.993 +/- 0.014 | 0.0 | 0
       4 | 0201830301 |  s: 0.998 +/- 0.011   d: 0.992 +/- 0.015 | 0.0 | 0
       5 | 0201830401 |  s: 0.996 +/- 0.011   d: 0.998 +/- 0.015 | 0.0 | 0
       6 | 0201830501 |  s: 1.000 +/- 0.009   d: 0.994 +/- 0.012 | 0.0 | 0
       7 | 0554170101 |  s: 0.996 +/- 0.004   d: 0.999 +/- 0.005 | 0.0 | 0
       8 | 0554170201 |  s: 1.001 +/- 0.004   d: 0.990 +/- 0.006 | 0.0 | 0
       9 | 0761220101 |  s: 1.002 +/- 0.005   d: 0.991 +/- 0.006 | 0.0 | 0


### Extract the Spectra


```python
os.chdir('%s/%s'%(base_dir, data_dir))
procs = []
for iobs, o in enumerate(obsids):
    os.chdir('%s/pn/spec'%o)
    if len(glob.glob('spec*grp')) == 0:
        time.sleep(1)
        p = subp.Popen(['/bin/bash', '-i', '-c', 
           'sasinit; xmm_spec.py pn_filtered.fits ds9.reg -o spec_%d > spec.log 2>&1'%(iobs+1)])
        procs.append(p)
    os.chdir('../../..')
# wait for the tasks to end
for p in procs: p.wait()
```

## Summary of PN Spectral Data


```python
os.chdir('%s/%s/'%(base_dir, data_dir))
print('{:5} | {:12} | {:10.8} | {:10.3} | {:10.3}'.format(
        'num', 'obsid', 'mjd', 'rate', 'exposure'))
spec_data = []
for iobs,o in enumerate(obsids):
    with pyfits.open('%s/pn/spec/spec_%d.grp'%(o, iobs+1)) as fp:
        exposure = fp[1].header['exposure']
        counts = fp[1].data.field('counts').sum()
        mjd = (atime.Time(fp[0].header['date_end']).mjd + 
               atime.Time(fp[0].header['date_obs']).mjd ) / 2
        spec_data.append([mjd, counts/exposure, exposure/1e3])
        text = '{:5} | {:12} | {:10.8} | {:10.3} | {:10.3}'.format(
                iobs+1, o, mjd, counts/exposure, exposure/1e3)
        print(text)
spec_data = np.array(spec_data)
```

    num   | obsid        | mjd        | rat        | exp       
        1 | 0013140101   |  51943.034 |       7.31 |       13.5
        2 | 0013140201   |  52283.816 |       12.0 |       9.94
        3 | 0201830201   |  53197.533 |       8.34 |       14.8
        4 | 0201830301   |  53201.043 |       8.02 |       14.0
        5 | 0201830401   |  53208.674 |       7.51 |       13.9
        6 | 0201830501   |  53224.964 |       11.9 |       14.0
        7 | 0554170101   |  54834.329 |       13.4 |       57.7
        8 | 0554170201   |  54674.832 |       12.1 |       62.9
        9 | 0761220101   |  57211.684 |       7.25 |       88.8



```python
## keep only exposures > 5ks
igood = np.argwhere(spec_data[:,2] >= 5)[:,0]
spec_obsids = obsids[igood]
spec_data = spec_data[igood]
print('There are %d spec observations'%len(spec_obsids))
print(', '.join(spec_obsids))
```

    There are 9 spec observations
    0013140101, 0013140201, 0201830201, 0201830301, 0201830401, 0201830501, 0554170101, 0554170201, 0761220101



```python
## save useful data for later #
os.chdir('%s/%s'%(base_dir, data_dir))
# save some useful data for other notebooks
np.savez('log/data.npz', obsids=obsids, spec_obsids=obsids, spec_data=spec_data)
```

---
---
## Process RGS Data


```python
os.chdir('%s/%s'%(base_dir, data_dir))
procs = []
for o in obsids:
    os.chdir(o)
    os.system('mkdir -p rgs')
    os.chdir('rgs')
    if len(glob.glob('spec_rgs*')) == 0:
        log_file = '../../log/%s_process_rgs.log'%o
        p = subp.Popen(['/bin/bash', '-i', '-c', 'sasinit; xmm_process rgs > %s 2>&1'%log_file])
        time.sleep(1)
        procs.append(p)
    os.chdir('../..')
# wait for the tasks to end
for p in procs: p.wait()
```


```python
obsids
```




    array(['0013140101', '0013140201', '0201830201', '0201830301',
           '0201830401', '0201830501', '0554170101', '0554170201',
           '0761220101'], dtype='<U10')



#### Rename the RGS spectra


```python
os.chdir('%s/%s'%(base_dir, data_dir))
for iobs,o in enumerate(obsids):
    os.chdir('%s/rgs'%o)
    os.system('rm *spec_rgs*')
    cmd = ('sasinit;rgscombine pha="spec_r1.src spec_r2.src" bkg="spec_r1.bgd spec_r2.bgd" '
           'rmf="spec_r1.rsp spec_r2.rsp" filepha="spec_rgs.src" filebkg="spec_rgs.bgd" '
           'filermf="spec_rgs.rsp"; grppha spec_rgs.src spec_rgs.grp "group min 20&exit"')
    subp.call(['/bin/bash', '-i', '-c', cmd])
    if len(glob.glob('spec_rgs_%d*'%iobs)) == 0:
        os.system('rename _rgs. _rgs_%d. spec_rgs*'%(iobs+1))
        root = 'spec_rgs_%d.'%(iobs+1) + '%s'
        with pyfits.open(root%'grp') as fp:
            fp[1].header['backfile'] = root%'bgd'
            fp[1].header['respfile'] = root%'rsp'
            os.system('rm tmp.grp > /dev/null 2>&1')
            fp.writeto('tmp.grp')
        os.system('mv %s _%s'%(root%'grp', root%'grp'))
        os.system('mv %s tmp.grp'%(root%'grp'))
        cmd = ('export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'
              'ftgrouppha tmp.grp %s opt respfile=%s')%(root%'grp', root%'rsp')
        subp.call(['/bin/bash', '-i', '-c', cmd])
        # use ogrppha.py
        cmd = ('ogrppha.py %s %s -s 4')%(root%'grp', root%'grp2')
        subp.call(['/bin/bash', '-i', '-c', cmd])
    os.chdir('../..')

```

---
## Light curve Extraction
Note that we use a rate of 0.6 for the filtering, higher than the standard 0.4 used for the spectral analysis


```python
os.chdir('%s/%s'%(base_dir, data_dir))
for iobs, o in enumerate(obsids):
    os.chdir('%s/pn'%o)
    os.system('mkdir -p lc')
    os.chdir('lc')
    os.system('ln -s ../../odf/ccf.cif ccf.cif >/dev/null 2>&1')
    if not os.path.exists('pn_filtered.fits'):
        # use region from spec, but with no central region extracion #
        reg_lines = open('../spec/ds9.reg').readlines()
        g = re.match('\\(.*\\)\\n', reg_lines[-1])
        dum_rad = reg_lines[-1].split(',')[2]
        reg_lines[-1] = reg_lines[-1].replace(',' + dum_rad + ',', ',0,')
        reg_text = ''.join(reg_lines)
        with open('ds9.reg', 'w') as fp: fp.write(reg_text)
        
        p = subp.call(['/bin/bash', '-i', '-c', 
                       'sasinit; xmm_filter.py ../pn.fits pn --std --stdR 0.7'])
    os.chdir('../../..')
```

#### Extraction function


```python
def _extract_lc(lcdir, ebins, dt, extra=''):
    os.chdir('%s/%s'%(base_dir, data_dir))
    procs = []
    for o in obsids:
        os.system('mkdir -p %s/pn/lc/%s'%(o, lcdir))
        os.chdir('%s/pn/lc/%s'%(o, lcdir))
        if len(glob.glob('lc_{:03g}_*.fits'.format(dt))) != 3*(len(ebins.split())-1): 
            cmd = ('sasinit;xmm_lc.py ../pn_filtered.fits ../ds9.reg'
                   ' -e "%s" -t %g %s >lc.log 2>&1')%(ebins, dt, extra)
            time.sleep(0.5)
            p = subp.Popen(['/bin/bash', '-i', '-c', cmd])
            procs.append(p)
        os.chdir('../../../..')
    # wait for the tasks to end
    for p in procs: p.wait()
```


```python
def _extract_lc__p(lcdir, ebins, dt, extra=''):
    """Run in parallel"""
    os.chdir('%s/%s'%(base_dir, data_dir))
    eBins = [' '.join(x) for x in zip(ebins.split()[:-1], ebins.split()[1:])]
    procs = []
    for o in obsids:
        os.system('mkdir -p %s/pn/lc/%s'%(o, lcdir))
        os.chdir('%s/pn/lc/%s'%(o, lcdir))
        if len(glob.glob('lc_{:03g}_*.fits'.format(dt))) != 3*(len(ebins.split())-1):
            for ib, eb in enumerate(eBins):
                wdir = 'tmp_%d'%(ib+1)
                os.system('mkdir -p %s'%wdir); os.chdir(wdir)
                cmd = ('sasinit;xmm_lc.py ../../pn_filtered.fits ../../ds9.reg'
                   ' -e "%s" -t %g %s >lc.log 2>&1')%(eb, dt, extra)
                p = subp.Popen(['/bin/bash', '-i', '-c', cmd])
                procs.append(p)
                os.chdir('..')
                if len(procs) >= 30:
                    for p in procs: p.wait()
                    procs = []
        os.chdir('../../../..')
    # wait for the tasks to end
    for p in procs: p.wait()
    for o in obsids:
        os.chdir('%s/pn/lc/%s'%(o, lcdir))
        for ib, eb in enumerate(eBins):
            wdir = 'tmp_%d'%(ib+1); os.chdir(wdir)
            os.system('rename __1 __%d *'%(ib+1)); os.system('mv lc_* ..')
            os.chdir('..')
            os.system('rm -r %s >/dev/null 2>&1'%wdir)
        os.chdir('../../../..')
```

### 1 bin in 2-10 keV: `1b`
with absolute corrections


```python
lcdir, ebins, dt = '1a', '2 10', 128
_extract_lc(lcdir, ebins, dt, '--abscorr')
```

### 3 bin in 2-10 keV: `4a -> 2 6 7 10`


```python
lcdir, ebins, dt = '3a', '2 6 7 10', 128
_extract_lc__p(lcdir, ebins, dt, '--abscorr')
```

### 4 bin in 2-10 keV: `4a -> 2 4 6 7 10`


```python
lcdir, ebins, dt = '4a', '2 4 6 7 10', 128
_extract_lc(lcdir, ebins, dt, '--abscorr')
```

### 8 bins: `8a -> ' '.join(['{:2.2g}'.format(x) for x in logspace(log10(2), log10(10), 9)])`


```python
lcdir, ebins, dt = '8a', '2 2.4 3 3.7 4.5 5.5 6.7 8.2 10', 128
_extract_lc(lcdir, ebins, dt, '--abscorr')
```

### 8 bins in PHA


```python
chans_from_avg_spec([3,4,5,6,7,8,9], np.linspace(2,10,9), suff='_8ch')
```


```python
lcdir, ebins, dt = '8ch', ('400 600 800 1000 1200 1400 1600 1800 2000'), 128
_extract_lc(lcdir, ebins, dt, '--chans')
```

### 8 bins in PHA (en in log)


```python
chans_from_avg_spec([3,4,5,6,7,8,9], np.logspace(np.log10(2),1,9), suff='_8lch')
```

    files spec_add_8lch.??? created.
    chans: 400 490 599 732 895 1094 1338 1636 2000



```python
lcdir, ebins, dt = '8lch', ('400 490 599 732 895 1094 1338 1636 2000'), 128
_extract_lc__p(lcdir, ebins, dt, '--chans')
```

---
### 12 bins: `12a -> ' '.join(['{:2.2g}'.format(x) for x in logspace(log10(2), log10(10), 13)])`


```python
lcdir, ebins, dt = '12a', '2 2.3 2.6  3 3.4 3.9 4.5 5.1 5.8 6.7 7.6 8.7 10', 128
_extract_lc(lcdir, ebins, dt, '--abscorr')
```

### 12 bins in PHA


```python
chans_from_avg_spec([3,4,5,6,7,8,9], np.linspace(2,10,13), suff='_12ch')
```


```python
lcdir, ebins, dt = '12ch', ('400 534 667 800 934 1067 1200 1334 1467 1600 1734 1867 2000'), 128
_extract_lc(lcdir, ebins, dt, '--chans')
```

### 12 bins in PHA (en in log)


```python
chans_from_avg_spec([3,4,5,6,7,8,9], np.logspace(np.log10(2),1,13), suff='_12lch')
```

    files spec_add_12lch.??? created.
    chans: 400 458 524 599 684 783 895 1023 1170 1338 1530 1749 2000



```python
lcdir, ebins, dt = '12lch', ('400 458 524 599 684 783 895 1023 1170 1338 1530 1749 2000'), 128
_extract_lc__p(lcdir, ebins, dt, '--chans')
```

---
### 16 bins: `16a -> ' '.join(['{:2.2g}'.format(x) for x in logspace(log10(2), log10(10), 17)])`


```python
lcdir, ebins, dt = '16a', '2 2.2 2.4 2.7 3 3.3 3.7 4 4.5 4.9 5.5 6 6.7 7.4 8.2 9 10', 128
_extract_lc(lcdir, ebins, dt, '--abscorr')
```

### 16 bins in PHA


```python
chans_from_avg_spec([3,4,5,6,7,8,9], np.linspace(2,10,17), suff='_16ch')
```


```python
lcdir, ebins, dt = '16ch', ('400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 '
                            '1700 1800 1900 2000'), 128
_extract_lc(lcdir, ebins, dt, '--chans')
```

### 16 bins in PHA (en in log)


```python
chans_from_avg_spec([3,4,5,6,7,8,9], np.logspace(np.log10(2),1,17), suff='_16lch')
```

    files spec_add_16lch.??? created.
    chans: 400 443 490 541 599 662 732 809 895 990 1094 1210 1338 1480 1636 1809 2000



```python
lcdir, ebins, dt = '16lch', ('400 443 490 541 599 662 732 809 895 990 1094 1210 1338 '
                             '1480 1636 1809 2000'), 128
_extract_lc__p(lcdir, ebins, dt, '--chans')
```

---
### 24 bins: `24a -> ' '.join(['{:2.2g}'.format(x) for x in logspace(log10(2), log10(10), 25)])`


```python
lcdir, ebins, dt = '24a', ('2 2.1 2.3 2.4 2.6 2.8  3 3.2 3.4 3.7 3.9 4.2 4.5 4.8 5.1 5.5 '
                           '5.8 6.3 6.7 7.2 7.6 8.2 8.7 9.4 10'), 128
_extract_lc(lcdir, ebins, dt, '--abscorr')
```

### 24 bins in PHA


```python
chans_from_avg_spec([3,4,5,6,7,8,9], np.linspace(2,10,25), suff='_24ch')
```


```python
lcdir, ebins, dt = '24ch', ('400 467 534 600 667 734 800 867 934 1000 1067 1134 1200 1267 '
                            '1334 1400 1467 1534 1600 1667 1734 1800 1867 1934 2000'), 128
_extract_lc(lcdir, ebins, dt, '--chans')
```

### 24 bins in PHA (en in log)


```python
chans_from_avg_spec([3,4,5,6,7,8,9], np.logspace(np.log10(2),1,25), suff='_24lch')
```

    files spec_add_24lch.??? created.
    chans: 400 428 458 490 524 560 599 640 684 732 783 837 895 957 1023 1094 1170 1251 1338 1431 1530 1636 1749 1871 2000



```python
lcdir, ebins, dt = '24lch', ('400 428 458 490 524 560 599 640 684 732 783 837 895 957 '
                             '1023 1094 1170 1251 1338 1431 1530 1636 1749 1871 2000'), 128
_extract_lc__p(lcdir, ebins, dt, '--chans')
```

---
### 32 bins: `32a -> ' '.join(['{:2.2g}'.format(x) for x in logspace(log10(2), log10(10), 33)])`


```python
lcdir, ebins, dt = '32a', ('2 2.1 2.2 2.3 2.4 2.6 2.7 2.8  3 3.1 3.3 3.5 3.7 3.8 4 4.3 '
                           '4.5 4.7 4.9 5.2 5.5 5.8 6 6.4 6.7  7 7.4 7.8 8.2 8.6 9 9.5 10'), 128
_extract_lc(lcdir, ebins, dt, '--abscorr')
```

### 32 bins in PHA


```python
chans_from_avg_spec([3,4,5,6,7,8,9], np.linspace(2,10,33), suff='_32ch')
```

    files spec_add_32ch.??? created.
    chans: 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1850 1900 1950 2000



```python
lcdir, ebins, dt = '32ch', ('400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 '
                            '1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 '
                            '1750 1800 1850 1900 1950 2000'), 128
_extract_lc(lcdir, ebins, dt, '--chans')
```

### 32 bins in PHA (en in log)


```python
chans_from_avg_spec([3,4,5,6,7,8,9], np.logspace(np.log10(2),1,33), suff='_32lch')
```

    files spec_add_32lch.??? created.
    chans: 400 421 443 466 490 515 541 569 599 629 662 696 732 770 809 851 895 941 990 1041 1094 1151 1210 1272 1338 1407 1480 1556 1636 1720 1809 1902 2000



```python
lcdir, ebins, dt = '32lch', ('400 421 443 466 490 515 541 569 599 629 662 696 732 770 '
                             '809 851 895 941 990 1041 1094 1151 1210 1272 1338 1407 '
                             '1480 1556 1636 1720 1809 1902 2000'), 128
_extract_lc__p(lcdir, ebins, dt, '--chans')
```

### 48 bins in PHA (en in log)


```python
chans_from_avg_spec([3,4,5,6,7,8,9], np.logspace(np.log10(2),1,49), suff='_48lch')
```

    files spec_add_48lch.??? created.
    chans: 400 414 428 443 458 473 490 506 524 541 560 579 599 619 640 662 684 708 732 757 783 809 837 865 895 925 957 990 1023 1058 1094 1132 1170 1210 1251 1294 1338 1384 1431 1480 1530 1582 1636 1692 1749 1809 1871 1935 2000



```python
lcdir, ebins, dt = '48lch', ('400 414 428 443 458 473 490 506 524 541 560 579 599 619 640 662 '
                             '684 708 732 757 783 809 837 865 895 925 957 990 1023 1058 1094 '
                             '1132 1170 1210 1251 1294 1338 1384 1431 1480 1530 1582 1636 1692 '
                             '1749 1809 1871 1935 2000'), 128
_extract_lc__p(lcdir, ebins, dt, '--chans')
```


```python

```


```python

```


```python

```

## 12 bins in PHA with equal counts/bin in `spec_9`


```python
os.chdir('%s/%s_spec'%(base_dir, data_dir))
az.misc.spec_2_ebins('spec_9.grp', 12, ebounds=[2,10])
en = np.loadtxt('energy.dat')[:,1]
chans_from_avg_spec([3,4,5,6,7,8,9], en, suff='_12cch')
print(' '.join(['%4.4g'%x for x in en]))
```

    Background File ./spec_9.bgd ...
    Response File ./spec_9.rmf ...
    Results written to energy.dat
    Results written to energy.dat.snr
    Results written to energy.dat.log
    files spec_add_12cch.??? created.
    chans: 400 454 524 592 661 735 816 906 1009 1130 1271 1465 1999
    [2.00377989 2.27287483 2.62256002 2.96251011 3.30749989 3.67750001
     4.08249998 4.53250027 5.04750013 5.65250015 6.35750008 7.32749987
     9.99749947]



```python
lcdir, ebins, dt = '12cch', ('400 454 524 592 661 735 816 906 1009 1130 1271 1465 1999'), 128
_extract_lc(lcdir, ebins, dt, '--chans')
```

## 24 bins in PHA with equal counts/bin in `spec_9`


```python
os.chdir('%s/%s_spec'%(base_dir, data_dir))
az.misc.spec_2_ebins('spec_9.grp', 24, ebounds=[2,10])
en = np.loadtxt('energy.dat')[:,1]
chans_from_avg_spec([3,4,5,6,7,8,9], en, suff='_24cch')
print(' '.join(['%4.4g'%x for x in en]))
```

    Background File ./spec_9.bgd ...
    Response File ./spec_9.rmf ...
    Results written to energy.dat
    Results written to energy.dat.snr
    Results written to energy.dat.log
    files spec_add_24cch.??? created.
    chans: 400 426 454 489 524 558 592 626 661 698 735 774 816 859 906 956 1009 1067 1130 1201 1271 1352 1465 1644 1999
    2.004 2.133 2.273 2.448 2.623 2.793 2.963 3.133 3.307 3.492 3.678 3.872 4.082 4.298 4.533 4.783 5.048 5.338 5.653 6.008 6.358 6.762 7.327 8.222 9.997



```python
lcdir, ebins, dt = '24cch', ('400 426 454 489 524 558 592 626 661 698 735 774 816 859 '
                             '906 956 1009 1067 1130 1201 1271 1352 1465 1644 1999'), 128
_extract_lc(lcdir, ebins, dt, '--chans')
```

## 32 bins in PHA with equal counts/bin in `spec_9`


```python
os.chdir('%s/%s_spec'%(base_dir, data_dir))
az.misc.spec_2_ebins('spec_9.grp', 32, ebounds=[2,10])
en = np.loadtxt('energy.dat')[:,1]
chans_from_avg_spec([3,4,5,6,7,8,9], en, suff='_32cch')
print(' '.join(['%4.4g'%x for x in en]))
```

    Background File ./spec_9.bgd ...
    Response File ./spec_9.rmf ...
    Results written to energy.dat
    Results written to energy.dat.snr
    Results written to energy.dat.log
    files spec_add_32cch.??? created.
    chans: 400 420 440 463 489 515 541 567 592 618 643 670 698 726 755 785 816 848 882 918 956 995 1037 1083 1130 1183 1238 1288 1352 1430 1546 1707 1999
    2.004 2.103 2.203 2.318 2.448 2.578 2.708 2.838 2.963 3.092 3.217 3.352 3.492 3.633 3.778 3.928 4.082 4.242 4.412 4.592 4.783 4.977 5.188 5.418 5.653 5.918 6.193 6.443 6.762 7.153 7.733 8.538 9.997



```python
lcdir, ebins, dt = '32cch', ('400 420 440 463 489 515 541 567 592 618 643 670 698 726 '
                             '755 785 816 848 882 918 956 995 1037 1083 1130 1183 1238 '
                             '1288 1352 1430 1546 1707 1999'), 128
_extract_lc(lcdir, ebins, dt, '--chans')
```

## 48 bins in PHA with equal counts/bin in `spec_9`
Force a minimum of 30 channels per bin


```python
os.chdir('%s/%s_spec'%(base_dir, data_dir))
az.misc.spec_2_ebins('spec_9.grp', 48, ebounds=[2,10])
edata = np.loadtxt('energy.dat')
en, chans = edata[:,1], edata[:,0]
ic, ich = 1, [0]
while ic < len(chans):
    if (chans[ic]-chans[ich[-1]]) >= 30: ich.append(ic)
    ic += 1
en = en[ich]
chans_from_avg_spec([3,4,5,6,7,8,9], en, suff='_48cch')
print(' '.join(['%4.4g'%x for x in en]))
```

    Background File ./spec_9.bgd ...
    Response File ./spec_9.rmf ...
    Results written to energy.dat
    Results written to energy.dat.snr
    Results written to energy.dat.log
    files spec_add_48cch.??? created.
    chans: 400 440 471 507 541 575 609 643 679 716 755 795 837 882 931 982 1037 1067 1098 1130 1165 1201 1238 1271 1307 1352 1401 1465 1546 1644 1781 1999
    2.004 2.203 2.358 2.538 2.708 2.878 3.048 3.217 3.398 3.582 3.778 3.977 4.188 4.412 4.658 4.912 5.188 5.338 5.492 5.653 5.827 6.008 6.193 6.358 6.537 6.762 7.008 7.327 7.733 8.222 8.907 9.997



```python
lcdir, ebins, dt = '48cch', ('400 440 471 507 541 575 609 643 679 716 755 795 837 882 931 '
                             '982 1037 1067 1098 1130 1165 1201 1238 1271 1307 1352 1401 1465 '
                             '1546 1644 1781 1999'), 128
_extract_lc(lcdir, ebins, dt, '--chans')
```


```python

```

### 16 bins: `16a -> ' '.join(['{:2.2g}'.format(x) for x in np.linspace(2,10,17)])`
apply absolute corrections


```python
lcdir, ebins, dt = '16a', ('2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10'), 128
_extract_lc(lcdir, ebins, dt, '--abscorr')
```

### 24 bins: `24a -> ' '.join(['{:2.2g}'.format(x) for x in np.linspace(2,10,25)])`
apply absolute corrections


```python
lcdir, ebins, dt = '24a', ('2 2.1 2.3 2.4 2.6 2.8  3 3.2 3.4 3.7 3.9 4.2 4.5 4.8 5.1 5.5 '
                           '5.8 6.3 6.7 7.2 7.6 8.2 8.7 9.4 10'), 128
_extract_lc(lcdir, ebins, dt, '--abscorr')
```

<br /><br />

---
## NuSTAR DATA


```python
base_dir
```




    '/home/abzoghbi/data/ngc5506/spec_timing_paper'




```python
os.chdir(base_dir)
data_dir = 'data/nustar'
os.system('mkdir -p %s'%data_dir)
obsids = np.array(['60061323002'])
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
```

    There are 1 observations
    60061323002



```python
# get the data #
os.chdir('%s/%s'%(base_dir, data_dir))
ftp = FTP('legacy.gsfc.nasa.gov', 'anonymous', 'anonymous@gmail.com')
ftp.cwd('nustar/data/obs')
failed = []
for o in obsids:
    ftp.cwd('%s/%s'%(o[1:3], o[0]))
    tar_file = '%s.tar'%o
    # download file only if not already downloaded
    if not os.path.exists(tar_file):
        try:
            ftp.retrbinary('RETR %s'%tar_file ,open(tar_file, 'wb').write)
        except:
            print('failed downloading %s'%o)
            os.system('rm %s >/dev/null 2>&1'%tar_file)
            failed.append(o)
    ftp.cwd('../..')
```


```python
for f in failed:
    if f in obsids:
        obsids = np.delete(obsids, np.argwhere(obsids==f)[0,0])
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
nobs = len(obsids)
```

    There are 1 observations
    60061323002


### Process the NuSTAR data
We use our shell script `nustar_process`.


```python
os.chdir('%s/%s'%(base_dir, data_dir))
os.system('mkdir -p log')
procs = []
for o in obsids:
    if not os.path.exists(o):
        os.system('tar -xf %s.tar'%o)
    if not os.path.exists('%s_p'%o):
        # download large files by http
        
        log_file = 'log/%s_process.log'%o
        cmd = ('export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'
               'nustar_process %s > %s 2>&1'%(o, log_file))
        proc = subp.Popen(['/bin/bash', '-i', '-c', cmd])
        procs.append(proc)
        time.sleep(0.2)

# wait for the tasks to end
for p in procs: p.wait()
```

### Spectral Extraction
- Use `nustar_spec.py` 
- Region size: 150''


```python
os.chdir('%s/%s'%(base_dir, data_dir))
exists = os.path.exists
obsids = np.sort(obsids)
procs = []
for iobs,o in enumerate(obsids):
    print('-- obs %s --'%o)
    os.chdir('%s_p/'%o)
    os.system('mkdir -p spec')
    os.chdir('spec')
    if len(glob.glob('spec*grp')) != 2:
        # check if we have a saved region file, or a temporary region file
        # for faster loading
        saved_reg = '../../log/%s_src.reg'%o
        if exists(saved_reg):
            os.system('cp %s src.reg'%saved_reg)
            os.system('cp %s bgd.reg'%(saved_reg.replace('_src.', '_bgd.')))
            region = ''
        else:
            region = '--create_region'

        cmd = ('export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'
              'nustar_spec.py -o spec_%d %s'%(iobs+1, region))
        p = subp.Popen(['/bin/bash', '-i', '-c', cmd])
        procs.append(p)
        if not exists(saved_reg):
            os.system('cp src.reg %s'%saved_reg) 
            os.system('cp bgd.reg %s'%(saved_reg.replace('_src.', '_bgd.')))
        time.sleep(0.3)
    os.chdir('../..')
# wait for the tasks to end
for p in procs: p.wait() 
```

    -- obs 60061323002 --



```python
## group the spectra #
os.chdir('%s/%s'%(base_dir, data_dir))
for iobs,o in enumerate(obsids):
    os.chdir('%s_p/spec'%o)
    cmd = ('rm *grp; ogrppha.py spec_{0}_a_sr.pha spec_{0}_a.grp -f 3 -s 6;'
           'ogrppha.py spec_{0}_b_sr.pha spec_{0}_b.grp -f 3 -s 6').format(iobs+1)
    subp.call(['/bin/bash', '-i', '-c', cmd])
    os.chdir('../..')
```


```python
os.chdir('%s/%s/'%(base_dir, data_dir))
print('{:5} | {:12} | {:10.8} | {:10.8} | {:10.3} | {:10.3}'.format(
        'num', 'obsid', 'mjd_s', 'mjd_e', 'rate', 'exposure'))
spec_data = []
for iobs,o in enumerate(obsids):
    with pyfits.open('%s_p/spec/spec_%d_a.grp'%(o, iobs+1)) as fp:
        exposure = fp[1].header['exposure']
        counts = fp[1].data.field('counts').sum()
        tmid = np.array([fp[0].header['tstart'], fp[0].header['tstop']])
        mref = fp[0].header['mjdrefi'] + fp[0].header['mjdreff']
        mjd = tmid / (24*3600) + mref
        spec_data.append([mjd, counts/exposure, exposure/1e3])
        text = '{:5} | {:12} | {:10.8} | {:10.8} | {:10.3} | {:10.5}'.format(
                iobs+1, o, mjd[0], mjd[1], counts/exposure, exposure/1e3)
        print(text)
spec_data = np.array(spec_data)
```

    num   | obsid        | mjd_s      | mjd_e      | rat        | exp       
        1 | 60061323002  |  56748.996 |  56750.244 |       1.88 |     56.585



```python

```
