"""
Reading and writing functions

TODO:  see if klgfbRead is needed or not
"""
import os
from datetime import date
import pickle
import configparser
from fnmatch import fnmatch

import numpy as np

import ChiantiPy.tools.util as util
import ChiantiPy.tools.constants as const
import ChiantiPy.Gui as chgui
from  ChiantiPy.fortranformat import FortranRecordReader

today = date.today()
chRef = ['produced as a part of the \'CHIANTI\' atomic database for astrophysical spectroscopy', ' K. Dere (GMU) - ' + today.strftime('%Y %B %d')]


def abundanceRead(abundancename=None, verbose=False):
    """
    Read abundance file `abundancename` and return the abundance values relative to hydrogen

    Keyword Arguments
    -----------------

    abundancename:  `str`
        the name of the abundance file in the $XUVTOP/abundance directory to be read
        the default is an empty string and then the 'default' abundance values are read

    verbose:  `bool`
        if true, prints out some info

    """
    abundance = np.zeros((50),np.float64)
    xuvtop = os.environ["XUVTOP"]
    abundDir = os.path.join(xuvtop,'abundance')
    abundList = os.listdir(abundDir)
    if abundancename:
        if os.path.isfile(abundancename):
            abundanceFileName = abundancename
            # a specific abundance file name has been specified
        else:
            cnt = abundancename.count('.abund')
            if cnt == 0:
                abundancename += '.abund'
            if abundancename in abundList:
                abundanceFileName = os.path.join(abundDir,abundancename)

            else:
                # the user will select an abundance file
                abundLabel = 'ChiantiPy - Select an abundance file'
                mypick = chgui.gui.selectorDialog(abundList, label=abundLabel)
                abundName = mypick.selectedText[0]
                if abundName is None:
                    print((' no abundance file selected'))
                    return 0
                else:
                    abundanceFileName = os.path.join(abundDir,abundName)
    else:
        defaults = defaultsRead()
        abundanceFileName = os.path.join(abundDir, defaults['abundfile'] +'.abund')

    with open(abundanceFileName,'r') as inpt:
        s1 = inpt.readlines()
    nlines = 0
    idx = -1
    while idx <= 0:
        minChar = min([5, len(s1[nlines])])
        aline = s1[nlines][0:minChar]
        idx = aline.find('-1')
        nlines += 1
    nlines -= 1
    for line in range(nlines):
        z,ab,element = s1[line].split()[0:3]
        abundance[int(z)-1] = float(ab)
    gz = np.nonzero(abundance)
    abs = 10.**(abundance[gz]-abundance[0])
    abundance.put(gz,abs)
    abundanceRef = s1[nlines+1:]
    return {'abundancename':abundanceFileName,'abundance':abundance,'abundanceRef':abundanceRef}


def zion2name(z,ion, dielectronic=False):
    """
    Convert `Z` and `ion` to generic name, e.g. 26, 13 -> fe_13

    Parameters
    ----------
    z : `int`
        the nuclear charge, for example 26 for Fe XIV
    ion : `int`
        the ion stage, for example, 14 for Fe XIV

    Keyword Arguments
    -----------------

    dielectronic:  `bool`
        if True, created the name of a dielectronic ion, with a 'd' at the end
    Todo
    ----

    See if dielectronic is still appropriate

    Put in separate module to avoid multiple copies

    Notes
    -----
    A duplicate of the routine in `ChiantiPy.tools.util` but needed by masterList Info

    """
    if ion == 0:
        thisone = 0
    elif ion == z+2:
        thisone = 0
    elif (z-1 < len(const.El)) and (ion <= z+1):
        thisone = const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone += 'd'
    else:
        # this should not actually happen
        thisone = 0
    return thisone


def convertName(name):
    """
    Convert ion name string to Z and Ion and other interesting info

    Parameters
    ----------
    name :  `str`
        a generic name of an ion in the CHIANTI database,
        such as fe_14 for Fe XIV

    Todo
    ----

    Put in separate module to avoid multiple copies

    Notes
    -----
    A duplicate of the routine in `ChiantiPy.tools.util` but needed by masterList Info

    """
    s2 = name.split('_')
    els = s2[0].strip()
    i1 = const.El.index(els)+1
    ions = s2[1].strip()
    d = ions.find('d')
    if d >0 :
        dielectronic = True
        ions = ions.replace('d','')
    else: dielectronic = False
    higher = zion2name(int(i1), int(ions)+1)
    lower = zion2name(int(i1), int(ions)-1)
    return {'Z':int(i1),'Ion':int(ions),'Dielectronic':dielectronic, 'Element':els, 'higher':higher, 'lower':lower}


def autoRead(ions, filename=None, total=True, verbose=False):
    """
    Read CHIANTI autoionization rates from a .auto file.

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V
    filename : `str`
        Custom filename, will override that specified by `ions`
    elvlcname : `str`
        If specified, the lsj term labels are returned in the 'pretty1' and 'pretty2'
        keys of 'Wgfa' dict
    total : `bool`
        Return the summed level 2 autoionization rates in 'Auto'
    verbose : `bool`

    Returns
    -------
    Auto : `dict`
        Information read from the .wgfa file. The dictionary structure is
        {"lvl1", "lvl2", "avalue", "pretty1", "pretty2", "ref","ionS", "filename"}

    """
    #
    if filename:
        autoname = filename

    else:
        fname = util.ion2filename(ions)
        autoname = fname+'.auto'
    #
    with open(autoname,'r') as input:
        s1 = input.readlines()

    nwvl = 0
    ndata = 2
    while ndata > 1:
        s1a = s1[nwvl]
        s2 = s1a.split()
        ndata = len(s2)
        nwvl += 1
    nwvl -= 1
    if verbose:
        print((' nwvl = %10i ndata = %4i'%(nwvl, ndata)))
    lvl1 = [0]*nwvl
    lvl2 = [0]*nwvl
    avalue = [0.]*nwvl
    pretty1  =  ['']*nwvl
    pretty2  =  ['']*nwvl
    #
    if verbose:
        print((' nwvl  =  %10i'%(nwvl)))
    #
    wgfaFormat = '(2i7,e12.2,a30,3x,a30)'
    header_line = FortranRecordReader(wgfaFormat)
    for ivl in range(nwvl):
#        inpt=FortranLine(s1[ivl],wgfaFormat)
        inpt = header_line.read(s1[ivl])
        lvl1[ivl] = inpt[0]
        lvl2[ivl] = inpt[1]
        avalue[ivl] = inpt[2]
        pretty1[ivl] = inpt[3].strip()
        pretty2[ivl] = inpt[4].strip()

    ref = []
    for i in range(nwvl+1,len(s1)):
        s1a = s1[i]
        ref.append(s1a.strip())
    Auto = {"lvl1":lvl1, "lvl2":lvl2, "avalue":avalue, "ref":ref, 'ionS':ions, 'filename':autoname, 'pretty1':pretty1, 'pretty2':pretty2}
    if total:
        avalueLvl = [0.]*max(lvl2)
        for lvl,  av in enumerate(avalue):
            l2 = lvl2[lvl] - 1
            avalueLvl[l2] += av
        Auto['avalueLvl'] = np.asarray(avalueLvl)

    if verbose:
        pstring1 = '%5s %5s %12s %12s %20s - %20s'
        print(pstring1%('lvl1', 'lvl2', 'auto value', 'branch ratio', 'pretty1', 'pretty2'))
        pstring = '%5i %5i %12.2e %12.2e %20s - %20s'
        for ivl, avalue in enumerate(avalue):
            l1 = lvl1[ivl]
            l2 = lvl2[ivl]
            br = avalue/avalueLvl[l2-1]
            print(pstring%(l1, l2, avalue, br, pretty1[ivl], pretty2[ivl]))
    #
    return Auto


def autoWrite(info, filename = None, minBranch = None):
    """
    Write data to a CHIANTI .wgfa style file

    Parameters
    ----------
    info : `dict`
        Should contain the following:
        ionS, the Chianti style name of the ion such as c_4 for C IV,
        lvl1, the lower level, the ground level is 1,
        lvl2, the upper level, wvl, the wavelength (in Angstroms),
        avalue, the autoionization rate,
        pretty1, descriptive text of the lower level (optional),
        pretty2, descriptive text of the upper level (optiona),
        ref, reference text, a list of strings
    outfile : `str`
    minBranch : `~numpy.float64`
        The transition must have a branching ratio greater than the specified to be written to the file
    """
    #
#    gname = info['ionS']
    if filename:
        autoname = filename
    else:
        print(' output filename not specified, no file will be created')
        return
    print((' auto file name = ', autoname))
    if minBranch is None:
        minBranch = 0.
    else:
        info['ref'].append(' minimum branching ratio = %10.2e'%(minBranch))
    with open(autoname, 'w') as out:
        #ntrans = len(info['lvl1'])
        nlvl = max(info['lvl2'])
        totalAvalue = np.zeros(nlvl, np.float64)
        if 'pretty1' in info:
            pformat = '%7i%7i%12.2e%30s - %30s'
        else:
            pformat = '%%7i%7i%12.2e'
        for itrans, avalue in enumerate(info['avalue']):
            # for autoionization transitions, lvl1 can be less than zero???
            if abs(info['lvl1'][itrans]) > 0 and info['lvl2'][itrans] > 0:
                totalAvalue[info['lvl2'][itrans] -1] += avalue

        for itrans, avalue in enumerate(info['avalue']):
            if avalue > 0.:
                branch = avalue/totalAvalue[info['lvl2'][itrans] -1]
            else:
                branch = 0.
            if branch > minBranch and abs(info['lvl1'][itrans]) > 0 and info['lvl2'][itrans] > 0:
                if 'pretty1' in info:
                    lbl2 =  info['pretty2'][itrans]
                    pstring = pformat%(info['lvl1'][itrans], info['lvl2'][itrans], avalue, info['pretty1'][itrans].rjust(30), lbl2.ljust(30))
                    out.write(pstring+'\n')
                else:
                    pstring = pformat%(info['lvl1'][itrans], info['lvl2'][itrans], avalue)
                    out.write(pstring+'\n')
        out.write(' -1 \n')
        cnt = 0
        for one in info['ref']:
            if 'filename' in one:
                cnt += 1
        if cnt == 0:
            out.write('%filename:  ' + filename + '\n')
        for one in info['ref']:
            out.write(one+'\n')

def cireclvlRead(ions, filename=None, filetype='cilvl'):
    """
    Read Chianti cilvl, reclvl, or rrlvl files and return data

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V
    filename : `str`, optional
        Custom filename, will override that specified by `ions`
    filetype : `str`
        {'cilvl', 'reclvl', 'rrlvl'} Type of the file to read
    """
    if filename:
        fname = filename
    else:
        fname = util.ion2filename(ions)

    paramname = fname + '.' + filetype

    #print('paramname %s'%(paramname))
    if os.path.isfile(paramname):
        with open(paramname,'r') as input:
            lines = input.readlines()
            input.close()
    else:
        print('file:  %s  does not exist'%(paramname))
        return

    iline = 0
    idx = -1
    while idx < 0:
        aline = lines[iline][0:5]
        idx = aline.find('-1')
        iline += 1
    ndata = iline - 1
    ntrans = ndata//2
    #
    #
    # need to find the maximum number of temperatures, not all lines are the same
    #
    ntemp = np.zeros(ntrans, 'int32')
    iline = 0
    for jline in range(0, ndata, 2):
        dummy = lines[jline].replace(os.linesep, '').split()
        ntemp[iline] = len(dummy[4:])
        iline += 1
    maxNtemp = ntemp.max()
#   print ' maxNtemp = ', maxNtemp
    temp = np.zeros((ntrans,maxNtemp), np.float64)
    iline = 0
    for jline in range(0, ndata, 2):
        recdat = lines[jline].replace(os.linesep, '').split()
        shortT = np.asarray(recdat[4:], np.float64)
        # the result of the next statement is to continue to replicate t
        t = np.resize(shortT, maxNtemp)
        if filetype == 'rrlvl':
            temp[iline] = t
        else:
            temp[iline] = 10.**t
        iline += 1
    #
    lvl1 = np.zeros(ntrans, 'int64')
    lvl2 = np.zeros(ntrans, 'int64')
    ci = np.zeros((ntrans, maxNtemp), np.float64)
    #
    idat = 0
    for jline in range(1, ndata, 2):
        cidat = lines[jline].replace(os.linesep, '').split()
        shortCi = np.asarray(cidat[4:], np.float64)
        lvl1[idat] = int(cidat[2])
        lvl2[idat] = int(cidat[3])
        ci[idat] = np.resize(shortCi, maxNtemp)
        idat += 1
    return {'temperature':temp, 'ntemp':ntemp,'lvl1':lvl1, 'lvl2':lvl2, 'rate':ci,'ref':lines[ndata+1:], 'ionS':ions}


def defaultsRead(verbose=False):
    """
    Read in configuration from .chiantirc file or set defaults if one is not found.
    """
    initDefaults = {'abundfile': 'sun_photospheric_2015_scott','ioneqfile': 'chianti', 'wavelength': 'angstrom', 'flux': 'energy','gui':True}

    if "HOME" in os.environ.keys():
        if os.path.isdir(os.environ["HOME"]):
            home = os.environ["HOME"]
            fn1 = os.path.join(home, '.config', 'chiantirc')
            fn2 = os.path.join(home, '.chianti', 'chiantirc')
            rcfile = False
            if os.path.isfile(fn1):
                rcfile = fn1
            elif os.path.isfile(fn2):
                rcfile = fn2
    elif "PROFILEHOME" in os.environ.keys() and rcfile is False:
        if os.path.isdir(os.environ["PROFILEHOME"]):
            home = os.environ["PROFILEHOME"]
            fn1 = os.path.join(home, '.chianti', 'chiantirc')
            fn2 = os.path.join(home, '.config', 'chiantirc')
            rcfile = False
            if os.path.isfile(fn1):
                rcfile = fn1
            elif os.path.isfile(fn2):
                rcfile = fn2
    if verbose:
        print('rcfile:  %s'%(rcfile))

    if os.path.isfile(rcfile):
        print((' reading chiantirc file'))
        config = configparser.RawConfigParser(initDefaults)
        config.read(rcfile)
        defaults = {}
        for anitem in config.items('chianti'):
            defaults[anitem[0]] = anitem[1]
        if defaults['gui'].lower() in ('t', 'y', 'yes', 'on', 'true', '1', 1, True):
            defaults['gui'] = True
        elif defaults['gui'].lower() in ('f', 'n', 'no', 'off', 'false', '0', 0, False):
            defaults['gui'] = False
        defaults['rcfile'] = rcfile
    else:
        defaults = initDefaults
        if verbose:
            print((' chiantirc file can not be found'))
            print(' should be in $HOME/.config, $HOME/.chianti or $PROFILEHOME/.chianti')
            print((' using the following defaults'))
            for akey in list(defaults.keys()):
                print((' %s = %s'%(akey, defaults[akey])))
    return defaults


def diRead(ions, filename=None):
    """
    Read chianti direct ionization .params files and return data.

    Parameters
    ----------

    ions : `str`
        Ion, e.g. 'c_5' for C V

    filename : `str`, optional
        Custom filename, will override that specified by `ions`
    """
    #
    if filename:
        paramname = filename
    else:
        zion = util.convertName(ions)
        if zion['Z'] < zion['Ion']:
            print(' this is a bare nucleus that has no ionization rate')
            return {'errorMessage':'diParams do not exist for this ion'}
        #
        fname = util.ion2filename(ions)
        paramname = fname+'.diparams'
    #
    input = open(paramname,'r')
    #  need to read first line and see how many elements
    line1 = input.readline()
    indices = line1.split()
    iz = int(indices[0])
    ion = int(indices[1])
    nspl = indices[2]
    nfac = int(indices[3])
    neaev = int(indices[4])
    nspl = int(nspl)
#    format=FortranFormat(str(nspl+1)+'E10.2')
    header_line = FortranRecordReader(str(nspl+1)+'E10.2')
#    inpt =header_line.read(s1[i][0:115])

    #
    ev1 = np.zeros(nfac,np.float64)
    btf = np.zeros(nfac,np.float64)
    xsplom = np.zeros([nfac, nspl],np.float64)
    ysplom = np.zeros([nfac, nspl],np.float64)
    #
    for ifac in range(nfac):
        line = input.readline()
#        paramdat=FortranLine(line,format)
        paramdat = header_line.read(line)
        btf[ifac] = paramdat[0]
        xsplom[ifac] = paramdat[1:]
        line = input.readline()
#        paramdat=FortranLine(line,format)
        paramdat = header_line.read(line)
        ev1[ifac] = paramdat[0]
        ysplom[ifac] = paramdat[1:]
    if neaev:
        line = input.readline()
        eacoef = line.split()
        eaev = [float(avalue) for avalue in eacoef]
    else:
        eaev = 0.
    # skip the  -1
    line = input.readline()
    ref = []
    hdr = input.readlines()
    for ahdr in hdr:
        ref.append(ahdr.strip())
    input.close()
    info = {"iz":iz,"ion":ion,"nspl":nspl,"neaev":neaev, 'nfac':nfac}
    if neaev:
        info['eaev'] = eaev
    DiParams = {"info":info,"btf":btf,"ev1":ev1,"xsplom":xsplom,"ysplom":ysplom, 'eaev':eaev,"ref":ref}
    return DiParams


def drRead(ions, filename=None):
    """
    Read CHIANTI dielectronic recombination .drparams files
    if filename is set, then reads that file

    Parameters
    ----------

    ions : `str`
        Ion, e.g. 'c_5' for C V

    filename : `str`, optional
        Custom filename, will override that specified by `ions`
    """
    #
    #
    if filename:
        paramname = filename
    else:
        fname = util.ion2filename(ions)
        paramname = fname+'.drparams'
    if os.path.isfile(paramname):
        input = open(paramname,'r')
        #  need to read first line and see how many elements
        lines = input.readlines()
        input.close()
        drtype = int(lines[0])
        ref = lines[4:-1]
        #
        if drtype == 1:
            # a Badnell type
            nparams = len(lines[1].split())
            if nparams == 10:
                header_line =  FortranRecordReader('2i5,8e12.4')
            elif nparams == 11:
                header_line =  FortranRecordReader('2i5,9e12.4')
            inpt1 = header_line.read(lines[1])
            inpt2 = header_line.read(lines[2])
            eparams = np.asarray(inpt1[2:], np.float64)
            cparams = np.asarray(inpt2[2:], np.float64)
            DrParams = {'drtype':drtype, 'eparams':eparams,'cparams':cparams,  'ref':ref}
        elif drtype == 2:
            # shull type
#            fmt=FortranFormat('2i5,4e12.4')
            header_line =  FortranRecordReader('2i5,4e12.4')
            inpt1 = header_line.read(lines[1])
#            params=np.asarray(FortranLine(lines[1],fmt)[2:], np.float64)
            params = np.asarray(inpt1[2:], np.float64)
            DrParams = {'drtype':drtype, 'params':params, 'ref':ref}
        else:
            DrParams = None
            print((' for ion %5s unknown DR type = %5i' %(ions, drtype)))
    else:
        DrParams = None
    return DrParams


def eaRead(ions, filename=None):
    """
    Read a CHIANTI excitation-autoionization file and calculate the EA ionization rate data
    derived from splupsRead.

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V
    filename : `str`, optional
        Custom filename, will override that specified by `ions`
    """
    if filename:
        splupsname = filename
    else:
        zion = util.convertName(ions)
        if zion['Z'] < zion['Ion']:
            print(' this is a bare nucleus that has no ionization rate')
            return
        #
        fname = util.ion2filename(ions)
        splupsname = fname+'.easplups'
    if not os.path.exists(splupsname):
        print((' could not find file:  ', splupsname))
        return {"lvl1":-1}
    # there is splups/psplups data
    else:
        input = open(splupsname,'r')
        s1 = input.readlines()
        input.close()
        nsplups = 0
        ndata = 2
        while ndata > 1:
            s1a = s1[nsplups][:]
            s2 = s1a.split()
            ndata = len(s2)
            nsplups = nsplups+1
        nsplups = nsplups-1
        lvl1 = [0]*nsplups
        lvl2 = [0]*nsplups
        ttype = [0]*nsplups
        gf = [0.]*nsplups
        de = [0.]*nsplups
        cups = [0.]*nsplups
        nspl = [0]*nsplups
        splups = np.zeros((nsplups,9),np.float64)
#        splupsFormat1='(6x,3i3,8e10.0)'
#        splupsFormat2='(6x,3i3,12e10.0)'
        #
        header_line1 = FortranRecordReader('6x,3i3,8e10.0')
        header_line2 = FortranRecordReader('6x,3i3,12e10.0')


        for i in range(0,nsplups):
            try:
                inpt = header_line1.read(s1[i])
#                inpt=FortranLine(s1[i],splupsFormat1)
            except:
                inpt = header_line2.read(s1[i])
#                inpt=FortranLine(s1[i],splupsFormat2)
            lvl1[i] = inpt[0]
            lvl2[i] = inpt[1]
            ttype[i] = inpt[2]
            gf[i] = inpt[3]
            de[i] = inpt[4]
            cups[i] = inpt[5]
            if len(inpt)  > 11:
                nspl[i] = 9
                splups[i].put(list(range(9)),inpt[6:])
            else:
                nspl[i] = 5
                splups[i].put(list(range(5)),inpt[6:])
        #
        de = np.asarray(de,np.float64)
        deEv = const.ryd2Ev*de
        ref = []
        for i in range(nsplups+1,len(s1)):
            s1a = s1[i][:-1]
            ref.append(s1a.strip())
    return {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de, "deryd":de, "deEv":deEv, "cups":cups
                ,"nspl":nspl,"splups":splups,"ref":ref}


def elvlcRead(ions, filename=None, getExtended=False, verbose=False, useTh=True):
    """
    Reads the new format elvlc files.

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V
    filename : `str`, optional
        Custom filename, will override that specified by `ions`
    getExtended : `bool`
    verbose : `bool`
    useTh : `bool`
        If True, the theoretical values (ecmth and erydth) are inserted when
        an energy value for ecm or eryd is zero(=unknown)
    """
    #
    #
    '%7i%30s%5s%5i%5s%5.1f%15.3f%15.3f \n'
    #
    header_line = FortranRecordReader('i7,a30,a5,i5,a5,f5.1,2f15.3')
#    elvlcFormat  = FortranFormat(fstring)
    #
    #
    if filename:
        elvlname = filename
        bname = os.path.basename(filename)
        ions = bname.split('.')[0]
    else:
        fname = util.ion2filename(ions)
        elvlname = fname+'.elvlc'
    if not os.path.isfile(elvlname):
        print((' elvlc file does not exist:  %s'%(elvlname)))
        return {'status':0}
    status = 1
    input = open(elvlname,'r')
    s1 = input.readlines()
    input.close()
    nlvls = 0
    ndata = 2
    while ndata > 1:
        s1a = s1[nlvls][:-1]
        s2 = s1a.split()
        ndata = len(s2)
        nlvls = nlvls+1
    nlvls -= 1
    if verbose:
        print((' nlvls = %i'%(nlvls)))
    lvl = [0]*nlvls
    conf  =  [0]*nlvls
    term = [0]*nlvls
    label = [0]*nlvls
    spin = [0]*nlvls
    spd = [0]*nlvls
    l = ['']*nlvls
    j = [0.]*nlvls
    mult = [0.]*nlvls
    ecm = [0]*nlvls
    ecmth = [0]*nlvls
    pretty = [0]*nlvls
    if getExtended:
        extended = [' ']*nlvls
    for i in range(0,nlvls):
        if verbose:
            print((s1[i][0:115]))
#        inpt = FortranLine(s1[i][0:115],elvlcFormat)
        inpt = header_line.read(s1[i][0:115])
        lvl[i] = inpt[0]
        term[i] = inpt[1].strip()
        label[i] = inpt[2]
        spin[i] = inpt[3]
        spd[i] = inpt[4].strip()
        l[i] = const.Spd.index(spd[i])
        j[i] = inpt[5]
        mult[i] = 2.*inpt[5] + 1.
        ecm[i] = inpt[6]
        ecmth[i] = inpt[7]
        if ecm[i] < 0.:
            if useTh:
                ecm[i] = ecmth[i]
        stuff = term[i].strip() + ' %1i%1s%3.1f'%( spin[i], spd[i], j[i])
        pretty[i] = stuff.strip()
        if getExtended:
            cnt = s1[i].count(',')
            if cnt > 0:
                idx = s1[i].index(',')
                extended[i] = s1[i][idx+1:]
    eryd = [ecm[i]*const.invCm2ryd if ecm[i] >= 0. else -1. for i in range(nlvls)]
    erydth = [ecmth[i]*const.invCm2ryd if ecmth[i] >= 0. else -1. for i in range(nlvls)]
    ref = []
    # this should skip the last '-1' in the file
    for i in range(nlvls+1,len(s1)):
        s1a = s1[i]
        ref.append(s1a.strip())
    newref = [aref.strip() for aref in ref if aref.strip() != '-1']

    info = {"lvl":lvl,"conf":conf, "term":term,'label':label, "spin":spin, "spd":spd, "l":l, "j":j,
             'mult':mult, "ecm":ecm, 'eryd':eryd,'erydth':erydth, "ecmth":ecmth, "ref":newref,
             "pretty":pretty, 'status':status, 'filename':elvlname}
    if getExtended:
        info['extended'] = extended
    return info


def elvlcWrite(info, outfile=None, round=0, addLvl=0, includeRyd=False, includeEv=False):
    '''
    Write Chianti data to .elvlc file.

    Parameters
    ----------
    info : `dict`
        Information about the Chianti data to write. Should contain
        the following keys: ionS, the Chianti style name of the ion such as c_4
        term, a string showing the configuration
        spin, an integer of the spin of the state in LS coupling
        l, an integer of the angular momentum quantum number
        spd, an string for the alphabetic symbol of the angular momemtum, S, P, D, etc
        j, a floating point number, the total angular momentum
        ecm, the observed energy in inverse cm, if unknown, the value is 0.
        eryd, the observed energy in Rydbergs, if unknown, the value is 0.
        ecmth, the calculated energy from the scattering calculation, in inverse cm
        erydth, the calculated energy from the scattering calculation in Rydbergs
        ref, the references in the literature to the data in the input info
    outfile : `str`
        Output filename. ionS + '.elvlc' (in current directory) if None
    round : `int`
        input to 'np.round' to round input values to maintain the correct number of significant figures
    addLvl : `int`
        Add a constant value to the index of all levels
    includeRyd : `bool`
        If True, write the Rydberg energies in the extended area, delimited by a comma
    includeEv : `bool`
        If True, write the energies in eV in the extended area, delimited by a comma

    Notes
    -----
    For use with files created after elvlc format change in November 2012

    See Also
    --------
    ChiantiPy.tools.archival.elvlcWrite : Write .elvlc file using the old format.
    '''
    if outfile:
        elvlcName = outfile
    else:
        try:
            gname = info['ionS']
        except:
            print(' ''ionS'' not included in input dict')
            return
        elvlcName = gname + '.elvlc'
    print((' elvlc file name = ', elvlcName))
    #
#    if not info.has_key('ecmx'):
#        info['ecmx'] = np.zeros_like(info['ecm'])
#    if not info.has_key('erydx'):
#        info['erydx'] = np.zeros_like(info['eryd'])
    if 'label' not in info:
        nlvl = len(info['ecm'])
        info['label'] = [' ']*nlvl
    if 'eryd' not in info:
        info['eryd'] = [x*const.invCm2ryd  if x >= 0. else  -1. for x in info['ecm']]
    if 'erydth' not in info:
        info['erydth'] = [x*const.invCm2ryd  if x >= 0. else  -1. for x in info['ecmth']]
    if 'eV' not in info:
        info['eV'] = [x*const.invCm2Ev if x >= 0. else  -1. for x in info['ecm']]
    if 'eVth' not in info:
        info['eVth'] = [x*const.invCm2Ev if x >= 0. else  -1. for x in info['ecmth']]
   #
    out = open(elvlcName, 'w')
    pformat = '%7i%30s%5s%5i%5s%5.1f%15.3f%15.3f'
    for i,  aterm in enumerate(info['term']):
        thisTerm = aterm.ljust(29)
        thisLabel = info['label'][i].ljust(4)
        pstring = pformat%(i+1+addLvl, thisTerm, thisLabel, info['spin'][i], info['spd'][i],info['j'][i],
            np.round(info['ecm'][i], round), np.round(info['ecmth'][i], round))
        if includeRyd:
             pstring += ' , %15.8f , %15.8f'%(info['eryd'][i], info['erydth'][i])
        if includeEv:
             pstring += ' , %15.8f , %15.8f'%(info['eV'][i], info['eVth'][i])
        pstring += '\n'
        out.write(pstring)
    out.write(' -1\n')
    cnt = 0
    for one in info['ref']:
        if 'filename' in one:
            cnt += 1
    if cnt == 0:
        out.write('%filename:  ' + elvlcName + '\n')

    ref = info['ref']
    newref = [aref.strip() for aref in ref if aref.strip() != '-1']

    for aref in newref:
        out.write(aref + '\n')
    out.close()
    return

def demRead(demName=''):
    """
    Read emission measure file `emName` and return the temperatures, densities and emission measures

    Keyword Arguments
    -----------------

    demName:  `str`
        the name of the differential emission measure file to read in the $XUVTOP/dem directory

    Returns
    -------

    DEM:  `dict`
        keywords are `temperature`, `density`, `dem`, `em`, `dt`, `filename`

    """
    xuvtop = os.environ["XUVTOP"]
    demdir = os.path.join(xuvtop,'dem')
    demList = os.listdir(demdir)
    temp = []
    dens = []
    dem = []
    em = []
    if demName:
        # a specific abundance file name has been specified
        cnt = demName.count('.dem')
        if cnt == 0:
            demName += '.dem'
    if demName in demList:
        demFileName = os.path.join(xuvtop,'dem',demName)
    else:
        # the user will select an emission measure file
        demlabel = 'ChiantiPy - Select a differential emission measure file'
        mypick = chgui.gui.selectorDialog(demList, label=demlabel)
        demName = mypick.selectedText[0]
        if demName is None:
            print((' no DEM file selected'))
            return 0
        else:
            demName = mypick.selectedText[0]
            demFileName = os.path.join(xuvtop,'dem',demName)
    with open(demFileName,'r') as inpt:
        s1 = inpt.readlines()
    nlines = 0
    idx = -1
    while idx <= 0:
        minChar = min([5, len(s1[nlines])])
        aline = s1[nlines][0:minChar]
        idx = aline.find('-1')
        nlines += 1
    nlines -= 1
    for line in range(nlines):
        a, b = s1[line].split()
        temp.append(10.**float(a))
        dem.append(10.**float(b))

    temp = np.asarray(temp)
    dem = np.asarray(dem)


    if 'quiet' in demName:
        dens = 3.e+8*np.ones_like(temp)
    elif 'active' in demName:
        dens = 2.e+9*np.ones_like(temp)
    elif 'flare' in demName:
        dens = 1.e+10*np.ones_like(temp)
    elif 'prom' in demName:
        dens = 3.e+8*np.ones_like(temp)
    elif 'hole' in demName:
        dens = 1.e+8*np.ones_like(temp)
    else:
        dens = 3.e+8*np.ones_like(temp)

    #  dlnT = dt/T
    dt = np.log(10.**.1)*temp
    em = dem*dt

    ref = s1[nlines+1:]
    return{'temperature':temp, 'density':dens, 'dem':dem, 'em':em, 'dt':dt,
        'ref':ref, 'filename':demFileName}

def emRead(filename = None,  directory = None,  verbose=False):
    """
    Read emission measure file `emName` and return the temperatures, densities and emission measures

    Keyword Arguments
    -----------------

    filename:  `str`
        the name of the emission measure file to read in the $XUVTOP/em directory

    Returns
    -------

    EM:  `dict`
        keywords are `temperature`, `density`, `em`, `filename`

    """
    xuvtop = os.environ["XUVTOP"]

#    emList = os.listdir(emdir)
    temp = []
    dens = []
    em = []
    if directory is not None:
        emDir = directory
    else:
        emDir = os.path.join(xuvtop,'em')
    dirList = os.listdir(emDir)
    emList = []
    for one in dirList:
        if fnmatch(one,  '*.em'):
            emList.append(one)
    if verbose:
        for one in emList:
            print(one)
    if filename:
        # a specific abundance file name has been specified
        cnt = filename.count('.em')
        if cnt == 0:
            filename += '.em'
    if filename in emList:
        emFileName = os.path.join(emDir, filename)
    else:
        # the user will select an emission measure file
        emlabel = 'ChiantiPy - Select an emission measure file'
        mypick = chgui.gui.selectorDialog(emList, label=emlabel)
        emName = mypick.selectedText[0]
        if emName is None:
            print((' no EM file selected'))
            return 0
        else:
            emName = mypick.selectedText[0]
            emFileName = os.path.join(emDir,emName)
    with open(emFileName,'r') as inpt:
        s1 = inpt.readlines()
    nlines = 0
    idx = -1
    while idx <= 0:
        minChar = min([5, len(s1[nlines])])
        aline = s1[nlines][0:minChar]
        idx = aline.find('-1')
        nlines += 1
    nlines -= 1
    for line in range(nlines):
        a, b, c = s1[line].split()
        temp.append(float(a))
        dens.append(float(b))
        em.append(float(c))
    ref = []
    for aref in s1[nlines+1:]:
        ref.append(aref.strip())
#    ref = s1[nlines+1:]
    return{'temperature':np.asarray(temp), 'density':np.asarray(dens), 'em':np.asarray(em), 'ref':ref, 'filename':emFileName}


def fblvlRead(ions, filename=None, verbose=False):
    """
    Read a Chianti energy level file for calculating the
    free-bound continuum

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V
    filename : `str`, optional
        Custom filename, will override that specified by `ions`

    """
    fstring = 'i5,a20,2i5,a3,i5,2f20.3'
    header_line = FortranRecordReader(fstring)
    #
    if filename:
        fblvlName = filename
        bname = os.path.basename(filename)
        ions = bname.split('.')[0]
    else:
        fname = util.convertName(ions)['filename']
        fblvlName = fname + '.fblvl'
    if os.path.exists(fblvlName):
        input = open(fblvlName,'r')
        s1 = input.readlines()
        input.close()
        nlvls = 0
        ndata = 2
        while ndata > 1:
            s1a = s1[nlvls][:-1]
            s2 = s1a.split()
            ndata = len(s2)
            nlvls = nlvls+1
        nlvls -= 1
        if verbose:
            print((' nlvls = %5i'%(nlvls)))
        lvl = [0]*nlvls
        conf = [0]*nlvls
        pqn = [0]*nlvls
        l = [0]*nlvls
        spd = [0]*nlvls
        mult = [0]*nlvls
        ecm = [0]*nlvls
        ecmth = [0]*nlvls
        for i in range(0,nlvls):
            if verbose:
                print((s1[i]))
#            inpt=FortranLine(s1[i],elvlcFormat)
            inpt = header_line.read(s1[i])
            lvl[i] = int(inpt[0])
            conf[i] = inpt[1].strip()
            pqn[i] = int(inpt[2])
            l[i] = int(inpt[3])
            spd[i] = inpt[4].strip()
            mult[i] = int(inpt[5])
            if inpt[6] == 0.:
                ecm[i] = float(inpt[7])
            else:
                ecm[i] = float(inpt[6])
                ecmth[i] = float(inpt[7])
        ref = []
        for i in range(nlvls+1,len(s1)-1):
            s1a = s1[i][:-1]
            ref.append(s1a.strip())
        return {"lvl":lvl,"conf":conf,'pqn':pqn,"l":l,"spd":spd,"mult":mult,
            "ecm":ecm,'ecmth':ecmth, 'filename':fblvlName,  'ref':ref}
    else:
        return {'errorMessage':' fblvl file does not exist %s'%(fblvlName)}

def grndLevelsRead():
    ''' to read the grndLevels.dat file
    give the number of ground levels to sum over
    in populate and drPopulate
    '''
    filename = os.path.join(os.environ['XUVTOP'], 'ioneq', 'grndLevels.dat')
    if os.path.isfile(filename):
        with open(filename, 'r') as inpt:
            lines = inpt.readlines()
    else:
        print(' did not fine file %s in %s'%(os.path.split(filename)))
        return

    for i, aline in enumerate(lines):
        if '-1'  in aline:
            divider = i
    grndLevels = []
    for aline in lines[:divider]:
        grndLevels.append(int(aline.split()[1]))
    return grndLevels

def gffRead():
    """
    Read the free-free gaunt factors of [1]_.

    References
    ----------
    .. [1] Sutherland, R. S., 1998, MNRAS, `300, 321
        <http://adsabs.harvard.edu/abs/1998MNRAS.300..321S>`_

    Notes
    -----
    This function reads the file and reverses the values of g2 and u
    """
    xuvtop = os.environ['XUVTOP']
    fileName = os.path.join(xuvtop, 'continuum','gffgu.dat' )
    input = open(fileName)
    lines = input.readlines()
    input.close()
    #
    #  the 1d stuff below is to make it easy to use interp2d
    ngamma = 41
    nu = 81
    nvalues = ngamma*nu
    g2 = np.zeros(ngamma, np.float64)
    g21d = np.zeros(nvalues, np.float64)
    u = np.zeros(nu, np.float64)
    u1d = np.zeros(nvalues, np.float64)
    gff = np.zeros((ngamma, nu), np.float64)
    gff1d = np.zeros(nvalues, np.float64)
    #
    iline = 5
    ivalue = 0
    for ig2 in range(ngamma):
        for iu in range(nu):
            values = lines[iline].split()
            u[iu] = float(values[1])
            u1d[ivalue] = float(values[1])
            g2[ig2] = float(values[0])
            g21d[ivalue] = float(values[0])
            gff[ig2, iu] = float(values[2])
            gff1d[ivalue] = float(values[2])
            iline += 1
            ivalue += 1
    #
    return {'g2':g2, 'g21d':g21d,  'u':u, 'u1d':u1d,  'gff':gff,  'gff1d':gff1d}


def gffintRead():
    """
    Read the integrated free-free gaunt factors of [1]_.

    """
    xuvtop = os.environ['XUVTOP']
    fileName = os.path.join(xuvtop, 'continuum','gffint.dat' )
    input = open(fileName)
    lines = input.readlines()
    input.close()
    #
    ngamma = 41
    g2 = np.zeros(ngamma, np.float64)
    gffint = np.zeros(ngamma, np.float64)
    s1 = np.zeros(ngamma, np.float64)
    s2 = np.zeros(ngamma, np.float64)
    s3 = np.zeros(ngamma, np.float64)
    #
    ivalue = 0
    start = 4
    for iline in range(start,start+ngamma):
        values = lines[iline].split()
        g2[ivalue] = float(values[0])
        gffint[ivalue] = float(values[1])
        s1[ivalue] = float(values[2])
        s2[ivalue] = float(values[3])
        s3[ivalue] = float(values[4])
        ivalue += 1
    #
    return {'g2':g2, 'gffint':gffint, 's1':s1, 's2':s2, 's3':s3}


def itohRead():
    """
    Read in the free-free gaunt factors of [2]_.

    References
    ----------
    .. [2] Itoh, N. et al., 2000, ApJS, `128, 125
        <http://adsabs.harvard.edu/abs/2000ApJS..128..125I>`_
    """
    xuvtop = os.environ['XUVTOP']
    itohName = os.path.join(xuvtop, 'continuum', 'itoh.dat')
    input = open(itohName)
    lines = input.readlines()
    input.close()
    gff = np.zeros((30, 121), np.float64)
    for iline in range(30):
        gff[iline] = np.asarray(lines[iline].split(), np.float64)
    return {'itohCoef':gff}


#def klgfbRead():
#    """
#    Read CHIANTI files containing the free-bound gaunt factors for n=1-6 from [13]_.
#
#    Returns
#    -------
#    {'pe', 'klgfb'} : `dict`
#        Photon energy and the free-bound gaunt factors
#
#    References
#    ----------
#    .. [13] Karzas and Latter, 1961, ApJSS, `6, 167
#        <http://adsabs.harvard.edu/abs/1961ApJS....6..167K>`_
#    """
#    xuvtop = os.environ['XUVTOP']
#    fname = os.path.join(xuvtop, 'continuum', 'klgfb.dat')
#    input = open(fname)
#    lines = input.readlines()
#    input.close()
#    #
#    ngfb = int(lines[0].split()[0])
#    nume = int(lines[0].split()[1])
#
#    gfb = np.zeros((ngfb, ngfb, nume), np.float64)
#    nlines = len(lines)
##        print 'nlines, nume, ngfb = ', nlines,  nume, ngfb
#    pe = np.asarray(lines[1].split(), np.float64)
#    for iline in range(2, nlines):
#        data = lines[iline].split()
#        n = int(data[0])
#        l = int(data[1])
#        gfb[n-1, l] = np.array(data[2:], np.float64)
#    return {'pe':pe, 'klgfb':gfb}

def klgbfnRead(filename):
    """
    Read CHIANTI files containing the free-bound gaunt factors for n=1-6 from [14]_.
    This reads the corrected KL files in CHIANTI version 10+

    Parameters
    ----------
    filename : `str`
        the filename for the KL data in the continuum directory, such as klgfb_1.dat, n=1,6


    Returns
    -------
    [{'pe', 'klgfb'}] : `list`
        Photon energy and the bound-free gaunt factors

    References
    ----------
    .. [14] Karzas and Latter, 1961, ApJSS, `6, 167
        <http://adsabs.harvard.edu/abs/1961ApJS....6..167K>`_
    """
    kl = np.loadtxt(filename)
    klt = kl.transpose()
    pe = klt[0]
    pef = np.flip(pe)
    gf = klt[1:]
#    gft = gf.transpose()
    gff = np.fliplr(gf)
    return {'pe':pef,  'klgbf':gff, 'filename':filename}


def maoParsRead(filename = None):
    ''' to read the mao et al par2.dat file for calculating the
    ratio of rrloss to rrrecomb
    The electron energy-loss rate due to radiative recombination.
    Mao J., Kaastra J., Badnell N.R.
    <Astron. Astrophys. 599, A10 (2017)>
    =2017A&A...599A..10M
    1-  2  I2    ---     s         Isoelectronic sequence number
    4-  5  I2    ---     z         Atomic number
    7- 16  E10.4 ---     a0        Primary fitting parameter
    18- 27  E10.4 ---     b0        Primary fitting parameter
    29- 38  E10.4 ---     c0        Addiational fitting parameter
    40- 49  E10.4 ---     a1        Primary fitting parameter
    51- 60  E10.4 ---     b1        Primary fitting parameter
    62- 71  E10.4 ---     a2        Addiational fitting parameter
    73- 82  E10.4 ---     b2        Addiational fitting parameter
    84- 86  F3.1  ---     mdp       Maximum deviation in percent

    Keyword Arguments
    -----------------

    filename:  `str`

    Returns
    -------

    data: `dict`
    '''

    if filename is None:
        filename = os.path.join(os.environ['XUVTOP'], 'continuum', 'rrloss_mao_2017_pars.dat')
    else:
        if not os.path.isfile(filename):
            print(' filename for maoPars data is incorrect')
            return

    with open(filename, 'r') as inpt:
        lines = inpt.readlines()
    nlines = len(lines)
    # this number is 255
    seq = [0]*nlines
    z = [0]*nlines
    a0 = [0]*nlines
    b0 = [0]*nlines
    c0 = [0]*nlines
    a1 = [0]*nlines
    b1 = [0]*nlines
    a2 = [0]*nlines
    b2 = [0]*nlines
    mdp = [0]*nlines
    for idx, aline in enumerate(lines):
        sline = aline.split()
        seq[idx] = int(sline[0])
        z[idx] = int(sline[1])
        a0[idx] = float(sline[2])
        b0[idx] = float(sline[3])
        c0[idx] = float(sline[4])
        a1[idx] = float(sline[5])
        b1[idx] = float(sline[6])
        a2[idx] = float(sline[7])
        b2[idx] = float(sline[8])
        mdp[idx] = float(sline[9])
    return {'seq':seq, 'z':z, 'a0':a0, 'b0':b0, 'c0':c0, 'a1':a1, 'b1':b1, 'a2':a2, 'b2':b2, 'mdp':mdp}

def ioneqRead(ioneqName='', minIoneq=1.e-20, verbose=False):
    """
    Reads an ioneq file
    ionization equilibrium values less then minIoneq are returns as zeros

    Keyword Arguments
    -----------------

    ioneqName:  `str`
        reads the file in the $XUVTOP/ioneq directory, if a blank, the default is read

    minIoneq:  `float`
        sets values to zero if less the minIoneq

    verbose:  `bool`
        if true, prints into to the terminal
    Returns
    -------
    {'ioneqname','ioneqAll','ioneqTemperature','ioneqRef'} : `dict`
        Ionization equilibrium values and the reference to the literature
    """
    dir = os.environ["XUVTOP"]
    ioneqdir = os.path.join(dir,'ioneq')
    ioneqInfo = util.findFileTypes(ioneqdir, type = '*.ioneq')
    ioneqNames = ioneqInfo['fileName']
    shortNames = []
    for aname in ioneqNames:
        shortNames.append(os.path.splitext(aname)[0])

    ioneqFullNames = ioneqInfo['fullName']
    if ioneqName not in shortNames:
        # the user will select an ioneq file
        gIoneq = chgui.gui.selectorDialog(ioneqNames,label='Select one',  multiChoice=False)
        gIoneq_idx = gIoneq.selectedIndex
        if len(gIoneq_idx) > 0:
            ioneqFileName = ioneqFullNames[gIoneq_idx[0]]
        else:
            print(' no file chosen')
            return
    else:
        index = shortNames.index(ioneqName)
        ioneqFileName = ioneqFullNames[index]
    #
    with open(ioneqFileName,'r') as input:
        s1 = input.readlines()
    ntemp,nele = s1[0].split()
    if verbose:
        print((' ntemp, nele = %5i %5i'%(ntemp, nele)))
    nTemperature = int(ntemp)
    nElement = int(nele)
    #
    header_linet = FortranRecordReader(str(nTemperature)+'f6.2')
    ioneqTemperature = header_linet.read(s1[1])
    ioneqTemperature = np.asarray(ioneqTemperature[:],np.float64)
    ioneqTemperature = 10.**ioneqTemperature
    nlines = 0
    idx = -1
    while idx < 0:
        aline = s1[nlines][0:5]
        idx = aline.find('-1')
        nlines += 1
    nlines -= 1
    #
    #
    header_lineq = FortranRecordReader('2i3,'+str(nTemperature)+'e10.2')
    #
    ioneqAll = np.zeros((nElement,nElement+1,nTemperature),np.float64)
    for iline in range(2,nlines):
        out = header_lineq.read(s1[iline])
        iz = out[0]
        ion = out[1]
        ioneqAll[iz-1,ion-1].put(list(range(nTemperature)),np.asarray(out[2:],np.float64))
    ioneqAll = np.where(ioneqAll > minIoneq, ioneqAll, 0.)
    ioneqRef = []
    for one in s1[nlines+1:]:
        ioneqRef.append(one[:-1])  # gets rid of the \n
    del s1
    return {'ioneqname':ioneqName,'ioneqAll':ioneqAll,'ioneqTemperature':ioneqTemperature,'ioneqRef':ioneqRef}


def ipRead(verbose=False):
    """
    Reads the ionization potential file

    Returns
    -------
    ip : array-like
        Ionization potential (in eV)
    """
    topdir = os.environ["XUVTOP"]
    ipname = os.path.join(topdir, 'ip','chianti.ip')
    ipfile = open(ipname)
    data = ipfile.readlines()
    ipfile.close()
    nip = 0
    ndata = 2
    maxz = 0
    while ndata > 1:
        s1 = data[nip]
        s2 = s1.split()
        ndata = len(s2)
        nip = nip+1
        if int(s2[0]) > maxz:
            maxz = int(s2[0])
    if verbose:
        print((' maxz = %5i'%(maxz)))
    nip = nip-1
    ip = np.zeros((maxz, maxz + 1), np.float64)
    for aline in data[0:nip]:
        s2 = aline.split()
        iz = int(s2[0])
        ion = int(s2[1])
        ip[iz-1, ion-1] = float(s2[2])
    return ip*const.invCm2Ev


def masterListRead():
    """
    Read a CHIANTI masterlist file.

    Returns
    -------
    masterlist : `list`
        All ions in Chianti database
    """
    dir = os.environ["XUVTOP"]
    fname = os.path.join(dir,'masterlist','masterlist.ions')
    input = open(fname,'r')
    s1 = input.readlines()
    input.close()
    masterlist = []
    for i in range(0,len(s1)):
        s1a = s1[i][:-1]
        s2 = s1a.split(';')
        masterlist.append(s2[0].strip())
    return masterlist


def masterListInfo(force=False, verbose=False):
    """
    Get information about ions in the CHIANTI masterlist.

    Keyword Arguments
    -----------------

    force:  `bool`
        if true, recreates the masterListInfo file

    verbose:  `bool`
        if true, prints into to the terminal

    Returns
    -------
    masterListInfo : `dict`
        {'wmin', 'wmax', 'tmin', 'tmax'} Minimum and maximum wavelengths in
        the wgfa file. Minimum and maximum temperatures for which the
        ionization balance is nonzero.

    Notes
    -----
    This function speeds up multi-ion spectral calculations.
    The information is stored in a pickled file 'masterlist_ions.pkl'
    If the file is not found, one will be created.
    """
    dir = os.environ["XUVTOP"]
    infoPath = os.path.join(dir, 'masterlist')
    infoName = os.path.join(dir,'masterlist','masterlist_ions.pkl')
    #masterName=os.path.join(dir,'masterlist','masterlist.ions')
    #
    makeNew = force == True or not os.path.isfile(infoName)
#    if os.path.isfile(infoName):
    if not makeNew:
#       print ' file exists - ',  infoName
        pfile = open(infoName, 'rb')
        masterListInfo = pickle.load(pfile)
        pfile.close
    elif os.access(infoPath, os.W_OK):
        # the file does not exist but we have write access and will create it
        defaults = defaultsRead()
        print((' defaults = %s'%(str(defaults))))
        ioneqName = defaults['ioneqfile']
        ioneq = ioneqRead(ioneqName = ioneqName)
        masterList = masterListRead()
        masterListInfo = {}
        haveZ = [0]*31
        haveStage = np.zeros((31, 31), np.int64)
        haveDielectronic = np.zeros((31, 31), np.int64)
        for one in masterList:
            if verbose:
                print((' ion = %s'%(one)))
            ionInfo = convertName(one)
            z = ionInfo['Z']
            stage = ionInfo['Ion']
            haveZ[z] = 1
            dielectronic = ionInfo['Dielectronic']
            if dielectronic:
                haveDielectronic[z, stage] = 1
            else:
                haveStage[z, stage] = 1
            thisIoneq = ioneq['ioneqAll'][z- 1, stage - 1 + dielectronic]
            good = thisIoneq > 0.
            goodTemp = ioneq['ioneqTemperature'][good]
            tmin = float(goodTemp.min())
            tmax = float(goodTemp.max())
            vgood = thisIoneq == thisIoneq.max()
            vgoodTemp = float(ioneq['ioneqTemperature'][vgood][0])
            wgfa = wgfaRead(one)
            nZeros = wgfa['wvl'].count(0.)
            # two-photon transitions are denoted by a wavelength of zero (0.)
            while nZeros > 0:
                wgfa['wvl'].remove(0.)
                nZeros = wgfa['wvl'].count(0.)
            # unobserved lines are denoted with a negative wavelength
            wvl = np.abs(np.asarray(wgfa['wvl'], np.float64))
            wmin = float(wvl.min())
            wmax = float(wvl.max())
            masterListInfo[one] = {'wmin':wmin, 'wmax':wmax, 'tmin':tmin, 'tmax':tmax, 'tIoneqMax':vgoodTemp}
        masterListInfo['haveZ'] = haveZ
        masterListInfo['haveStage'] = haveStage
        masterListInfo['haveDielectronic'] = haveDielectronic
        #  now do the bare ions from H thru Zn
        #  these are only involved in the continuum
        for iz in range(1, 31):
            ions = zion2name(iz, iz+1)
            thisIoneq = ioneq['ioneqAll'][iz-1, iz]
            good = thisIoneq > 0.
            goodTemp = ioneq['ioneqTemperature'][good]
            tmin = float(goodTemp.min())
            tmax = float(goodTemp.max())
            wmin = 0.
            wmax = 1.e+30
            masterListInfo[ions] = {'wmin':wmin, 'wmax':wmax, 'tmin':tmin, 'tmax':tmax}
        pfile = open(infoName, 'wb')
        pickle.dump(masterListInfo, pfile)
        pfile.close
    else:
        # the file does not exist and we do NOT have write access to creat it
        # will just make an inefficient, useless version
        masterListInfo = {}
        for one in masterList:
            ionInfo = convertName(one)
            z = ionInfo['Z']
            stage = ionInfo['Ion']
            dielectronic = ionInfo['Dielectronic']
            wmin = 0.
            wmax = 1.e+30
            masterListInfo[one] = {'wmin':wmin, 'wmax':wmax, 'tmin':1.e+4, 'tmax':1.e+9}
        #  now do the bare ions from H thru Zn
        #  these are only involved in the continuum
        for iz in range(1, 31):
            ions = zion2name(iz, iz+1)
            wmin = 0.
            wmax = 1.e+30
            masterListInfo[ions] = {'wmin':wmin, 'wmax':wmax, 'tmin':1.e+4, 'tmax':1.e+9}
        pfile = open(infoName, 'wb')
        pickle.dump(masterListInfo, pfile)
        pfile.close
        masterListInfo = {'noInfo':'none'}
    return masterListInfo


def rrRead(ions, filename=None):
    """
    Read CHIANTI radiative recombination .rrparams files


    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V

    Keyword Arguments
    -----------------

    filename : `str`, optional
        Custom filename, will override that specified by `ions`

    Returns
    -------
    {'rrtype','params','ref'} : `dict`
    """
    #
    #
    if filename:
        paramname = filename
    else:
        fname = util.ion2filename(ions)
        paramname = fname+'.rrparams'
    if os.path.isfile(paramname):
        input = open(paramname,'r')
        #  need to read first line and see how many elements
        lines = input.readlines()
        input.close()
        rrtype = int(lines[0])
        ref = lines[3:-1]
        #
        if rrtype == 1:
            # a Badnell type
#            fmt=FortranFormat('3i5,e12.4,f10.5,2e12.4')
            header_line =  FortranRecordReader('3i5,e12.4,f10.5,2e12.4')
#            params=FortranLine(lines[1],fmt)
            params = header_line.read(lines[1])
            RrParams = {'rrtype':rrtype, 'params':params, 'ref':ref}
        elif rrtype == 2:
            # a Badnell type
#            fmt=FortranFormat('3i5,e12.4,f10.5,2e11.4,f10.5,e12.4')
            header_line =  FortranRecordReader('3i5,e12.4,f10.5,2e11.4,f10.5,e12.4')
#            params=FortranLine(lines[1],fmt)
            params = header_line.read(lines[1])
            RrParams = {'rrtype':rrtype, 'params':params, 'ref':ref}
        elif rrtype == 3:
            # a Shull type
#            fmt=FortranFormat('2i5,2e12.4')
            header_line =  FortranRecordReader('2i5,2e12.4')
#            params=FortranLine(lines[1],fmt)
            params  =  header_line.read(lines[1])
            RrParams={'rrtype':rrtype, 'params':params, 'ref':ref}
        else:
            RrParams = None
            print((' for ion %5s unknown RR type = %5i' %(ions, rrtype)))
        return RrParams
    else:
        return {'rrtype':-1}

def rrLossRead():
    ''' to read the Mao 2017 rr loss parameters [12]_

    References
    ----------
    .. [12] Mao J., Kaastra J., Badnell N.R., `2017 Astron. Astrophys. 599, A10
        <http://adsabs.harvard.edu/abs/2017A%26A...599A..10M>`_
    '''
    filename = os.path.join(os.environ['XUVTOP'], 'continuum', 'rrloss_mao_2017_pars.dat')
    inpt = open(filename, 'r')
    lines = inpt.readlines()
    inpt.close()
    iso = []
    z = []
    a0 = []
    b0 = []
    c0 = []
    a1 = []
    b1 = []
    a2 = []
    b2 = []
    mdp = []
    for aline in lines:
        iso.append(int(aline.split()[0]))
        z.append(int(aline.split()[1]))
        a0.append(float(aline.split()[2]))
        b0.append(float(aline.split()[3]))
        c0.append(float(aline.split()[4]))
        a1.append(float(aline.split()[5]))
        b1.append(float(aline.split()[6]))
        a2.append(float(aline.split()[7]))
        b2.append(float(aline.split()[8]))
        mdp.append(float(aline.split()[9]))

    return {'iso':iso, 'z':z, 'a0':a0, 'b0':b0, 'c0':c0, 'a1':a1, 'b1':b1, 'a2':a2, 'b2':b2, 'mdp':mdp}


def scupsRead(ions, filename=None, verbose=False):
    '''
    Read the new format v8 scups file containing the scaled temperature and upsilons from [8]_.

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V

    Keyword Arguments
    -----------------

    filename : `str`, optional
        Custom filename, will override that specified by `ions`
    verbose : `bool`
        if True, prints info to the terminal

    '''
    #
    if filename:
        scupsFileName = filename
        bname = os.path.basename(scupsFileName)
        ions = bname.split('.')[0]
    else:
        fname = util.ion2filename(ions)
        scupsFileName = fname+'.scups'
    if not os.path.isfile(scupsFileName):
        print((' elvlc file does not exist:  %s'%(scupsFileName)))
        return {'status':0}
    #status = 1
    #
    if os.path.isfile(scupsFileName):
        if verbose:
            print(' scupsFileName = %s'%(scupsFileName))
        inpt = open(scupsFileName)
        lines = inpt.readlines()
        inpt.close()
    else:
        print(('file does not exist: '+str(scupsFileName)))
        return {'errorMessage':'file does not exist' +str(scupsFileName)}
        return
    #ll = lines[1].split()
    #temp = np.asarray(ll[3:], np.float64)
    minusOne = 0
    counter = 0
    while not minusOne:
        if '-1' in lines[counter][:4]:
            minusOne = 1
        else:
            counter += 1
    ntrans = (counter)/3
    #print(' counter %10i ntrans %10i'%(counter, ntrans))
    lvl1 = []
    lvl2 = []
    de = []
    gf = []
    lim = []
    ttype = []
    cups = []
    ntemp = []
    btemp = []
    bscups = []
    counter = 0
    #print(' counter %10i ntrans %10i'%(counter, ntrans))
    # the int seems to be needed for Python3
    for itrans in range(int(ntrans)):
        if verbose:
            print((lines[counter]))
            print((lines[counter+1]))
            print((lines[counter+2]))
        ll1 = lines[counter].split()
        lvl1.append(int(ll1[0]))
        lvl2.append(int(ll1[1]))
        de.append(float(ll1[2]))
        gf.append(float(ll1[3]))
        lim.append(float(ll1[4]))
        ntemp.append(int(ll1[5]))
        ttype.append(int(ll1[6]))
        cups.append(float(ll1[7]))
        ll2 = lines[counter+1].split()
        ll3 = lines[counter+2].split()
#        print ' ll2 = ', ll2
#        print ' gf = ', ll2[2]
        btemp.append(np.asarray(ll2, np.float64))
        bscups.append(np.asarray(ll3, np.float64))
        counter += 3
    counter += 1
    ref = []
    for aline in lines[counter:-1]:
        ref.append(aline.strip('\n'))
    return {'ions':ions, 'lvl1':lvl1, 'lvl2':lvl2, 'de':de, 'gf':gf, 'lim':lim, 'ttype':ttype,'cups':cups,'ntemp':ntemp, 'btemp':btemp, 'bscups':bscups, 'ntrans':ntrans, 'ref':ref}


def splomRead(ions, ea=False, filename=None):
    """
    Read chianti .splom files

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V

    Keyword Arguments
    -----------------

    ea : `bool`
        if true, reads the  .easplom file
    filename : `str`, optional
        Custom filename, will override that specified by `ions`

    Returns
    -------
    {'lvl1', 'lvl2', 'ttype', 'gf', 'deryd', 'c', 'splom', 'ref'} : `dict`

    Notes
    -----
    Still needed for ionization cross sections
    """
    #
    if filename is None:
        fname = util.ion2filename(ions)
        if ea:
            splomname = fname+'.easplom'
        else:
            splomname = fname+'.splom'
    else:
        splomname = filename
    with open(splomname,'r') as inpt:
        lines = inpt.readlines()
    #  need to read first line and see how many elements
    #indices=line1[0:15]
    remainder = lines[0][16:]
    nom = remainder.split(' ')
#    format = FortranFormat('5i3,'+str(len(nom))+'E10.2')
    header_line = FortranRecordReader('5i3,'+str(len(nom))+'E10.2')
    #  go back to the beginning
    #  find the -1
    marker = []
    for iline,  aline in enumerate(lines):
        if len(aline.split()) == 1:
            marker.append(iline)
        else:
            pass
    iline = 0
    z = []
    stage = []
    lvl1 = []
    lvl2 = []
    ttype = []
    gf = []
    de = []
    f = []
    splom = []
    for iline in range(marker[0]):
#        splomdat = FortranLine(lines[iline],format)
        splomdat = header_line.read(lines[iline])
        z.append(int(splomdat[0]))
        stage.append(int(splomdat[1]))
        l1 = splomdat[2]
        l2 = splomdat[3]
        tt1 = splomdat[4]
        gf1 = splomdat[5]
        de1 = splomdat[6]
        f1 = splomdat[7]
        splom1 = splomdat[8:]
        lvl1.append(int(l1))
        lvl2.append(int(l2))
        ttype.append(int(tt1))
        gf.append(float(gf1))
        de.append(float(de1))
        f.append(float(f1))
        splom.append(splom1)
        iline = iline+1
    hdr = lines[iline+1:-1]
    ref = []
    for ahdr in hdr:
        ref.append(ahdr.strip())
    de = np.asarray(de,np.float64)
    deEv = de*const.ryd2Ev
    splomout = np.asarray(splom,np.float64)
    splomout = np.transpose(splomout)
    # note:  de is in Rydbergs
    splom = {"z":z, "stage":stage, "lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de, "deryd":de, "deEv":deEv,
        "c":f,"splom":splomout,"ref":ref}
    return  splom


def splupsRead(ions, filename=None, filetype='splups'):
    """
    Read a CHIANTI .splups file

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V

    Keyword Arguments
    -----------------
    filename : `str`, optional
        Custom filename, will override that specified by `ions`
    filetype : `str`, optional
        {`psplups`,`cisplups`,`splups`} Type of file to read

    Returns
    -------
    {'lvl1', 'lvl2', 'ttype', 'gf', 'de', 'cups', 'bsplups', 'ref'} : `dict`
    """
    #
    if filename:
        splupsname = filename
    else:
        fname = util.ion2filename(ions)
        splupsname = fname+'.'+filetype
    if not os.path.exists(splupsname):
        #TODO: raise exception here or just let the open() function do that for us
        return {'file not found':splupsname}
    # there is splups/psplups data
    else:
        input = open(splupsname,'r')
        s1 = input.readlines()
        input.close()
        nsplups = 0
        ndata = 2
        while ndata > 1:
            s1a = s1[nsplups][:]
            s2 = s1a.split()
            ndata = len(s2)
            nsplups = nsplups+1
        nsplups = nsplups-1
        lvl1 = [0]*nsplups
        lvl2 = [0]*nsplups
        ttype = [0]*nsplups
        gf = np.zeros(nsplups, np.float64)
        de = np.zeros(nsplups, np.float64)
        cups = np.zeros(nsplups, np.float64)
        nspl = [0]*nsplups
#        splups=np.zeros((nsplups,9),np.float64)
        splups = [0.]*nsplups
        if filetype == 'psplups':
#            splupsFormat1 = FortranFormat('3i3,8e10.3')
#            splupsFormat2 = FortranFormat('3i3,3e10.3')
            header_line = FortranRecordReader('3i3,3e10.3')
        else:
#            splupsFormat1='(6x,3i3,8e10.3)'
#            splupsFormat2 = FortranFormat('6x,3i3,3e10.3')
            header_line = FortranRecordReader('6x,3i3,3e10.3')        #
        for i in range(0,nsplups):
#            inpt=FortranLine(s1[i],splupsFormat2)
            inpt = header_line.read(s1[i])
            lvl1[i] = inpt[0]
            lvl2[i] = inpt[1]
            ttype[i] = inpt[2]
            gf[i] = inpt[3]
            de[i] = inpt[4]
            cups[i] = inpt[5]
            if filetype == 'psplups':
                as1 = s1[i][39:].rstrip()
            else:
                as1 = s1[i][45:].rstrip()
            nspl[i] = len(as1)//10
#            splupsFormat3 = FortranFormat(str(nspl[i])+'E10.2')
#            splupsFormat3 = '(' + str(nspl[i]) + 'e10.3' + ')'
            header_line3 = FortranRecordReader('%i'%(nspl[i])+'e10.3' )
#            inpt = FortranLine(as1, splupsFormat3)
            inpt = header_line3.read(as1)
            spl1 = np.asarray(inpt[:], np.float64)
            splups[i] = spl1
        #
        ref = []
        for i in range(nsplups+1,len(s1)):
            s1a = s1[i][:-1]
            ref.append(s1a.strip())
        if filetype == 'psplups':
#            self.Npsplups=nsplups
#            self.Psplups={"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
#                ,"nspl":nspl,"splups":splups,"ref":ref}
            return {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
                ,"nspl":nspl,"splups":splups,"ref":ref, 'filename':splupsname}
        else:
#            self.Splups={"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
#                ,"nspl":nspl,"splups":splups,"ref":ref}
            return {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
                ,"nspl":nspl,"splups":splups,"ref":ref, 'filename':splupsname}


def twophotonHRead():
    """
    Read the two-photon Einstein A values and distribution function for the H sequence.

    Returns
    -------
    {'y0', 'z0', 'avalue', 'asum', 'psi0'} : `dict`
    """
    xuvtop = os.environ['XUVTOP']
    fName = os.path.join(xuvtop, 'continuum', 'hseq_2photon.dat')
    dFile = open(fName, 'r')
    a = dFile.readline()
    y0 = np.asarray(a.split())
    a = dFile.readline()
    z0 = np.asarray(a.split())
    nz = 30
    avalue = np.zeros(nz, np.float64)
    asum = np.zeros(nz, np.float64)
    psi0 = np.zeros((nz, 17), np.float64)
    for iz in range(nz):
        a = dFile.readline().split()
        avalue[iz] = float(a[1])
        asum[iz] = float(a[2])
        psi = np.asarray(a[3:])
        psi0[iz] = psi
    dFile.close()
    return {'y0':y0, 'z0':z0, 'avalue':avalue, 'asum':asum, 'psi0':psi0.reshape(30, 17)}


def twophotonHeRead():
    """
    Read the two-photon Einstein A values and distribution function for the He sequence.

    Returns
    -------
    {'y0', 'avalue', 'psi0'} : `dict`
    """
    xuvtop = os.environ['XUVTOP']
    fName = os.path.join(xuvtop, 'continuum', 'heseq_2photon.dat')
    dFile = open(fName, 'r')
    a = dFile.readline()
    y0 = np.asarray(a.split())
    nz = 30
    avalue = np.zeros(nz, np.float64)
    psi0 = np.zeros((nz, 41), np.float64)
    for iz in range(1, nz):
        a = dFile.readline().split()
        avalue[iz] = float(a[1])
        psi = np.asarray(a[2:])
        psi0[iz] = psi
    dFile.close()
    return {'y0':y0, 'avalue':avalue, 'psi0':psi0.reshape(30, 41)}


def vernerRead():
    """
    Reads the photoionization cross-section data from [6]_.

    Returns
    -------
    {'pqn','l','eth','e0','sig0','ya','p', yw'} : `dict`
        `pqn` is the principal quantum number, `l` is the subshell orbital quantum number, `eth` (in eV) is the subshell ionization threshold energy; `sig0`, `ya`, `p`, and `yw` are all fit parameters used in calculating the total photoionization cross-section.

    References
    ----------
    .. [6] Verner & Yakovlev, 1995, A&AS, `109, 125 <http://adsabs.harvard.edu/abs/1995A%26AS..109..125V>`_
    """
    xuvtop = os.environ['XUVTOP']
    fname = os.path.join(xuvtop, 'continuum', 'verner_short.txt')
    input = open(fname)
    lines = input.readlines()
    input.close()
    #
    nlines = 465
    maxZ = 30+1
    maxNel = 30 +1# also equal max(stage)
    #
    #z = np.array(nlines,'int32')
    #nel = np.array(nlines,'int32')
    pqn = np.zeros((maxZ,maxNel),'int32')
    l = np.zeros((maxZ,maxNel),'int32')
    eth = np.zeros((maxZ,maxNel),np.float64)
    e0 = np.zeros((maxZ,maxNel),np.float64)
    sig0 = np.zeros((maxZ,maxNel),np.float64)
    ya = np.zeros((maxZ,maxNel),np.float64)
    p = np.zeros((maxZ,maxNel),np.float64)
    yw = np.zeros((maxZ,maxNel),np.float64)
    #
    fstring = 'i2,i3,i2,i2,6f11.3'
    header_line = FortranRecordReader(fstring)
#    vernerFormat=FortranFormat(fstring)
    #
    for iline in range(nlines):
#        out=FortranLine(lines[iline],vernerFormat)
        out = header_line.read(lines[iline])
        z = out[0]
        nel = out[1]
        stage = z - nel + 1
        pqn[z,stage] = out[2]
        l[z,stage] = out[3]
        eth[z,stage] = out[4]
        e0[z,stage] = out[5]
        sig0[z,stage] = out[6]
        ya[z,stage] = out[7]
        p[z,stage] = out[8]
        yw[z,stage] = out[9]
    #
    return {'pqn':pqn, 'l':l, 'eth':eth, 'e0':e0, 'sig0':sig0, 'ya':ya, 'p':p, 'yw':yw}


def versionRead():
    """
    Read the version number of the CHIANTI database
    """
    xuvtop = os.environ['XUVTOP']
    vFileName = os.path.join(xuvtop, 'VERSION')
    vFile = open(vFileName)
    versionStr = vFile.readline()
    vFile.close()
    return versionStr.strip()


def wgfaRead(ions, filename=None, elvlcname=0, total=False, verbose=False):
    """
    Read CHIANTI data from a .wgfa file.

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V

    Keyword Arguments
    -----------------

    filename : `str`
        Custom filename, will override that specified by `ions`
    elvlcname : `str`
        If specified, the lsj term labels are returned in the 'pretty1' and 'pretty2'
        keys of 'Wgfa' dict
    total : `bool`
        Return the summed level 2 avalue data in 'Wgfa'
    verbose : `bool`

    Returns
    -------
    Wgfa : `dict`
        Information read from the .wgfa file. The dictionary structure is
        {"lvl1","lvl2","wvl","gf","avalue","ref","ionS","filename"}

    See Also
    --------
    ChiantiPy.tools.archival.wgfaRead : Read .wgfa file with the old format.
    """
    #
    if filename:
        wgfaname = filename
        if not elvlcname:
            elvlcname = os.path.splitext(wgfaname)[0] + '.elvlc'
            if os.path.isfile(elvlcname):
                elvlc = elvlcRead('', elvlcname)
            else:
                elvlc = 0
        else:
            elvlc = elvlcRead('',elvlcname)

    else:
        fname = util.ion2filename(ions)
        wgfaname = fname+'.wgfa'
        elvlcname = fname + '.elvlc'
        if os.path.isfile(elvlcname):
            elvlc = elvlcRead('', elvlcname)
        else:
            elvlc = 0
    if verbose:
        if elvlc:
            print(' have elvlc data')
        else:
            print(' do not have elvlc data')
    #
    input = open(wgfaname,'r')
    s1 = input.readlines()
    input.close()
    nwvl = 0
    ndata = 2
    while ndata > 1:
        s1a = s1[nwvl]
        s2 = s1a.split()
        ndata = len(s2)
        nwvl += 1
    nwvl -= 1
    if verbose:
        print((' nwvl = %10i ndata = %4i'%(nwvl, ndata)))
    lvl1 = [0]*nwvl
    lvl2 = [0]*nwvl
    wvl = [0.]*nwvl
    gf = [0.]*nwvl
    avalue = [0.]*nwvl
    if elvlc:
        pretty1  =  ['']*nwvl
        pretty2  =  ['']*nwvl
    #
    if verbose:
        print((' nwvl  =  %10i'%(nwvl)))
    #
    wgfaFormat = '(2i5,f15.3,2e15.3)'
    header_line = FortranRecordReader(wgfaFormat)
    for ivl in range(nwvl):
        if verbose:
            print(' index %5i  %s'%(ivl, s1[ivl]))
#        inpt=FortranLine(s1[ivl],wgfaFormat)
        inpt = header_line.read(s1[ivl])
        lvl1[ivl] = inpt[0]
        lvl2[ivl] = inpt[1]
        wvl[ivl] = inpt[2]
        gf[ivl] = inpt[3]
        avalue[ivl] = inpt[4]
        if elvlc:
            idx1 = elvlc['lvl'].index(inpt[0])
            idx2 = elvlc['lvl'].index(inpt[1])
            pretty1[ivl] = elvlc['pretty'][idx1]
            pretty2[ivl] = elvlc['pretty'][idx2]

    ref = []
    # should skip the last '-1' in the file
    for i in range(nwvl+1,len(s1)):
        s1a = s1[i]
        ref.append(s1a.strip())
    newref = [aref.strip() for aref in ref if aref.strip() != '-1']
    Wgfa = {"lvl1":lvl1, "lvl2":lvl2, "wvl":wvl, "gf":gf, "avalue":avalue, "ref":newref, 'ionS':ions,
        'filename':wgfaname}
    if total:
        avalueLvl = [0.]*max(lvl2)
        for iwvl in range(nwvl):
            avalueLvl[lvl2[iwvl] -1] += avalue[iwvl]
        Wgfa['avalueLvl'] = np.asarray(avalueLvl)

    if elvlc:
        Wgfa['pretty1'] = pretty1
        Wgfa['pretty2'] = pretty2
    #
    return Wgfa


def wgfaWrite(info, filename = None, minBranch = 1.e-5, rightDigits = 4, maxLvl1 = None,  comment=None):
    """
    Write data to a CHIANTI .wgfa file

    Parameters
    ----------
    info : `dict`
        Should contain the following keys:
        ionS, the Chianti style name of the ion such as c_4 for C IV,
        lvl1, the lower level, the ground level is 1,
        lvl2, the upper level,
        wvl, the wavelength (in Angstroms),
        gf,the weighted oscillator strength,
        avalue, the A value,
        pretty1, descriptive text of the lower level (optional),
        pretty2, descriptive text of the upper level (optiona),
        ref, reference text, a list of strings

    filename : `str`

    minBranch : `~numpy.float64`
        The transition must have a branching ratio greater than the specified minBranch
            to be written to the file.  default = 1.e-5

    rightDigits :  `int`
        the number of digits to the right of the decimal point in the wavelength
        default = 5

    maxLvl1 : `int`
        the largest level to be written. default is None

    comment :  `str`
        add a comment to the reference section.  default = None
    """
    #
#    gname = info['ionS']
    if filename:
        wgfaname = filename
    else:
        print(' output filename not specified, no file will be created')
        return
#        wgfaname = gname + '.wgfa'
    print((' wgfa file name = ', wgfaname))
    if minBranch > 0.:
        info['ref'].append(' minimum branching ratio = %10.2e'%(minBranch))
    if comment is not None:
        info['ref'].append(comment)

    info['ref'].extend(chRef)
    with open(wgfaname, 'w') as outpt:
        #ntrans = len(info['lvl1'])
        nlvl = max(info['lvl2'])
        totalAvalue = np.zeros(nlvl, np.float64)
        if 'pretty1' in info:
            pformat = '%5i%5i%15.' + str(rightDigits) + 'f%15.3e%15.3e%30s - %30s'
        else:
            pformat = '%5i%5i%15.' + str(rightDigits) + 'f%15.3e%15.3e'
        for itrans, avalue in enumerate(info['avalue']):
            # for autoionization transitions, lvl1 can be less than zero
            if abs(info['lvl1'][itrans]) > 0 and info['lvl2'][itrans] > 0:
                totalAvalue[info['lvl2'][itrans] -1] += avalue

        for itrans, avalue in enumerate(info['avalue']):
            if info['wvl'][itrans] == 0.:
                branch = 1.
            elif avalue > 0.:
                branch = avalue/totalAvalue[info['lvl2'][itrans] -1]
            else:
                branch = 0.
            test1 = branch > minBranch
            test2 = abs(info['lvl1'][itrans]) > 0
            test3 = info['lvl2'][itrans] > 0
            if maxLvl1:
                test4 = info['lvl1'][itrans] <= maxLvl1
            else:
                test4 = True
            if test1 and test2 and test3 and test4:
                if 'pretty1' in info:
                    # generally only useful with NIST data
                    if 'transType' in info:
                        if info['transType'][itrans] != '':
                            lbl2 = info['pretty2']+'  ' + info['transType'][itrans]
                    else:
                        lbl2 =  info['pretty2'][itrans]
                    pstring = pformat%(info['lvl1'][itrans], info['lvl2'][itrans], info['wvl'][itrans], info['gf'][itrans], avalue, info['pretty1'][itrans].rjust(30), lbl2.ljust(30))
                    outpt.write(pstring+'\n')
                else:
                    pstring = pformat%(info['lvl1'][itrans], info['lvl2'][itrans], info['wvl'][itrans], info['gf'][itrans], avalue)
                    outpt.write(pstring+'\n')
        outpt.write(' -1 \n')
        # if filename in ref, do not repeat
        cnt = 0
        for one in info['ref']:
            if 'filename' in one:
                cnt += 1
        if cnt == 0:
            outpt.write('%filename:  ' + wgfaname + '\n')
        ref = info['ref']
        newref = [aref.strip() for aref in ref if aref.strip() != '-1']
        for one in newref:
            outpt.write(one+'\n')
        outpt.write(today.strftime('%Y %B %d') +'\n')
        outpt.write(' -1 \n')

