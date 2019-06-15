"""
Reading and writing functions
"""
import os
from datetime import date
import fnmatch
import pickle
import configparser

import numpy as np

import ChiantiPy.tools.util as util
import ChiantiPy.tools.constants as const
import ChiantiPy.Gui as chgui
from  ChiantiPy.fortranformat import FortranRecordReader

today = date.today()
def abundanceRead(abundancename=''):
    """
    Read abundance file `abundancename` and return the abundance values relative to hydrogen
    """
    abundance = np.zeros((50),np.float64)
    xuvtop = os.environ["XUVTOP"]
    if abundancename:
        # a specific abundance file name has been specified
        abundancefile = os.path.join(xuvtop,'abundance',abundancename+'.abund')
    else:
        # the user will select an abundance file
        abundir = os.path.join(xuvtop,'abundance')
        abundlabel = 'ChiantiPy - Select an abundance file'
        #fname = chianti.gui.chpicker(abundir, filter='*.abund', label=abundlabel)
        fname = chgui.gui.chpicker(abundir, filter='*.abund', label=abundlabel)
        if fname is None:
            print((' no abundance file selected'))
            return 0
        else:
            abundancefile = os.path.join(abundir, fname)
            abundancefilename = os.path.basename(fname)
            abundancename,ext = os.path.splitext(abundancefilename)
#    else:
#        # the default abundance file will be used
#        abundancename=self.Defaults['abundfile']
#        fname=os.path.join(xuvtop,'abundance',abundancename+'.abund')
    input = open(abundancefile,'r')
    s1 = input.readlines()
    input.close()
    nlines = 0
    idx = -1
    while idx <= 0:
        minChar = min([5, len(s1[nlines])])
        aline = s1[nlines][0:minChar]
        idx = aline.find('-1')
        nlines += 1
    nlines -= 1
    for line in range(nlines):
        z,ab,element = s1[line].split()
        abundance[int(z)-1] = float(ab)
    gz = np.nonzero(abundance)
    abs = 10.**(abundance[gz]-abundance[0])
    abundance.put(gz,abs)
    abundanceRef = s1[nlines+1:]
    return {'abundancename':abundancename,'abundance':abundance,'abundanceRef':abundanceRef}


def zion2name(z,ion, dielectronic=False):
    """
    Convert `Z` and `ion` to generic name, e.g. 26, 13 -> fe_13

    Parameters
    ----------
    z : `int`
        the nuclear charge, for example 26 for Fe XIV
    ion : `int`
        the ion stage, for example, 14 for Fe XIV


    Notes
    -----
    A duplicate of the routine in `ChiantiPy.tools.util` but needed by masterList Info
    TODO: Put in separate module to avoid multiple copies
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


    Notes
    -----
    A duplicate of the routine in `ChiantiPy.tools.util` but needed by masterList Info
    TODO: Put in separate module to avoid multiple copies
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
    input = open(autoname,'r')
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
        for iwvl in range(nwvl):
            avalueLvl[lvl2[iwvl] -1] += avalue[iwvl]
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


def autoWrite(info, outfile = None, minBranch = None):
    """
    Write data to a CHIANTI .wgfa file

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
    if outfile:
        autoname = outfile
    else:
        print(' output filename not specified, no file will be created')
        return
    print((' auto file name = ', autoname))
    if minBranch == None:
        minBranch = 0.
    else:
        info['ref'].append(' minimum branching ratio = %10.2e'%(minBranch))
    out = open(autoname, 'w')
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
    out.write('%filename:  ' + autoname + '\n')
    for one in info['ref']:
        out.write(one+'\n')
    out.close()

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
    initDefaults = {'abundfile': 'sun_photospheric_2015_scott','ioneqfile': 'chianti', 'wavelength': 'angstrom', 'flux': 'energy','gui':False}
    rcfile = os.path.join(os.environ['HOME'],'.chianti/chiantirc')
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
    else:
        defaults = initDefaults
        if verbose:
            print((' chiantirc file (/HOME/.chianti/chiantirc) does not exist'))
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
    hdr = input.readlines()
    input.close()
    info = {"iz":iz,"ion":ion,"nspl":nspl,"neaev":neaev, 'nfac':nfac}
    if neaev:
        info['eaev'] = eaev
    DiParams = {"info":info,"btf":btf,"ev1":ev1,"xsplom":xsplom,"ysplom":ysplom, 'eaev':eaev,"ref":hdr}
    return DiParams


def drRead(ions, filename=None):
    """
    Read CHIANTI dielectronic recombination .drparams files
    if filename is set, then reads that file
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
            if len(inpt)  > 13:
                nspl[i] = 9
                splups[i].put(list(range(9)),inpt[6:])
            else:
                nspl[i] = 5
                splups[i].put(list(range(5)),inpt[6:])
        #
        ref = []
        for i in range(nsplups+1,len(s1)):
            s1a = s1[i][:-1]
            ref.append(s1a.strip())
    return {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
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
    info = {"lvl":lvl,"conf":conf, "term":term,'label':label, "spin":spin, "spd":spd, "l":l, "j":j,
             'mult':mult, "ecm":ecm, 'eryd':eryd,'erydth':erydth, "ecmth":ecmth, "ref":ref,
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
        Output filename. ionS+'.elvlc' (in current directory) if None
    round : `int`
        input to `np.round' to round input values to maintain the correct number of significant figures
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
        info['eryd'] = [x*const.invCm2ryd  if x > 0. else  -1. for x in info['ecm']]
    if 'erydth 'not in info:
        info['erydth'] = [x*const.invCm2ryd  if x > 0. else  -1. for x in info['ecmth']]
    if 'eV' not in info:
        info['eV'] = [x*const.invCm2Ev if x > 0. else  -1. for x in info['ecm']]
    if 'eVth 'not in info:
        info['eVth'] = [x*const.invCm2Ev if x > 0. else  -1. for x in info['ecmth']]
   #
    out = open(elvlcName, 'w')
    for i,  aterm in enumerate(info['term']):
        thisTerm = aterm.ljust(29)
        thisLabel = info['label'][i].ljust(4)
#        print, ' len of thisTerm = ', len(thisTerm)
        pstring = '%7i%30s%5s%5i%5s%5.1f%15.3f%15.3f'%(i+1+addLvl, thisTerm, thisLabel, info['spin'][i], info['spd'][i],info['j'][i],  np.round(info['ecm'][i], round), np.round(info['ecmth'][i], round))
        if includeRyd:
             pstring += ' , %15.8f , %15.8f'%(info['eryd'][i], info['erydth'][i])
        if includeEv:
             pstring += ' , %15.8f , %15.8f'%(info['eV'][i], info['eVth'][i])
        pstring += '\n'
        out.write(pstring)
    out.write(' -1\n')
    out.write('%filename:  ' + os.path.split(elvlcName)[1] + '\n')
    for aref in info['ref']:
        out.write(aref + '\n')
    out.close()
    return


def fblvlRead(ions, filename=None, verbose=False):
    """
    Read a Chianti energy level file for calculating the
    free-bound continuum
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

    References
    ----------
    .. [1] Sutherland, R. S., 1998, MNRAS, `300, 321
        <http://adsabs.harvard.edu/abs/1998MNRAS.300..321S>`_
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
    Read in the free-free gaunt factors of [1]_.

    References
    ----------
    .. [1] Itoh, N. et al., 2000, ApJS, `128, 125
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


def klgfbRead():
    """
    Read CHIANTI files containing the free-bound gaunt factors for n=1-6 from [1]_.

    Returns
    -------
    {'pe', 'klgfb'} : `dict`
        Photon energy and the free-bound gaunt factors

    References
    ----------
    .. [1] Karzas and Latter, 1961, ApJSS, `6, 167
        <http://adsabs.harvard.edu/abs/1961ApJS....6..167K>`_
    """
    xuvtop = os.environ['XUVTOP']
    fname = os.path.join(xuvtop, 'continuum', 'klgfb.dat')
    input = open(fname)
    lines = input.readlines()
    input.close()
    #
    ngfb = int(lines[0].split()[0])
    nume = int(lines[0].split()[1])

    gfb = np.zeros((ngfb, ngfb, nume), np.float64)
    nlines = len(lines)
#        print 'nlines, nume, ngfb = ', nlines,  nume, ngfb
    pe = np.asarray(lines[1].split(), np.float64)
    for iline in range(2, nlines):
        data = lines[iline].split()
        n = int(data[0])
        l = int(data[1])
        gfb[n-1, l] = np.array(data[2:], np.float64)
    return {'pe':pe, 'klgfb':gfb}


def ioneqRead(ioneqName='', minIoneq=1.e-20, verbose=False):
    """
    Reads an ioneq file
    ionization equilibrium values less then minIoneq are returns as zeros
    Returns
    -------
    {'ioneqname','ioneqAll','ioneqTemperature','ioneqRef'} : `dict`
        Ionization equilibrium values and the reference to the literature
    """
    dir = os.environ["XUVTOP"]
    ioneqdir = os.path.join(dir,'ioneq')
    ioneqNames = util.listRootNames(ioneqdir)
    if ioneqName not in ioneqNames:
        # the user will select an ioneq file
        choice = chgui.gui.chpicker(ioneqdir, label='Select a single ioneq file')
        if choice.rootName in ioneqNames:
            fname = choice.fileName
            ioneqName = choice.rootName
#        fname1 = choice.baseName
#        fname1 = chgui.gui.chpicker(ioneqdir,filter='*.ioneq',label = 'Select an Ionization Equilibrium file')
#        fname = os.path.join(ioneqdir, fname1)
        if fname == None:
            print(' no ioneq file selected')
            return False
        else:
            ioneqfilename = os.path.basename(fname)
            ioneqname,ext = os.path.splitext(ioneqfilename)
    else:
        filelist = util.listFiles(ioneqdir)
        idx = ioneqNames.index(ioneqName)
        fname = filelist[idx]
#        fname = os.path.join(dir,'ioneq',ioneqname+'.ioneq')
#        newlist = fnmatch.filter(filelist, '*.ioneq')
#        baselist = []
#        for one in newlist:
#            baselist.append(os.path.basename(one))
#        cnt = baselist.count(ioneqname+'.ioneq')
#        if cnt == 0:
#            print((' ioneq file not found:  ', fname))
#            print(' the following files do exist: ')
#            for one in newlist:
#                print((os.path.basename(one)))
#            return
#        elif cnt == 1:
#            idx = baselist.index(ioneqname+'.ioneq')
#            if verbose:
#                print((' file exists:  ', newlist[idx]))
#            fname = newlist[idx]
#        elif cnt > 1:
#            print((' found more than one ioneq file', fname))
#            return
    #
    input = open(fname,'r')
    s1 = input.readlines()
    input.close()
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
#    ioneqformat=FortranFormat('2i3,'+str(nTemperature)+'e10.2')
    header_lineq = FortranRecordReader('2i3,'+str(nTemperature)+'e10.2')
    #
    ioneqAll = np.zeros((nElement,nElement+1,nTemperature),np.float64)
    for iline in range(2,nlines):
#        out=FortranLine(s1[iline],ioneqformat)
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
    ip = np.zeros((maxz, maxz), np.float64)
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
        haveStage = np.zeros((31, 31), 'Int32')
        haveDielectronic = np.zeros((31, 31), 'Int32')
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


def photoxRead(ions):
    """
    Read CHIANTI photoionization .photox files

    Returns
    -------
    {'lvl1', 'lvl2', 'energy', 'cross', 'ref'} : `dict`
        Energy (in Rydbergs) and cross section (in :math:`\mathrm{cm}^{-2}`)

    Notes
    -----
    The photox files are not in any released version of the CHIANTI database.
    """
    #
    zion = util.convertName(ions)
    if zion['Z'] < zion['Ion']:
        print((' this is a bare nucleus that has no ionization rate'))
        return
    #
    fname = util.ion2filename(ions)
    paramname = fname+'.photox'
    input = open(paramname,'r')
    lines = input.readlines()
    input.close
    # get number of energies
#    neng = int(lines[0][0:6])
    dataEnd = 0
    lvl1 = []
    lvl2 = []
    energy = []
    cross = []
    icounter = 0
    while not dataEnd:
        lvl11 = int(lines[icounter][:8])
        lvl21 = int(lines[icounter][8:15])
        ener = lines[icounter][15:].split()
        energy1 = np.asarray(ener, np.float64)
        #
        icounter += 1
        irsl = int(lines[icounter][:8])
        ind0 = int(lines[icounter][8:15])
        if irsl != lvl11 or ind0 != lvl21:
            # this only happens if the file was written incorrectly
            print((' lvl1, lvl2 = %7i %7i'%(lvl11, lvl21)))
            print((' irsl, indo = %7i %7i'%(irsl,  ind0)))
            return
        crs = lines[icounter][15:].split()
        cross1 = np.asarray(crs, np.float64)
        lvl1.append(lvl11)
        lvl2.append(lvl21)
        energy.append(energy1)
        cross.append(cross1)
        icounter += 1
        dataEnd = lines[icounter].count('-1')
    ref = lines[icounter+1:-1]
    cross = np.asarray(cross, np.float64)
    energy = np.asarray(energy, np.float64)
    return {'lvl1':lvl1, 'lvl2':lvl2,'energy':energy, 'cross':cross,  'ref':ref}


def rrRead(ions, filename=None):
    """
    Read CHIANTI radiative recombination .rrparams files

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
        ref = lines[3:-2]
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
    ''' to read the Mao 2017 rr loss parameters

    References
    ----------
    .. [1] Mao J., Kaastra J., Badnell N.R., `2017 Astron. Astrophys. 599, A10
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
    Read the new format v8 scups file containing the scaled temperature and upsilons from [1]_.

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V
    filename : `str`, optional
        Custom filename, will override that specified by `ions`
    verbose : `bool`

    References
    ----------
    .. [1] Burgess, A. and Tully, J. A., 1992, A&A, `254, 436 <http://adsabs.harvard.edu/abs/1992A%26A...254..436B>`_
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
    ea : `bool`
        Read .easplom file
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
    if type(filename) == type(None):
        fname = util.ion2filename(ions)
        if ea:
            splomname = fname+'.easplom'
        else:
            splomname = fname+'.splom'
    else:
        splomname = filename
    input = open(splomname,'r')
    #  need to read first line and see how many elements
    line1 = input.readline()
    #indices=line1[0:15]
    remainder = line1[16:]
    nom = remainder.split(' ')
#    format = FortranFormat('5i3,'+str(len(nom))+'E10.2')
    header_line = FortranRecordReader('5i3,'+str(len(nom))+'E10.2')
    #  go back to the beginning
    input.seek(0)
    lines = input.readlines()
    data = 5
    iline = 0
    lvl1 = []
    lvl2 = []
    ttype = []
    gf = []
    de = []
    f = []
    splom = []
    #ntrans = 0
    while data > 1:
#        splomdat = FortranLine(lines[iline],format)
        splomdat = header_line.read(lines[iline])
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
        data = len(lines[iline].split(' ',2))
    hdr = lines[iline+1:-1]
    de = np.asarray(de,np.float64)
    splomout = np.asarray(splom,np.float64)
    splomout = np.transpose(splomout)
    input.close()
    # note:  de is in Rydbergs
    splom = {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"deryd":de,"c":f
        ,"splom":splomout,"ref":hdr}
    return  splom


def splupsRead(ions, filename=None, filetype='splups'):
    """
    Read a CHIANTI .splups file

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V
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
            header_line3 = FortranRecordReader(str(nspl[i])+'e10.3' )
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

def trRead(ionS):
    ''' read the files containing total recombination rates .trparams
    '''
    stuff = util.convertName(ionS)
    filename = stuff['filename']
    trname = filename + '.trparams'
    if os.path.exists(trname):
        temperature = []
        rate = []
        inpt = open(trname)
        lines = inpt.readlines()
        ndata = int(lines[0])
        inpt.close()
        for jline in range(1, ndata+1):
            dummy = lines[jline].replace(os.linesep, '').split()
            temperature.append(float(dummy[0]))
            rate.append(float(dummy[1]))
        return {'temperature':np.asarray(temperature, np.float64), 'rate':np.asarray(rate, np.float64)}
    else:
        return 'file does not exist'



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
    Reads the photoionization cross-section data from [1]_.

    Returns
    -------
    {'pqn','l','eth','e0','sig0','ya','p', yw'} : `dict`
        `pqn` is the principal quantum number, `l` is the subshell orbital quantum number, `eth` (in eV) is the subshell ionization threshold energy; `sig0`, `ya`, `p`, and `yw` are all fit parameters used in calculating the total photoionization cross-section.

    References
    ----------
    .. [1] Verner & Yakovlev, 1995, A&AS, `109, 125 <http://adsabs.harvard.edu/abs/1995A%26AS..109..125V>`_
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
    Wgfa = {"lvl1":lvl1,"lvl2":lvl2,"wvl":wvl,"gf":gf,"avalue":avalue,"ref":ref, 'ionS':ions, 'filename':wgfaname}
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


def wgfaWrite(info, outfile = None, minBranch = 1.e-5, rightDigits = 4, maxLvl1 = None):
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
    outfile : `str`
    minBranch : `~numpy.float64`
        The transition must have a branching ratio greater than the specified minBranchto be written to the file
    """
    #
#    gname = info['ionS']
    if outfile:
        wgfaname = outfile
    else:
        print(' output filename not specified, no file will be created')
        return
#        wgfaname = gname + '.wgfa'
    print((' wgfa file name = ', wgfaname))
    if minBranch > 0.:
        info['ref'].append(' minimum branching ratio = %10.2e'%(minBranch))
    out = open(wgfaname, 'w')
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
                out.write(pstring+'\n')
            else:
                pstring = pformat%(info['lvl1'][itrans], info['lvl2'][itrans], info['wvl'][itrans], info['gf'][itrans], avalue)
                out.write(pstring+'\n')
    out.write(' -1\n')
    out.write('%filename:  ' + wgfaname + '\n')
    for one in info['ref']:
        out.write(one+'\n')
    out.write(today.strftime('%Y %B %d') +'\n')
    out.write(' -1 \n')
    out.close()
