"""
Functions for reading pre-v8 CHIANTI files
"""

import os

#from .FortranFormat import *
from  ChiantiPy.fortranformat import FortranRecordReader
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.util as util
    #
    # -------------------------------------------------------------------------------------
    #
def elvlcRead(ions, filename = 0, verbose=0,  useTh=0):
    """
    read a chianti energy level file and returns
    {"lvl":lvl,"conf":conf,"term":term,"spin":spin,"l":l,"spd":spd,"j":j
    ,"mult":mult,"ecm":ecm,"eryd":eryd,"ecmth":ecmth,"erydth":erydth,"ref":ref,"pretty":pretty, 'ionS':ions}
    if a energy value for ecm or eryd is zero(=unknown), the theoretical values
    (ecmth and erydth) are inserted
    """
    #
#    (i3,i6,a15,2i3,a3,f4.1,i3,f15.3,f15.6,f15.3,f15.6,f15.3,f15.6)
#    fstring='i3,i6,a15,i3,i3,a3,f4.1,i3,4f15.2'
#    elvlcFormat=FortranFormat(fstring)
#    header_line = FortranRecordReader('i3,i6,a15,i3,i3,a3,i3,f4.1,f15.3,f15.6,f15.3,f15.6')
    header_line = FortranRecordReader('i3,i6,a15,i3,i3,a3,f4.1,i3,f15.3,f15.6,f15.3,f15.6') #',f4.1')

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
    conf = [0]*nlvls
    term = [0]*nlvls
    spin = [0]*nlvls
    l = [0]*nlvls
    spd = [0]*nlvls
    j = [0]*nlvls
    mult = [0]*nlvls
    ecm = [0]*nlvls
    eryd = [0]*nlvls
    ecmth = [0]*nlvls
    erydth = [0]*nlvls
    pretty = [0]*nlvls
    label = []
    for i in range(0,nlvls):
        if verbose:
            print((s1[i][0:115]))
#        inpt = FortranLine(s1[i][0:115],elvlcFormat)
        inpt  =  header_line.read(s1[i])
        lvl[i] = inpt[0]
        conf[i] = inpt[1]
        label.append(str(inpt[1]))
        term[i] = inpt[2].strip()
        spin[i] = inpt[3]
        l[i] = inpt[4]
        spd[i] = inpt[5].strip()
        j[i] = inpt[6]
        mult[i] = inpt[7]
        ecm[i] = inpt[8]
        eryd[i] = inpt[9]
        ecmth[i] = inpt[10]
        erydth[i] = inpt[11]
        if ecm[i] == 0.:
            if useTh:
                ecm[i] = ecmth[i]
                eryd[i] = erydth[i]
        stuff  =  term[i].strip() + ' %1i%1s%3.1f'%( spin[i], spd[i], j[i])
        pretty[i] = stuff.strip()
    ref = []
    for i in range(nlvls+1,len(s1)-1):
        s1a = s1[i][:-1]
        ref.append(s1a.strip())
#    self.const.Elvlc = {"lvl":lvl,"conf":conf,"term":term,"spin":spin,"l":l,"spd":spd,"j":j
#            ,"mult":mult,"ecm":ecm,"eryd":eryd,"ecmth":ecmth,"erydth":erydth,"ref":ref}
    return {"lvl":lvl,"conf":conf,"label":label,"term":term,"spin":spin,"l":l,"spd":spd,"j":j
            ,"mult":mult,"ecm":ecm,"eryd":eryd,"ecmth":ecmth,"erydth":erydth,"ref":ref,"pretty":pretty, 'ionS':ions, 'status':status}
    #
    # -------------------------------------------------------------------------------------
    #
def elvlcWrite(info, outfile=None, addLvl=0):
    '''
    Write Chianti data to .elvlc file.

    Parameters
    ----------
    info : `dict`
        Information about the Chianti data to write. Should contain
        the following keys: ionS, the Chianti style name of the ion such as c_4
        conf, an integer denoting the configuration - not too essential
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
    addLvl : `int`
        Add a constant value to the index of all levels

    Notes
    -----
    For use with files created before elvlc format change in November 2012

    See Also
    --------
    ChiantiPy.tools.io.elvlcWrite : Write .elvlc file using the new format.
    '''
    if outfile:
        elvlcName = outfile
    else:
        elvlcName = info['ionS'] + '.elvlc'
    print((' elvlc file name = %s'%(elvlcName)))
    out = open(elvlcName, 'w')
    for i,  conf in enumerate(info['conf']):
        mult = int(2.*info['j'][i]+1.)
        thisTerm = info['term'][i].ljust(14)
        pstring = '%3i%6s%15s%3i%3i%2s%5.1f%3i%15.3f%15.6f%15.3f%15.6f \n'%(i+1+addLvl, conf, thisTerm, info['spin'][i], info['l'][i], info['spd'][i], info['j'][i], mult, info['ecm'][i], info['eryd'][i], info['ecmth'][i], info['erydth'][i])
    #i3,a6,a15,2i3,a2,f5.1,i3,f15.3,f15.6,f15.3,f15.6
        out.write(pstring)
    out.write(' -1\n')
    out.write('%filename:  ' + elvlcName + '\n')
    for one in info['ref']:
        out.write(one+'\n')
    out.write(' -1\n')
    out.close()

def wgfaRead(ions, filename=None, elvlcname=-1, total=False, verbose=False):
    """
    Read CHIANTI data from a .wgfa file.

    Parameters
    ----------
    ions : `str`
        Ion, e.g. 'c_5' for C V
    filename : `str`
        Custom filename, will override that specified by `ions`
    elvlcname : `str`
        If specified, the lsj term labels are returned in the `pretty1` and `pretty2` keys of `Wgfa`
    total : `bool`
        Return the level 2 avalue data in `Wgfa`
    verbose : `bool`

    Returns
    -------
    Wgfa : `dict`
        Information read from the .wgfa file. The dictionary structure is {"lvl1","lvl2","wvl","gf","avalue","ref","ionS","filename"}

    Notes
    -----
    This is text-wise not different than the v8 version except that it uses the
    archival elvlcRead in `~ChiantiPy.tools.archival` though this has now been commented out. Can this routine be removed? Should the elvlcRead routine be uncommented?

    See Also
    --------
    ChiantiPy.tools.io.wgfaRead : Read .wgfa file with the new format.
    """
    #
    if filename:
        wgfaname = filename
        if elvlcname < 0:
            elvlcname = 0
            elvlc = 0
        elif not elvlcname:
            elvlcname = os.path.splitext(wgfaname)[0] + '.elvlc'
            if os.path.isfile(elvlcname):
                elvlc = elvlcRead('', elvlcname)
            else:
                elvlc = 0
        else:
            elvlc = elvlcRead('',elvlcname)

    else:
        fname=util.ion2filename(ions)
        wgfaname=fname+'.wgfa'
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
    input=open(wgfaname,'r')
    s1=input.readlines()
    input.close()
    nwvl=0
    ndata=2
    while ndata > 1:
        s1a=s1[nwvl]
        s2=s1a.split()
        ndata=len(s2)
        nwvl += 1
    nwvl -= 1
    if verbose:
        print(' nwvl = %10i ndata = %4i'%(nwvl, ndata))
    lvl1=[0]*nwvl
    lvl2=[0]*nwvl
    wvl=[0.]*nwvl
    gf=[0.]*nwvl
    avalue=[0.]*nwvl
    if elvlc:
        pretty1 = ['']*nwvl
        pretty2 = ['']*nwvl
    #
    if verbose:
        print(' nwvl = %10i'%(nwvl))
    #
    wgfaFormat='(2i5,f15.3,2e15.3)'
    for ivl in range(nwvl):
        inpt=FortranLine(s1[ivl],wgfaFormat)
        lvl1[ivl]=inpt[0]
        lvl2[ivl]=inpt[1]
        wvl[ivl]=inpt[2]
        gf[ivl]=inpt[3]
        avalue[ivl]=inpt[4]
        if elvlc:
            pretty1[ivl] = elvlc['pretty'][inpt[0] - 1]
            pretty2[ivl] = elvlc['pretty'][inpt[1] - 1]

    ref=[]
    # should skip the last '-1' in the file
    for i in range(nwvl+1,len(s1) -1):
        s1a=s1[i][:-1]
        ref.append(s1a.strip())
    Wgfa={"lvl1":lvl1,"lvl2":lvl2,"wvl":wvl,"gf":gf,"avalue":avalue,"ref":ref, 'ionS':ions, 'filename':wgfaname}
    if total:
        avalueLvl = [0.]*max(lvl2)
        for iwvl in range(nwvl):
            avalueLvl[lvl2[iwvl] -1] += avalue[iwvl]
        Wgfa['avalueLvl'] = avalueLvl

    if elvlc:
        Wgfa['pretty1'] = pretty1
        Wgfa['pretty2'] = pretty2
    #
    return Wgfa
    #
    # --------------------------------------
    #
