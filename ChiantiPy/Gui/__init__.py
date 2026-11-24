"""
Select GUI package
"""

#
import os
import configparser

initDefaults = {'abundfile': 'sun_photospheric_2021_asplund','ioneqfile': 'chianti', 'wavelength': 'angstrom',
    'flux': 'energy','gui':True}
#check chiantirc for gui selection
#rcfile1=os.path.join(os.environ['HOME'],'.config/chiantirc')
#rcfile2=os.path.join(os.environ['HOME'],'.chianti/chiantirc')
#rcparse=configparser.ConfigParser()
#rcparse.read(rcfile)
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

    print((' chiantirc file can not be found'))
    print(' should be in $HOME/.config, $HOME/.chianti or $PROFILEHOME/.chianti')
    print((' using the following defaults'))
    for akey in list(defaults.keys()):
        print((' %s = %s'%(akey, defaults[akey])))


#check for available gui
hasPyQt6=False
try:
    import PyQt6
    hasPyQt6 = True
    print(' found PyQt6 widgets')
    del PyQt6
except ImportError:
    print(' did not find PyQt6')
#set gui
if hasPyQt6 and defaults['gui']:
    from .gui_qt6 import gui
    print(' using PyQt6 widgets')
else:
    from .gui_cl import gui
    print(' using CLI for selections')
