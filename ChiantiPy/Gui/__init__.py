'''
trying to get a correct gui package
'''

import os
try:
    #Python 3
    import configparser
except ImportError:
    #Python 2
    import ConfigParser as configparser

#check chiantirc for gui selection
rcfile=os.path.join(os.environ['HOME'],'.chianti/chiantirc')
rcparse=configparser.ConfigParser()
rcparse.read(rcfile)
try:
    if rcparse.get('chianti','gui').lower() == 'true':
        use_gui=True
    else:
        use_gui=False
except (KeyError,configparser.NoSectionError) as e:
    #default to true if section/field don't exist
    use_gui=True

#check for available gui
hasWx=False
hasPyQt4=False
try:
    import PyQt4
    hasPyQt4 = True
    print(' found PyQt4 widgets')
    del PyQt4
except ImportError:
    try:
        import wx
        hasWx = True
        print(' found Wx widgets')
        del wx
    except ImportError:
        print(' using cli')

#set gui
if hasPyQt4 and use_gui:
    from .gui_qt import gui
    print(' using PyQt4 widgets')
elif hasWx and use_gui:
    from .gui_wx import gui
    print(' using Wx widgets')
else:
    from .gui_cl import gui
    print(' using CLI for selections')
