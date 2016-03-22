'''
trying to get a correct gui package
'''

#use chiantirc options to set gui
import os
try:
    #Python 3
    import configparser
except ImportError:
    #Python 2
    import ConfigParser as configparser
    
rcfile=os.path.join(os.environ['HOME'],'.chianti/chiantirc')
rcparse=configparser.ConfigParser()
rcparse.read(rcfile)
try:
    if rcparse['chianti']['gui'].lower() == 'true':
        use_gui=True
    else:
        use_gui=False
except KeyError:
    #default to true if section/field don't exist
    use_gui=True

hasWx=False
hasPyQt4=False
try: 
    import PyQt4
    hasPyQt4 = True
    print(' found PyQt4 widgets')
    del PyQt4
except:
    try:
        import wx
        hasWx = True
        print(' found Wx widgets')
        del wx
    except:
        print(' using cli')
#
if hasPyQt4 and use_gui:
    from .gui_qt import gui 
    print(' using PyQt4 widgets')
elif hasWx and use_gui:
    from .gui_wx import gui
    print(' using Wx widgets')
else:
    from .gui_cl import gui 
    print(' using CLI for selections')

