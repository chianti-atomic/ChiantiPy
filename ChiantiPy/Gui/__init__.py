"""
Select GUI package
"""

#
import os
import configparser

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
hasPyQt5=False
try:
    import PyQt5
    hasPyQt5 = True
    print(' found PyQt5 widgets')
    del PyQt5
except ImportError:
    print(' using cli')

#set gui
if hasPyQt5 and use_gui:
    from .gui_qt5 import gui
    print(' using PyQt5 widgets')
else:
    from .gui_cl import gui
    print(' using CLI for selections')
