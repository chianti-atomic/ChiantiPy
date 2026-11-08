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
hasPyQt6=False
try:
    import PyQt6
    hasPyQt6 = True
    print(' found PyQt6 widgets')
    del PyQt6
except ImportError:
    print(' using cli')

#set gui
if hasPyQt6 and use_gui:
    from .gui_qt6 import gui
    print(' using PyQt6 widgets')
else:
    from .gui_cl import gui
    print(' using CLI for selections')
