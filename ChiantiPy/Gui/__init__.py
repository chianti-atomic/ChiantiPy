'''
trying to get a correct gui package
'''

#
hasWx=0
hasPyQt4=0
try: 
    import PyQt4
    hasPyQt4 = 1
    print(' found PyQt4 widgets')
    del PyQt4
except:
    try:
        import wx
        hasWx = 1
        print(' found Wx widgets')
        del wx
    except:
        print(' using cli')
#
if hasPyQt4:
    from .gui_qt import gui 
    print(' using PyQt4 widgets')
elif hasWx:
    from .gui_wx import gui
    print(' using Wx widgets')
else:
    from .gui_cl import gui 
    print(' using CLI for selections')

