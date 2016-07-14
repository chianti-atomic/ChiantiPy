'''
wxWidget selection dialogs.
'''

' wxWidget selection dialogs'

import wx
from Chianti.Gui.gui_wx.ui import *

def chpicker(dir,filter='All files (*.*)|*.*',label='ChiantiPy'):
    '''
    Select a filename using a wx gui dialog.
    '''
    app=wx.App()
    a=wx.FileDialog(None)
    a.SetMessage(label)
    a.SetWildcard(filter)
    a.SetDirectory(dir)
    a.ShowModal()
    a.name = a.GetFilename()
    a.Destroy()
    return a.name

class selectorDialog:
    '''
    Make a single or multiple selection from a list of items.

    expects the input of an array of items, will select one or more
    '''
    def __init__(self, items,label='Your list',title='Select One'):
        a = wx.App()
        a = wx.MultiChoiceDialog(None,label,title,items)
        a.ShowModal()
        self.selectedIndex = a.GetSelections()
        selectedItems=[]
        for one in self.selectedIndex:
            selectedItems.append(items[one])
        self.selectedItems = selectedItems
        a.Destroy

class choice2Dialog:
    '''Make a single or multiple selection from a list of items and another
    single or multiple selection from the same list.

    Useful for picking numerators and denominators.

    expects the input of an array of items, will select one or more from both widgets.
    '''
    def __init__(self, items,  label=None):
#       app = wx.PySimpleApp(0)
        app = wx.App()
        dlg = ui_choice2Dialog(None, -1, "")
        nItems = len(items)
        for i in range(nItems):
            dlg.numListBox.Insert(str(items[i]),i)
            dlg.denListBox.Insert(str(items[i]),i)
#       app.SetTopWindow(dialog_1)
        app.SetTopWindow(dlg)
        if dlg.ShowModal() == wx.ID_OK:
            self.numIndex = dlg.numListBox.GetSelections()
#           print ' numIndex = ', self.numIndex
            self.denIndex = dlg.denListBox.GetSelections()
            self.numText = []
            self.denText = []
            for one in self.numIndex:
#               print ' one = ', one
#               self.numText.append(items[self.numIndex[one]])
                self.numText.append(items[one])
            for one in self.denIndex:
#               self.denText.append(items[self.denIndex[one]])
                self.denText.append(items[one])
        dlg.Destroy()


