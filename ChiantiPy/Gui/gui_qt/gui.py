'''
PyQt4 widget selection dialogs
'''
import sys,  os
from PyQt4 import QtGui
import ChiantiPy
from ChiantiPy.Gui.gui_qt.ui import *
' qt4 selection dialogs'

def chpicker(dir, filter='*.*', label='ChiantiPy'):
    '''Select a filename using a Qt gui dialog.'''
    app=QtGui.QApplication(sys.argv)
    a=QtGui.QFileDialog()
    a.setDirectory(dir)
    a.setFilter(filter)
#    mylabel=QtCore.QString('some label')
#    a.setLabelText(mylabel)
    a.setWindowTitle(label)
    a.setModal(True)
    a.exec_()
    qfilename=a.selectedFiles()
    return str(qfilename[0])
    #
#
class selectorDialog(QtGui.QDialog):
    '''Make a single or multiple selection from a list of items.

    expects the input of an array of items, will select one or more
    '''
    def __init__(self, items, label=None ,  parent=None, multiChoice=True):
#       if using the Qt4Agg backend for matplotlib, the following line needs to be comment out
#        app=QtGui.QApplication(sys.argv)
        QtGui.QDialog.__init__(self)
        self.ui = Ui_selectorDialogForm()
        self.ui.setupUi(self)
        if multiChoice:
            self.ui.listWidget.setSelectionMode(QtGui.QListWidget.MultiSelection)
        else:
            self.ui.listWidget.setSelectionMode(QtGui.QListWidget.SingleSelection)
        if label == None:
            self.setWindowTitle('ChiantiPy')
        else:
            self.setWindowTitle('ChiantiPy - '+label)
        imagefile = os.path.join(ChiantiPy.__path__[0], "images/chianti2.png")
        self.setWindowIcon(QtGui.QIcon(imagefile))
        for anitem in items:
#            print ' item = ', anitem, QtCore.QString(anitem)
            self.ui.listWidget.addItem(str(anitem))
        self.exec_()

    def accept(self):
#        print ' selector button pushed'
        nitems = self.ui.listWidget.count()
#        print ' nitems = ', nitems
        self.selectedIndex=[]
        self.selectedText=[]
        for i in range(nitems):
            anitem = self.ui.listWidget.item(i)
#            print 'selected? = ', anitem.isSelected()
            if anitem.isSelected():
#                print ' item = ' , str(anitem.text())
                self.selectedText.append(str(anitem.text()))
                self.selectedIndex.append(i)
        self.done(1)

    def reject(self):
#        print ' cancel button pushed'
        self.selectedIndex = None
        self.selectedText = None
        self.done(1)
    #
#from choice2DialogForm import *
#
class choice2Dialog(QtGui.QDialog):
    '''Make a single or multiple selection from a list of items and another
    single or multiple selection from the same list.

    Useful for picking numerators and denominators.

    expects the input of an array of items, will select one or more from both widgets.'''
    def __init__(self, items,  label=None ,  parent=None):
#       if using the Qt4Agg backend for matplotlib, the following line needs to be comment out
#        app=QtGui.QApplication(sys.argv)
        QtGui.QDialog.__init__(self)
#        app=QtGui.QApplication(sys.argv)
        self.ui = Ui_choice2DialogForm()
        self.ui.setupUi(self)
        if label == None:
            self.setWindowTitle('ChiantiPy')
        else:
            self.setWindowTitle('ChiantiPy - '+label)
        self.setWindowIcon(QtGui.QIcon('images/chianti2.png'))
        for anitem in items:
#            print ' item = ', anitem, QtCore.QString(anitem)
            self.ui.numListWidget.addItem(str(anitem))
            self.ui.denListWidget.addItem(str(anitem))
        self.exec_()
        #
    def accept(self):
        nitems = self.ui.numListWidget.count()
        self.numIndex=[]
        self.numText=[]
        for i in range(nitems):
            anitem = self.ui.numListWidget.item(i)
#            print 'selected? = ', anitem.isSelected()
            if anitem.isSelected():
#                print ' item = ' , str(anitem.text())
                self.numText.append(str(anitem.text()))
                self.numIndex.append(i)
        self.denIndex=[]
        self.denText=[]
        for i in range(nitems):
            anitem = self.ui.denListWidget.item(i)
#            print 'selected? = ', anitem.isSelected()
            if anitem.isSelected():
#                print ' item = ' , str(anitem.text())
                self.denText.append(str(anitem.text()))
                self.denIndex.append(i)
        self.done(1)

    def reject(self):
        print(' cancel button pushed')
        self.done(1)
