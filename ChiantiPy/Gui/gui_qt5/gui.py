'''
PyQt5 widget selection dialogs
'''
import os
from PyQt5 import QtGui, QtWidgets

#from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog
#from PyQt5.QtGui import QIcon

import ChiantiPy
from ChiantiPy.Gui.gui_qt5.ui import *
''' qt5 selection dialogs
'''

class chpicker(QtWidgets.QWidget):
    ''' dialog to select a single file name under the directory
    code largely taken from pythonspot.com
    '''
    def __init__(self, dir, label):
        super().__init__()
        self.title = 'PyQt5 file dialogs - pythonspot.com'
        self.left = 10
        self.top = 10
        self.width = 640
        self.height = 480
        self.dir = dir
        self.label = label
        self.fileName = None
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.openFileNameDialog()
#        self.openFileNamesDialog()
#        self.saveFileDialog()

        self.show()
        self.close()

    def openFileNameDialog(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, self.label, self.dir,"All Files (*);;Python Files (*.py)", options=options)
#        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", self.dir,"All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            print(fileName)
        self.fileName = fileName
        self.baseName = os.path.split(fileName)[1]
        self.rootName = os.path.splitext(self.baseName)[0]

#    def openFileNamesDialog(self):
#        options = QtWidgets.QFileDialog.Options()
#        options |= QtWidgets.QFileDialog.DontUseNativeDialog
#        files, _ = QtWidgets.QFileDialog.getOpenFileNames(self,"QFileDialog.getOpenFileNames()", "","All Files (*);;Python Files (*.py)", options=options)
#        if files:
#            print(files)

#    def saveFileDialog(self):
#        options = QFileDialog.Options()
#        options |= QFileDialog.DontUseNativeDialog
#        fileName, _ = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.txt)", options=options)
#        if fileName:
#            print(fileName)


#def chpicker(dir, filter='*.*', label='ChiantiPy'):
#    '''Select a filename using a Qt gui dialog.
#    '''
##    app=QtWidgets.QApplication(sys.argv)
#    a=QtWidgets.QFileDialog.getOpenFileName()
#    a.SetFileMode(1)
#    a.setDirectory(dir)
#    # a.setFilter(filter)
##    mylabel=QtCore.QString('some label')
##    a.setLabelText(mylabel)
#    a.setWindowTitle(label)
#    a.setModal(True)
#    a.exec_()
#    qfilename=a.selectedFiles()
#    return str(qfilename[0])
    #
#
class selectorDialog(QtWidgets.QDialog):
    '''Make a single or multiple selection from a list of items.

    expects the input of an array of items, will select one or more'''
    def __init__(self, items, label=None ,  parent=None, multiChoice=True):
#       if using the Qt4Agg backend for matplotlib, the following line needs to be comment out
#        app=QtGui.QApplication(sys.argv)
        QtWidgets.QDialog.__init__(self)
        self.ui = Ui_selectorDialogForm()
        self.ui.setupUi(self)
        if multiChoice:
            self.ui.listWidget.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
        else:
            self.ui.listWidget.setSelectionMode(QtWidgets.QListWidget.SingleSelection)
        if label is None:
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
class choice2Dialog(QtWidgets.QDialog):
    '''Make a single or multiple selection from a list of items and another
    single or multiple selection from the same list.

    Useful for picking numerators and denominators.

    expects the input of an array of items, will select one or more from both widgets.'''
    def __init__(self, items,  label=None ,  parent=None, multi=True):
#       if using the Qt5Agg backend for matplotlib, the following line needs to be comment out
#        app=QtGui.QApplication(sys.argv)
        QtWidgets.QDialog.__init__(self)
#        app=QtGui.QApplication(sys.argv)
        self.ui = Ui_choice2DialogForm()
        self.ui.setupUi(self)
        if label is None:
            self.setWindowTitle('ChiantiPy')
        else:
            self.setWindowTitle('ChiantiPy - '+label)
        self.setWindowIcon(QtGui.QIcon('images/chianti2.png'))
        if multi:
            self.ui.numListWidget.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
            self.ui.denListWidget.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
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
