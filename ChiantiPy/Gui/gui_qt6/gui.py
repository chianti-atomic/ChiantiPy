'''
PyQt6 widget selection dialogs
'''
import os
from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QLabel,
    QListWidget,
    QVBoxLayout,
    QAbstractItemView,
    QWidget
)
from PyQt6.QtGui import QFont

#import ChiantiPy
#from ChiantiPy.Gui.gui_qt6.ui import *
''' qt6 selection dialogs
'''

class chpicker(QWidget):
    ''' dialog to select a single file name under the directory
    code largely taken from pythonspot.com
    '''
    def __init__(self, dir, label):
        super().__init__()
        self.title = 'PyQt5 file dialogs - pythonspot.com'
        self.left = 10
        self.top = 10
        self.width = 740  # was 640
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

class selectorDialog(QDialog):
    def __init__(self, items, label=None,  parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select or Enter Item")
        self.resize(640, 480)

        # 1. Create the widgets
        if label is not None:
            self.label = QLabel(label)
        else:
            self.label = QLabel("Select a line or use ctrl c to select several")
#        self.text_box = QLineEdit()
#        self.text_box.setPlaceholderText("Type something here...")

        self.list_widget = QListWidget()
        self.list_widget.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection)
        font = QFont("DejaVu Sans Mono", 14)
        self.list_widget.setFont(font)
        self.list_widget.addItems(items)

        # Optional: Update text box when a list item is clicked
#        self.list_widget.currentTextChanged.connect(self.text_box.setText)

        # 2. Create Accept (OK) and Reject (Cancel) buttons
        self.button_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)

        # 3. Arrange widgets in a vertical layout
        layout = QVBoxLayout()
        layout.addWidget(self.label)
#        layout.addWidget(self.text_box)
        layout.addWidget(self.list_widget)
        layout.addWidget(self.button_box)
        self.setLayout(layout)
        self.exec()

    def get_values(self):
        """Helper method to retrieve the entered text and selected list item."""
#        selected_item = self.list_widget.currentItem()
        selected_items = self.list_widget.selectedItems()
#        self.selected_items = selected_items
#        selected_texts = [item.text() for item in selected_items]
        outpt = {}
        for i,  anitem in enumerate(selected_items):
            outpt['%i'%(i)] = anitem.text()
        return outpt

