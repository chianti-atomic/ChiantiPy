# Form implementation generated from reading ui file 'choice2DialogForm.ui'
#
# Created: Tue Mar 17 13:52:16 2009
#      by: PyQt4 UI code generator 4.4.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
#
class Ui_selectorDialogForm(object):
    def setupUi(self, selectorDialogForm):
        selectorDialogForm.setObjectName("selectorDialogForm")
        selectorDialogForm.setWindowModality(QtCore.Qt.WindowModal)
        selectorDialogForm.setModal(False)
#        selectorDialogForm.setWindowModality(QtCore.Qt.ApplicationModal)
        selectorDialogForm.resize(500,300)
        self.buttonBox = QtWidgets.QDialogButtonBox(selectorDialogForm)
        self.buttonBox.setGeometry(QtCore.QRect(30,240,341,32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.listWidget = QtWidgets.QListWidget(selectorDialogForm)
        self.listWidget.setGeometry(QtCore.QRect(70,20,360,200))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.listWidget.setFont(font)
        self.listWidget.setObjectName("listWidget")
        self.retranslateUi(selectorDialogForm)
        self.buttonBox.accepted.connect(selectorDialogForm.accept)
        self.buttonBox.rejected.connect(selectorDialogForm.reject)
#        QtCore.QMetaObject.connectSlotsByName(selectorDialogForm)

    def retranslateUi(self, selectorDialogForm):
        selectorDialogForm.setWindowTitle(QtWidgets.QApplication.translate("selectorDialogForm", "Dialog", None))

#
class Ui_choice2DialogForm(object):
    def setupUi(self, choice2DialogForm):
        choice2DialogForm.setObjectName("choice2DialogForm")
        choice2DialogForm.setWindowModality(QtCore.Qt.WindowModal)
        choice2DialogForm.resize(543,368)
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setWeight(75)
        font.setBold(True)
        choice2DialogForm.setFont(font)
        choice2DialogForm.setModal(False)
        self.widget = QtWidgets.QWidget(choice2DialogForm)
        self.widget.setGeometry(QtCore.QRect(10,10,522,341))
        self.widget.setObjectName("widget")
        self.gridLayout = QtWidgets.QGridLayout(self.widget)
        self.gridLayout.setSizeConstraint(QtWidgets.QLayout.SetMinimumSize)
        # self.gridLayout.setMargin(10)
        self.gridLayout.setObjectName("gridLayout")
        self.numLabel = QtWidgets.QLabel(self.widget)
        self.numLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.numLabel.setObjectName("numLabel")
        self.gridLayout.addWidget(self.numLabel,0,0,1,1)
        self.denLabel = QtWidgets.QLabel(self.widget)
        self.denLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.denLabel.setObjectName("denLabel")
        self.gridLayout.addWidget(self.denLabel,0,1,1,1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.numListWidget = QtWidgets.QListWidget(self.widget)
        # self.numListWidget.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.numListWidget.setObjectName("numListWidget")
        self.horizontalLayout.addWidget(self.numListWidget)
        self.denListWidget = QtWidgets.QListWidget(self.widget)
        # self.denListWidget.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.denListWidget.setObjectName("denListWidget")
        self.horizontalLayout.addWidget(self.denListWidget)
        self.gridLayout.addLayout(self.horizontalLayout,1,0,1,2)
        spacerItem = QtWidgets.QSpacerItem(168,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem,2,0,1,1)
        self.buttonBox = QtWidgets.QDialogButtonBox(self.widget)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox,2,1,1,1)

        self.retranslateUi(choice2DialogForm)
        self.buttonBox.rejected.connect(choice2DialogForm.reject)
        self.buttonBox.accepted.connect(choice2DialogForm.accept)
#        QtCore.QMetaObject.connectSlotsByName(choice2DialogForm)

    def retranslateUi(self, choice2DialogForm):
        choice2DialogForm.setWindowTitle(QtWidgets.QApplication.translate("choice2DialogForm", "Pick numerators and denominators", None))
        self.numLabel.setText(QtWidgets.QApplication.translate("choice2DialogForm", "Numerator", None))
        self.denLabel.setText(QtWidgets.QApplication.translate("choice2DialogForm", "Denominator", None))

