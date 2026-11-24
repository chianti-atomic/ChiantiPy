'''
Command line selection dialogs.
'''
import os
import fnmatch

'command-line selection dialogs'


def chpicker(path,  filter='*.*', label=None):
    '''Select a filename from using a command line dialog.

    the label keyword is included for consistency but does nothing'''
    if not os.path.isdir(path):
        path=os.curdir
#        if pattern is None:
#            pattern="All files (*.*)"
    allNames=os.listdir(path)
    names = fnmatch.filter(allNames, filter)
    print(' - make a selection from one of these - ')
    for i, one in enumerate(names):
        print('%6i %s '%(i, one))
    raw = input(' type the index of your selection >>. ')
    fileName = os.path.join(path, names[int(raw)])
    return fileName


    #
#
class selectorDialog:
    '''Make a single or multiple selection from a list of items.

    expects the input of an array of items, will select one or more
    the label and parent keywords are for consistency with other modules but do nothing'''
    def __init__(self, items, label=None ,  parent=None, multiChoice=False):
        #
        print(' - make a selection from these - ')
        for i, one in enumerate(items):
            print('%6i %s '%( i, one))
        if multiChoice:
            print(' type the comma-separated index/indices of your selection and hit return')
        else:
            print(' type the index of your selection and hit return')
        raw = input('>>> ')
        #
        sraw = list(raw.split(','))
        self.selectedIndex = []
        self.selectedText = []
        for one in sraw:
            self.selectedIndex.append(int(one))
            self.selectedText.append(items[int(one)])
    #
