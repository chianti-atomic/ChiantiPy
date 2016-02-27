from distutils.core import setup
setup(name = 'ChiantiPy',
    description = 'a Python interface to the CHIANTI atomic database for astrophysical spectroscopy',
    long_description = open('README').read(),
    version = '0.7.0',
    author = 'Ken Dere',
    author_email = 'kdere@gmu.edu',
    url = 'http://chiantipy.sourceforge.net',
    download_url = 'http://sourceforge.net/projects/chiantipy',
#    package_dir = {'chianti':''},
    packages = ['chianti','chianti.core','chianti.fortranformat','chianti.Gui','chianti.Gui.gui_qt','chianti.Gui.gui_wx','chianti.Gui.gui_cl'],
#    py_modules = ['FortranFormat', 'chianti.constants', 'data', 'filters', 'ionrec', 'mputil', 'chianti.sources', 'util', 'version'],
    classifiers = [
    'Development Status :: 3 - Alpha',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Intended Audience :: End Users/Desktop',
    'License :: OSI Approved :: ISC License (ISCL)',
    'Operating System :: POSIX :: Linux',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Astronomy',
    'Topic :: Scientific/Engineering :: Physics'
    ]
    )
