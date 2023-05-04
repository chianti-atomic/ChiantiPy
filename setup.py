from setuptools import setup, find_packages



with open("README.rst", "r", encoding='utf-8') as fh:
    long_description = fh.read()


setup(name = 'ChiantiPy',
    description = 'a Python interface to the CHIANTI atomic database for astrophysical spectroscopy',
    long_description = long_description,
#    long_description_content_type = text/rst,
    version = '0.15.1',
    author = 'Ken Dere',
    author_email = 'kdere@gmu.edu',
    url = 'https://github.com/chianti-atomic/ChiantiPy',
    packages = find_packages(),
    classifiers = [
    'Development Status :: 3 - Alpha',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Intended Audience :: End Users/Desktop',
    'License :: OSI Approved :: ISC License (ISCL)',
    'Operating System :: POSIX :: Linux',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Astronomy',
    'Topic :: Scientific/Engineering :: Physics'
    ]
    )
