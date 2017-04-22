from setuptools import setup

setup(
  
  name='sabl_mpl',

  packages=['sabl_mpl'],

  author='Andrew Jewett',

  author_email='jewett.aij@gmail.com',

  url='https://jensenlab.caltech.edu',

  # BSD 3-Clause License:
  # - http://choosealicense.com/licenses/bsd-3-clause
  # - http://opensource.org/licenses/BSD-3-Clause

  license='BSD',

  classifiers=['Development Status :: 3 Alpha',
               'License :: OSI Approved :: BSD License',
               'Environment :: Console',
               'Operating System :: MacOS :: MacOS X',
               'Operating System :: POSIX :: Linux',
               'Operating System :: Microsoft :: Windows'
  ],

  entry_points={
      'console_scripts': ['sabl.py=sabl_mpl.sabl:main']
      },
  zip_safe=True,
  include_package_data=True
)
