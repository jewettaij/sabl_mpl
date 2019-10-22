from setuptools import setup

setup(
  
  name='sabl_mpl',

  packages=['sabl_mpl'],

  author='Andrew Jewett',

  author_email='jewett.aij@gmail.com',

  url='https://jensenlab.caltech.edu',

  install_requires=[
      'numpy',
      'matplotlib',
  ],
    
  license='MIT',

  classifiers=['Development Status :: 3 - Alpha',
               'License :: OSI Approved :: MIT License',
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
