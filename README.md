sabl_mpl
===========

##  Usage

    sabl.py [-p w] [-discard h r] [-nspline n] [-alpha a] image_file.png 

##  Description

"sabl_mpl.py" is a simple spline drawing program which uses the matplotlib gui.  It was intended to be used to measure the curvature of cell membranes in electron microscopy images.  It has some features to make it convenient to repetitively measure many distances and angles in any 2-D image manually.  It requires python and matplotlib to be installed in advance.
(This software was featured in Yao et al. (2017), DOI 10.15252/embj.201696235)

Images showing how to use it are located in the "docs" directory, as well as detailed usage and installation instructions.

(NOTE: You must click using the MIDDLE mouse button to draw and make measurements)

##  Note:

This is a graphical program which must be invoked from a terminal.  (When doing so, do not use "&".)

## Requirements

This program requires python 2.6 or later.
It also requires that the "matplotlib" python module is installed.  As with all linux and unix programs, this script must be copied to a directory which is located in your PATH.  (For more information see "Installation Details" at the end of this file.)

## Installation Instructions

There are two ways to install sabl_mpl:

## Installation using pip

If you are familiar with pip, then run the following command from within the directory where this README file is located:

    pip install .

Make sure that your default pip install bin directory is in your PATH.  (This is usually something like ~/.local/bin/ or ~/anaconda3/bin/.  If you have installed anaconda, this will be done for you automatically.)  Later, you can uninstall sabl_mpl using:

    pip uninstall sabl_mpl


## Manual Installation method:

Alternatively, you can edit your PATH variable manually to include
the subdirectory where the sabl.py script is located,
as well as the subdirectory where most of the python scripts are located.
Suppose the directory with this README file is named ``sabl_mpl''
and is located in your home directory:

If you use the bash shell, typically you would edit your 
`~/.profile`, `~/.bash_profile` or `~/.bashrc` files 
to contain the following lines:

    export PATH="$PATH:$HOME/sabl_mpl/sabl_mpl"

## License

sabl_mpl is available under the terms of the open-source 3-clause BSD 
license.  (See `LICENSE.md`.)

