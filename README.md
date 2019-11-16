[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![GitHub repo size](https://img.shields.io/github/repo-size/jewettaij/sabl_mpl)]()



sabl_mpl
===========

![](./doc/example_multiple_measurements_LR.png)

## WARNING: DOES NOT RUN ON MODERN PYTHON DISTRIBUTIONS

As of 2019-11-15, *sabl_mpl* still works with python v2.7 and matplotlib v1.3.1.

However this program does *not* work with more recent versions of
matplotlib, including 
[v2.0.0](https://github.com/jewettaij/sabl_mpl/issues/1),
and
[v3.0.1](https://github.com/jewettaij/sabl_mpl/issues/2),

*There is nothing I can do to fix this,
 because the bugs are in matplotlib itself.*

As a work-around, try installing an old version of python or locate an
apple Mac.  (Some modern apple laptops still ship with python 2.7 and
matplotlib 1.3.1.)

*(My apologies to anyone who uses this program.
  Over the last couple years I have watched matplotlib's GUI
  code deteriorate significantly.
  I recommend steering clear of using generic matplotlib buttons and
  widgets for any future projects that require a GUI.)


##  Usage

    sabl.py [-p w] [-discard h r] [-nspline n] [-alpha a] image_file.png 


##  Description

"sabl_mpl.py" is a simple spline drawing program which uses the matplotlib gui.  It was intended to be used to measure the curvature of cell membranes and other features in electron microscopy images.  It is a convenient way to repetitively measure many distances and angles in a 2-D image *manually*.  It requires python and matplotlib to be installed in advance.
(This software was featured in Yao et al. EMBO J. (2017), DOI 10.15252/embj.201696235)

Images showing how to use it are located in the "docs" directory, as well as detailed usage and installation instructions.

(NOTE: You must click using the MIDDLE mouse button to draw and make measurements)


##  Note:

This is a graphical program which ***must*** be invoked from a terminal.  (This is because I was too lazy to write a proper GUI.  This program will report measurements to the terminal where it was invoked.  So when invoking, I suggest using "sabl_mpl.py", *not* "sabl_mpl.py &".)


## Requirements

This program requires *python 2.6* or later.

It also requires the *matplotlib* python module.

The *cv2* python module is recommended but not required.
*(If the "cv2" module is installed, then sabl_mpl will be
able to read a wider variety of image formats.)*


## Installation Instructions

There are two ways to install sabl_mpl:


### Installation using pip

If you are familiar with pip, then run the following command from within the directory where this README file is located:

    pip install .

Make sure that your default pip install bin directory is in your PATH.  (This is usually something like ~/.local/bin/ or ~/anaconda3/bin/.  If you have installed anaconda, this will be done for you automatically.)  Later, you can uninstall sabl_mpl using:

    pip uninstall sabl_mpl


### Manual Installation method:

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

sabl_mpl is available under the terms of the [MIT License](LICENSE.md).

