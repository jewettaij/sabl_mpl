#!/usr/bin/env python
# Author: Andrew Jewett (jewett.aij at g mail)
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2020, Scripps Research
# All rights reserved.

"""
   A crude command line utility which reads the coordinates of 3 points,
   finds the circle which passes through them, and returns the center
   and the radius of that circle.

   Usage:

     curvature3pts.py < coordinates.txt

   Details:

     This program reads a text file from the standard input in multi-column
   format, with N numbers on each line delimited by spaces. (Usually N=2 or 3.)
   This file should contain exactly 3 lines (corresponding to the 3 points
   on the circle).

   This program works in an arbitrary number of dimensions.
   (The number of dimensions is specified by the number of columns in the file.)
   
"""

from math import *
import sys
import numpy as np


def CircleFrom3Points2D(r1, r2, r3):
    """
    3 points pass through a circle.  Find the center of that circle
    and its radius.  3 eqns (below) with 3 unknowns (x0, y0, r)
       (x1 - x0)^2 + (y1 - y0)^2  =  r^2
       (x2 - x0)^2 + (y2 - y0)^2  =  r^2
       (x3 - x0)^2 + (y3 - y0)^2  =  r^2
    Solve for (x0, y0) using A * (x0, y0) = B where:
    """
    B = np.array([r2[0]**2 - r1[0]**2 + 
                  r2[1]**2 - r1[1]**2,
                  r3[0]**2 - r2[0]**2 + 
                  r3[1]**2 - r2[1]**2])
    A = np.array([[2.0 * (r2[0] - r1[0]),
                   2.0 * (r2[1] - r1[1])],
                  [2.0 * (r3[0] - r2[0]),
                   2.0 * (r3[1] - r2[1])]])
    x0, y0 = np.linalg.solve(A,B)
    r = sqrt((r1[0] - x0)**2 + (r1[1] - y0)**2)
    return r, x0, y0



def CircleFrom3Points(r1, r2, r3):
    """ 
    This is the N-dimensional generalization of CircleFrom3Points2D().
    (It works in 3D and also higher dimensions.
     This function is not necessary for "sabl.py" which is a 2D program.
     Consequently, I never got around to testing it carefully. 
     Hopefully this function works, but test it first.  -A 2020-6-10)
    """

    # Decompose this into a 2D problem using Graham-Schmidt decomposition
    # of the original vectors into the basis defined by va=r1-r2 and vb=r3-r2.
    # Then apply CircleFrom3Points2D() to find the radius of curvature
    # and the central point.

    va = r1-r2
    vb = r3-r2

    ea = va / np.linalg.norm(va)
    eb = vb - np.inner(vb,ea)*ea
    eb /= np.linalg.norm(eb)

    # Now express the vectors r1-r2, r2-r2, and r3-r2
    # in the basis formed by unit vectors ea and eb.
    # The resutling _r1, _r2, _r3 vectors are 2D vectors.
    
    _r1 = np.array([np.inner(va, ea), np.inner(va, eb)])
    _r2 = np.array(r2-r2)  # (this should be the zero vector)
    _r3 = np.array([np.inner(vb, ea), np.inner(vb, eb)])

    # Now invoke "CircleFrom3Points2D()" to calculate the radius and center
    # of the circle in this 2D coordinate system
    r, x0, y0 =  CircleFrom3Points2D(_r1, _r2, _r3)

    # Now convert x0, y0 back into the original coordinate system
    r0 = r2 + x0*ea + y0*eb

    # Now return the results to the caller
    return r, r0



def main():
    lines = sys.stdin.readlines()

    # r1 r2 and r3 are arrays containing the xyz coordinates of the 3 points
    r1 = np.array(list(map(float, lines[0].strip().split())))  #x,y,z of point 1
    r2 = np.array(list(map(float, lines[1].strip().split())))  #x,y,z of point 2
    r3 = np.array(list(map(float, lines[2].strip().split())))  #x,y,z of point 3

    radius, center = CircleFrom3Points(r1, r2, r3)

    sys.stdout.write('circle_center =')
    for i in range(0, center.size):
        sys.stdout.write(' '+str(center[i]))
    sys.stdout.write('\n')

    sys.stdout.write('radius = '+str(radius)+'\n')

if __name__ == '__main__':
    main()
