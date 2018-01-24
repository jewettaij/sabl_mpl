#!/usr/bin/env python

# sabl.py: "spline analysis of bumps and lumps..."
#
# As of v1.3.1, matplotlib is currently too unstable and feature-poor for me
# me to recommend using it for an interactive program much fancier than this.
# If you need a more sophisticated graphical user interface (with text-entry
# fields, or reliable check-buttons, for example), then search google to find
# out how to combine matplotlib with TkInter or PyQt, or avoid matplotlib
# entirely.
#
# Incidentally, the following web pages were helpful:
#http://matplotlib.org/users/event_handling.html
#http://stackoverflow.com/questions/21352580/matplotlib-plotting-numerous-disconnected-line-segments-with-different-colors
#
# -Andrew 2016-5-10


from math import *
import sys
import numpy as np
import random
from matplotlib import pyplot as plt
from matplotlib import patches as patches
from matplotlib.widgets import Button, RadioButtons, CheckButtons
import bisect
from collections import deque

g_program_name = __file__.split('/')[-1]
g_date_str     = '2018-1-23'
g_version_str  = '0.49'
g_filename_in  = ''



class InputError(Exception):
    """ A generic exception object containing a string for error reporting.
        (Raising this exception implies that the caller has provided
         a faulty input file or argument.)

    """

    def __init__(self, err_msg):
        self.err_msg = err_msg
    def __str__(self):
        return self.err_msg
    def __repr__(self):
        return str(self)



def Distance(ra, rb):
    """ The distance between two points in n-dimensional space (n=2) """

    rsqr = 0.0
    for d in range(0, len(ra)):
        rsqr += (rb[d]-ra[d])**2
    return sqrt(rsqr)


def DistanceSqr(ra, rb):
    rsqr = 0.0
    for d in range(0, len(ra)):
        rsqr += (rb[d]-ra[d])**2
    return rsqr


def NormalizeInPlace(rv):
    """ Replace the contents of a vector with a new vector of length 1,
        pointing in the same direction 

    """

    r = 0.0
    for d in range(0, len(rv)):
        r += rv[d]**2
    r = sqrt(r)
    for d in range(0, len(rv)):
        rv[d] /= r


def CalcAngle2D(ra, rb):
    """ Calculate the angle between a pair of 2-dimensional vectors """
    ra_dot_rb   = ra[0]*rb[0] + ra[1]*rb[1]
    ra_cross_rb = ra[0]*rb[1] - ra[1]*rb[0]
    return atan2(ra_cross_rb, ra_dot_rb)



def FindClosest2D(xclick, yclick, xlist, ylist):
    min_distsq = 0.0
    i_min_dist = 0
    for i in range(0, len(xlist)):
        distsq = DistanceSqr((xclick, yclick),
                             (xlist[i], ylist[i]))
        if ((distsq < min_distsq) or (i == 0)):
            min_distsq = distsq
            i_min_dist = i
    return i_min_dist


def CircleFrom3Points2D(r1, r2, r3):
    # 3 points pass through a circle.  Find the center of that circle
    # and its radius.  3 eqns (below) with 3 unknowns (x0, y0, r)
    #    (x1 - x0)^2 + (y1 - y0)^2  =  r^2
    #    (x2 - x0)^2 + (y2 - y0)^2  =  r^2
    #    (x3 - x0)^2 + (y3 - y0)^2  =  r^2
    # Solve for (x0, y0) using A * (x0, y0) = B where:
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




def BetweenInclusive(a, b, c):
    return (a <= b <= c) or (c <= b <= a)


def BetweenExclusive(a, b, c):
    return (a < b < c) or (c < b < a)


## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver stolen from:
# https://gist.github.com/ofan666/1875903

def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    Solves a (non-cyclic) system of equations of the form:

      a_i*x_{i-1} + b_i*x_i + c_i*x_{i+1} = d_i
      where a_1=0, and c_n=0

    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

    '''

    nf = len(a)     # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d))     # copy the array
    for it in xrange(1, nf):
        mc = ac[it]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]

    xc = ac
    xc[-1] = dc[-1]/bc[-1]

    for il in xrange(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    del bc, cc, dc  # delete variables from memory

    return xc




def SplineEval(t, t_interval, c3a, c3b, c1a, c1b):
    #sys.stderr.write('t = ' + str(t) + '\n')
    #sys.stderr.write('t_interval = ' + str(t_interval) + '\n')
    #sys.stderr.write('c3a = ' + str(c3a) + '\n')
    #sys.stderr.write('c3b = ' + str(c3b) + '\n')
    #sys.stderr.write('c1a = ' + str(c1a) + '\n')
    #sys.stderr.write('c1b = ' + str(c1b) + '\n')
    ta = t
    tb = t_interval - t
    return (c3a*ta*ta + c1a)*ta + (c3b*tb*tb + c1b)*tb



def SplineEvalD1(t, t_interval, c3a, c3b, c1a, c1b):
    ta = t
    tb = t_interval - t
    return c3a*ta*ta*3.0 + c1a - (c3b*tb*tb*3.0 + c1b)



def SplineEvalD2(t, t_interval, c3a, c3b, c1a, c1b):
    ta = t
    tb = t_interval - t
    return c3a*ta*6.0 + c3b*tb*6.0



def SplineInterpEval(t, c3a, c3b, c1a, c1b, t_control):
    i = bisect.bisect(t_control, t) - 1
    if i < 0: i = 0
    if i >= len(t_control)-1: i = len(t_control)-2
    return SplineEval(t-t_control[i], t_control[i+1]-t_control[i],
                      c3a[i], c3b[i], c1a[i], c1b[i])


def SplineInterpEvalD1(t, c3a, c3b, c1a, c1b, t_control):
    i = bisect.bisect(t_control, t) - 1
    if i < 0: i = 0
    if i >= len(t_control)-1: i = len(t_control)-2
    return SplineEvalD1(t-t_control[i], t_control[i+1]-t_control[i],
                        c3a[i], c3b[i], c1a[i], c1b[i])



def SplineInterpEvalD2(t, c3a, c3b, c1a, c1b, t_control):
    i = bisect.bisect(t_control, t) - 1
    if i < 0: i = 0
    if i >= len(t_control)-1: i = len(t_control)-2
    return SplineEvalD2(t-t_control[i], t_control[i+1]-t_control[i],
                        c3a[i], c3b[i], c1a[i], c1b[i])



def SplineCurvature2D(t, t_interval, c3a, c3b, c1a, c1b):
    # first derivatives
    x1,y1 = SplineEvalD1(t, t_interval, c3a, c3b, c1a, c1b)
    # second derivatives
    x2,y2 = SplineEvalD2(t, t_interval, c3a, c3b, c1a, c1b)
    #print('t',t)
    #print('t_control',t_control)
    #print('i',i)
    #print('x\',y\'', (x1,y1))
    #print('x\'\',y\'\'', (x2,y2))
    denom = pow((x1*x1 + y1*y1), 1.5)
    curvature = (x1*y2 - x2*y1) / denom
    return curvature



def SplineInterpCurvature2D(t, c3a, c3b, c1a, c1b, t_control):
    i = bisect.bisect(t_control, t) - 1
    if i < 0: i = 0
    if i >= len(t_control)-1: i = len(t_control)-2
    return SplineCurvature2D(t-t_control[i], t_control[i+1]-t_control[i],
                             c3a[i], c3b[i], c1a[i], c1b[i])





def SplineInterpEvalMany(t_many, c3a, c3b, c1a, c1b, t_control):
    N = len(t_many)
    D = len(c3a[0])
    x_id = np.zeros((N,D))
    for i in range(0, N):
        # Note: Each evaluation in this loop requires O(ln(N)) time due to use
        #       of "bisect" (see the implementation of SplineEval()).
        #       A more efficient implementation would exploit the fact that
        #       the entries in t_many are robably in sorted order.  If the times
        #       in t_many are evenly distributed between each intervals between 
        #       control points, then it's more efficient to use 
        #       SplineInterpEvalNtimesPerPoint() instead.
        x_id[i] = SplineInterpEval(t_many[i], c3a, c3b, c1a, c1b, t_control)
    return x_id


def __SplineMany__(NSUB, t_interval, c3a, c3b, c1a, c1b):
    assert(len(c3a) == len(c3b) == len(c1a) == len(c1b))
    D = len(c3a)
    r = np.zeros((NSUB, D))   # r is a list of vectors of dimensionality D (D=2)
    t = np.zeros(NSUB)

    for i in range(0, NSUB):
        t[i] = i * (t_interval/NSUB)
        r[i] = SplineEval(t[i], t_interval, c3a, c3b, c1a, c1b)
        #Note: r[i] is a vector (numpy array), so this is completely equivalent:
        #for d in range(0, D):
        #    r[i][d] = SplineEval(t, t_i, t_ip1, c3a[d], c3b[d], c1a[d], c1b[d])

    return r, t


def SplineInterpEvalNtimesPerPoint(NSUB, c3a, c3b, c1a, c1b, t_control):
    assert(len(c3a) == len(c3b) == len(c1a) == len(c1b) == len(t_control))
    Ncontrol = len(c3a)
    assert(len(c3a[0]) == len(c3b[0]) == len(c1a[0]) == len(c1b[0]))
    D = len(c3a[0])
    r = np.zeros((1+(Ncontrol-1)*NSUB, D))
    t = np.zeros(1+(Ncontrol-1)*NSUB)

    for i in range(0, Ncontrol-1):
        tmp = __SplineMany__(NSUB, t_control[i+1]-t_control[i],
                             c3a[i], c3b[i], c1a[i], c1b[i])
        #r[i*NSUB:(i+1)*NSUB] = tmp[0]
        for j in range(0, NSUB):
            r[i*NSUB + j] = tmp[0][j]
            t[i*NSUB + j] = t_control[i] + tmp[1][j]

    r[-1] = SplineEval(t_control[Ncontrol-1]-t_control[Ncontrol-2], 
                       t_control[Ncontrol-1]-t_control[Ncontrol-2],
                       c3a[Ncontrol-2], c3b[Ncontrol-2], 
                       c1a[Ncontrol-2], c1b[Ncontrol-2])
    t[-1] = t_control[-1]
    return r, t






def CalcNaturalCubicSplineCoeffs(r, alpha_exponent):
    N = len(r)
    assert(N >= 2)
    D = len(r[0])
    t_control = np.zeros(N)

    # Some of the variables used in this program ("h_i", "b_i", "e_i", "v_i")
    # were used by N. J. Salamon in his class notes located here:
    # http://www.esm.psu.edu/courses/emch407/njs
    # http://www.esm.psu.edu/courses/emch407/njs/notes02/ch4_3.doc

    # I define a few additional variables:
    #
    # e_d2rdt2[d][i] is the second derivative of the spline at the ith 
    #               control point (and in the dth direction)
    #
    # Once we have figured out e_d2rdt2[d][i], we can calculate the spline
    # anywhere using:
    #   SplineEval(t, h_i, 
    #              e_{i+1} / (6*h_i), 
    #              e_i / (6*h_{i+1}),
    #              r_{i+1}/h_i - e_{i+1}*h_i,
    #              r_i/h_i - e_i*h_i)

    e_d2rdt2 = np.zeros((D,N))

    # We want to solve this system of equations for e_1, e_2, ..., e_{n-2}:
    #    (note: e_i is shorthand for e_d2rdt2[d][i])
    #
    # h_{i-1}*e_{i-1} + u_i*e_i + h_{i+1} * e_{i+1}  =  v_i
    #   where h_i, u_i and v_i are shorthand for:
    #
    # h_i = change in "time" parameter for the ith interval, which is usually:
    #     = |r_{i+1} - r_i |^alpha    (alpha varies depending on settings)
    #   and
    # u_i = 2 * (h_{i-1} + h_i )
    # v_i = 6 * (b_i     - b_{i-1})
    # b_i =  (r_i - r_{i-1}) / h_i
    #
    # ...subject to the constraints that the curve is straight at the endpoints:
    # e_0 = 0       <-- first control point  (indexing begins at 0)
    # e_{n-1} = 0   <-- this is the last control point  (indexing begins at 0)

    h_dt = np.zeros(N-1)   # h_dt[i] is the i'th time interval 
                           # in the parameterization
    b_drdt = np.zeros((N-1,D))  # b_drdt is a discrete version of the derivative

    t_total=0.0
    for i in range(0, N-1):
        t_control[i] = t_total
        sqr_distance_i_ip1 = 0.0
        for d in range(0, D):
            sqr_distance_i_ip1 += (r[i+1][d] - r[i][d])**2
        h_dt[i] = sqr_distance_i_ip1**(0.5*alpha_exponent)
        for d in range(0, D):
            b_drdt[i][d] = (r[i+1][d] - r[i][d]) / h_dt[i]
        t_total += h_dt[i]
    t_control[N-1] = t_total

    a_coeff = np.zeros(N)
    b_coeff = np.zeros(N)
    c_coeff = np.zeros(N)
    d_coeff = np.zeros((D,N))

    for i in range(1, N-1):
        # h_dt[i] is the difference in "time" in the parametric curve 
        # between pairs of control points.  If "alpha_exponent" is 0
        # then the time interval between control points is uniform.  (Typically,
        # "alpha_exponent" is 0.5.)
        a_coeff[i]    =      h_dt[i-1]
        b_coeff[i]    = 2.0*(h_dt[i-1] + h_dt[i])
        c_coeff[i]    =      h_dt[i]
        for d in range(0, D):
            d_coeff[d][i] = 6.0*(b_drdt[i][d] - b_drdt[i-1][d])

    a_coeff[0] = 0.0
    b_coeff[0] = 1.0
    c_coeff[0] = 0.0
    for d in range(0, D):
        d_coeff[d][0] = 0.0
    a_coeff[N-1] = 0.0
    b_coeff[N-1] = 1.0
    c_coeff[N-1] = 0.0
    for d in range(0, D):
        d_coeff[d][N-1] = 0.0

    if N >= 3:
        for d in range(0, D):
            e_d2rdt2[d] = TDMAsolver(a_coeff, b_coeff, c_coeff, d_coeff[d])

    # alternately, if that fails, try the matrix inverter that comes with numpy:
    #M = np.zeros((N,N))
    #for i in range(0,N):
    #    if i-1>=0:
    #        M[i][i-1] = a_coeff[i]
    #    M[i][i]   = b_coeff[i]
    #    if i+1 < N:
    #        M[i][i+1] = c_coeff[i]
    #print('M', M)
    #for d in range(0, D):
    #    e_d2rdt2[d] = np.linalg.solve(M, d_coeff[d])

    e_d2rdt2 = e_d2rdt2.transpose()
    c3a = np.zeros((N, D))
    c3b = np.zeros((N, D))
    c1a = np.zeros((N, D))
    c1b = np.zeros((N, D))
    # faster to precompute these coefficients in advance:
    for i in range(0, N-1):
        c3a[i] = e_d2rdt2[i+1] / (6*h_dt[i])
        c3b[i] = e_d2rdt2[i]   / (6*h_dt[i])
        c1a[i] = r[i+1]/h_dt[i]  -  e_d2rdt2[i+1]*h_dt[i]/6.0
        c1b[i] = r[i]/h_dt[i]    -  e_d2rdt2[i]*h_dt[i]/6.0
    # c3a = e_{i+1} / (6*h_i)
    # c3b =     e_i / (6*h_{i+1})
    # c1a = r_{i+1}/h_i - e_{i+1}*h_i
    # c1b =  r_i / h_i - e_i*h_i

    # Return these spline coefficients to the caller.
    # We can use these to quickly evaluate the spline repetatively later on

    return c3a, c3b, c1a, c1b, t_control



def FindInflectionPoints(NSUB, c3a, c3b, c1a, c1b, t_control):
    # Here I divide each interval between control points in the spline
    # into NSUB smaller intervals.  I return the locations in this list 
    # where the (magnitude of) the spline curvature has a local maxima.
    # A more advanced strategy would be to solve for these maxima algebraically
    # (or numerically using binary subdivision) to give an exact answer.
    # I can do this, if we need to.
    assert(len(c3a) == len(c3b) == len(c1a) == len(c1b) == len(t_control))
    Ncontrol = len(c3a)

    t_inflections = []
    kcurve_prev = 0.0
    t_prev = 0.0
    for i in range(0, Ncontrol-1):
        delta_t = (t_control[i+1] - t_control[i]) / NSUB
        for j in range(1, NSUB+1):
            t = t_control[i] + j*delta_t
            kcurve = SplineCurvature2D(t - t_control[i],
                                       t_control[i+1] - t_control[i], 
                                       c3a[i], c3b[i], c1a[i], c1b[i])
            if kcurve * kcurve_prev < 0:
                t_inflections.append(0.5*(t + t_prev))
            t_prev = t
            kcurve_prev = kcurve
    return t_inflections



def FindCurvatureMaxima(NSUB, c3a, c3b, c1a, c1b, t_control):

    # Here I divide each interval between control points in the spline
    # into NSUB smaller intervals.  I return the locations
    # in this list where the spline curvature changes sign.
    # A more advanced strategy would be to solve for these maxima algebraically
    # (or numerically using binary subdivision) to give an exact answer.
    # I can do this, if we need to.
    assert(len(c3a) == len(c3b) == len(c1a) == len(c1b) == len(t_control))
    Ncontrol = len(c3a)

    t_maxima = []
    k_max_curvatures = []
    kcurve_prev2 = 0.0
    kcurve_prev1 = 0.0
    t_prev2 = 0.0
    t_prev1 = 0.0

    for i in range(0, Ncontrol-1):
        delta_t = (t_control[i+1] - t_control[i]) / NSUB
        for j in range(1, NSUB+1):
            t = t_control[i] + j*delta_t
            kcurve = SplineCurvature2D(t - t_control[i],
                                       t_control[i+1] - t_control[i], 
                                       c3a[i], c3b[i], c1a[i], c1b[i])
            abs_kcurve = abs(kcurve)
            #sys.stderr.write('t='+str(t)+', kcurve='+str(kcurve)+'\n')
            if (kcurve_prev2 < kcurve_prev1) and (abs_kcurve < kcurve_prev1):
                t_maxima.append(t_prev1)
                k_max_curvatures.append(kcurve)

            t_prev2 = t_prev1
            t_prev1 = t
            kcurve_prev2 = kcurve_prev1
            kcurve_prev1 = abs_kcurve

    curve_maxima_directions = np.array([[0.0, 0.0] 
                                        for i in range(0, len(t_maxima))])

    for i in range(0, len(t_maxima)):
        t = t_maxima[i]
        tangent_vect = SplineInterpEvalD1(t,
                                          c3a, c3b, 
                                          c1a, c1b, 
                                          t_control)
        NormalizeInPlace(tangent_vect)

        # normal_vect is the tangent_vect rotated by 90 degrees
        if k_max_curvatures[i] >= 0.0:
            normal_vect = np.array([tangent_vect[1], -tangent_vect[0]])
        else:
            normal_vect = np.array([-tangent_vect[1], tangent_vect[0]])

        curve_maxima_directions[i] = normal_vect

    #print(k_max_curvatures)
    #print(curve_maxima_directions)
    return t_maxima, k_max_curvatures, curve_maxima_directions



def IdentifyShallowIntervals(min_bump_height,
                             bump_width,
                             k_max_curvatures,
                             curve_normal_directions,
                             x_interp,
                             y_interp,
                             t_interp,
                             x_candidates,
                             y_candidates,
                             t_candidates):

    assert(len(x_interp) == len(y_interp) == len(t_interp))
    N = len(x_interp)

    contour_length = 0.0
    contour_length_vs_i = np.array([0.0 for i in range(0, N)])
    for i in range(1, N):
        segment_length = sqrt((x_interp[i]-x_interp[i-1])**2 +
                              (y_interp[i]-y_interp[i-1])**2)
        contour_length += segment_length
        contour_length_vs_i[i] = contour_length

    is_shallow = [True for i_maxima in range(0, len(t_candidates))]
        

    for i_a in range(0, N-1):
        x_a = x_interp[i_a]
        y_a = y_interp[i_a]
        t_a = t_interp[i_a]

        i_b = i_a + 1
        while ((i_b < N) and 
               (contour_length_vs_i[i_b]-contour_length_vs_i[i_a]<=bump_width)):

            x_b = x_interp[i_b]
            y_b = y_interp[i_b]
            t_b = t_interp[i_b]

            xab_dir = (x_b - x_a)
            yab_dir = (y_b - y_a)
            length  = sqrt(xab_dir*xab_dir + yab_dir*yab_dir)
            xab_dir /= length
            yab_dir /= length
            if i_b == 99:
                pass

            for i in range(i_a, i_b):
                x = x_interp[i]
                y = y_interp[i]
                t = t_interp[i]

                x_minus_xa = x - x_a
                y_minus_ya = y - y_a
                distance_to_line = abs(xab_dir * y_minus_ya -
                                       yab_dir * x_minus_xa)

                #print(i_a, i_b, y_a, y_b, y_minus_ya, y - y_b, distance_to_line)
                if (distance_to_line >= min_bump_height):
                    # Then find all the points of interest 
                    # that lie within this interval

                    for i_candidate in range(0, len(t_candidates)):
                        t_candidate = t_candidates[i_candidate]
                        if (t_a <= t_candidate) and (t_candidate <= t_b):

                        #if (x_minus_xa*curve_normal_directions[i_candidate][0] +
                        #    y_minus_ya*curve_normal_directions[i_candidate][1]
                        #    > 0.0):
 
                           #print(y_minus_ya, curve_normal_directions[i_candidate])
                            is_shallow[i_candidate] = False

            i_b += 1
        #sys.stderr.write("i_a = "+str(i_a)+"\n")
    return is_shallow




def BumpHeightExceedsThreshold(t,
                               min_bump_height,
                               bump_width,
                               curve_normal_direction,
                               x_interp,
                               y_interp,
                               t_interp):
    assert(len(x_interp) == len(y_interp) ==len(t_interp))
    N = len(t_interp)
    bump_width_half  = bump_width * 0.5
    #print("bump_width_half, min_bump_height = ",
    #       bump_width_half, min_bump_height)
    #criteria_satisfied_a = False
    i = bisect.bisect(t_interp, t) - 1
    i_a = i
    contour_length = 0.0
    x = x0 = x_prev = x_interp[i]
    y = y0 = y_prev = y_interp[i]
    while (contour_length < bump_width_half) and (i_a > 0):
        i_a -= 1
        x = x_interp[i_a]
        y = y_interp[i_a]
        contour_length += sqrt((x-x_prev)*(x-x_prev) + (y-y_prev)*(y-y_prev))
        #if ((x0-x)*curve_normal_direction[0] + (y0-y)*curve_normal_direction[1]
        #    >= min_bump_height):
        #    criteria_satisfied_a = True
        #    break
        x_prev = x
        y_prev = y
    x_a = x
    y_a = y

    #criteria_satisfied_b = False
    i_b = i
    contour_length = 0.0
    x = x_prev = x_interp[i]
    y = y_prev = y_interp[i]
    while (contour_length < bump_width_half) and (i_b+1 < N):
        i_b += 1
        x = x_interp[i_b]
        y = y_interp[i_b]
        contour_length += sqrt((x-x_prev)*(x-x_prev) + (y-y_prev)*(y-y_prev))
        #if ((x0-x)*curve_normal_direction[0] + (y0-y)*curve_normal_direction[1]
        #    >= min_bump_height):
        #    criteria_satisfied_b = True
        #    break
        x_prev = x
        y_prev = y
    x_b = x
    y_b = y

    #if (criteria_satisfied_a and criteria_satisfied_b):
    #    return True

    if i_a == i_b:
        assert(bump_width == 0)
        return True
    
    I_a = i_a
    while I_a < i:
        X_a = x_interp[I_a]
        Y_a = y_interp[I_a]

        I_b = i_b
        while I_b > i:
            X_b = x_interp[I_b]
            Y_b = y_interp[I_b]

            x_direction = (X_b - X_a)
            y_direction = (Y_b - Y_a)
            direction_length = sqrt(x_direction*x_direction + y_direction*y_direction)
            x_direction /= direction_length
            y_direction /= direction_length

            x_displace = x0 - X_a
            y_displace = y0 - Y_a
            distance_to_line = abs(x_direction*y_displace -
                                   y_direction*x_displace)
            #print(Y_a, Y_b, y0, distance_to_line)
            #print(curve_normal_direction[0], curve_normal_direction[1])
            if (distance_to_line >= min_bump_height):
                distance_along_line = (x_displace*x_direction + y_displace*y_direction)
                nearest_point_on_line = (X_a + distance_along_line*x_direction,
                                         Y_a + distance_along_line*y_direction)
                displacement_to_line = (x0 - nearest_point_on_line[0],
                                        y0 - nearest_point_on_line[0])

                if (displacement_to_line[0] * curve_normal_direction[0] +
                    displacement_to_line[0] * curve_normal_direction[1] > 0.0):

                    #FOLLOWING LINES ARE FOR DEBUGGING (a specific dataset)
                    if abs(y0 - 4470.07)<0.1:
                        print((X_a, Y_a),(X_b,Y_b), distance_to_line)
                        print(curve_normal_direction[0], curve_normal_direction[1])
                        print(nearest_point_on_line[0],nearest_point_on_line[1])
                        print(sqrt(displacement_to_line[0]**2 +
                                displacement_to_line[1]**2))

                    return True

            I_b -= 1
        I_a += 1

    return False




def IdentifyShallowPoints(min_bump_height,
                          bump_width,
                          k_max_curvatures,
                          curve_normal_directions,
                          x_interp,
                          y_interp,
                          t_interp,
                          x_candidates,
                          y_candidates,
                          t_candidates):

    is_shallow = [True for i in range(0, len(t_candidates))]

    for i in range(0, len(t_candidates)):
        if BumpHeightExceedsThreshold(t_candidates[i],
                                      min_bump_height,
                                      bump_width,
                                      curve_normal_directions[i],
                                      x_interp,
                                      y_interp,
                                      t_interp):
            is_shallow[i] = False

    return is_shallow







class LineSegmentBuilder:

    def __init__(self, ax,
                 set_x_clickable = [],
                 set_y_clickable = []):
        self.x_clickable = set_x_clickable
        self.y_clickable = set_y_clickable
        self.num_markers = 3
        self.marker = '.' # default marker

        self.object_crds_x = [] # x,y coordinates of the most recently
        self.object_crds_y = [] # selected point (not yet part of a line)
        self.objects_crds_x = [] # a list of coords for pairs of points
        self.objects_crds_y = [] # that already belong to lines
        self.distances = [] # A list of distances between these pairs of points
        self.color=(0.0, 0.585, 1.0)
        self.linewidth=2.0

        # The remaining data members are specific to matplotlib:
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.object_ax_lines = []
        #self.connect()

        self.mouse_cid       = None
        self.key_press_cid   = None
        self.key_release_cid = None


    def connect(self):
        self.mouse_cid    = self.canvas.mpl_connect('button_press_event',
                                                    self.handle_mouse_press)
        self.key_press_cid = self.canvas.mpl_connect('key_press_event',
                                                     self.handle_key_press)
        self.key_release_cid = self.canvas.mpl_connect('key_release_event',
                                                       self.handle_key_release)

    def handle_mouse_press(self, event):
        if event.inaxes!=self.ax:
            return
        if event.button != 2:
            return

        if len(self.x_clickable) == 0:
            x = event.xdata
            y = event.ydata
        else:
            i_min_dist = FindClosest2D(event.xdata, event.ydata, self.x_clickable, self.y_clickable)
            x = self.x_clickable[i_min_dist]
            y = self.y_clickable[i_min_dist]

        # Disregard whenever the user clicks twice in the same location:
        if ((len(self.object_crds_x) != 0) and
            (x == self.object_crds_x[-1]) and
            (y == self.object_crds_y[-1])):
            return

        self.object_crds_x.append(x)
        self.object_crds_y.append(y)

        if len(self.object_crds_x) == 1:
            self.current_object, = self.ax.plot(self.object_crds_x,
                                                self.object_crds_y,
                                                marker=self.marker,
                                                color=self.color,
                                                linewidth=self.linewidth)
            self.current_object.set_zorder(50)
            self.object_ax_lines.append(self.current_object)
            self.canvas.draw()
        if len(self.object_crds_x) == 2:
            assert(len(self.object_crds_y) == 2)
            self.objects_crds_x.append(self.object_crds_x)
            self.objects_crds_y.append(self.object_crds_y)
            if self.num_markers >= 2:
                x0 = self.object_crds_x[0]
                y0 = self.object_crds_y[0]
                xincr = ((self.object_crds_x[1] - self.object_crds_x[0]) / 
                         (self.num_markers - 1))
                yincr = ((self.object_crds_y[1] - self.object_crds_y[0]) / 
                         (self.num_markers - 1))
                self.current_object.set_data([x0 + i*xincr for i in range(0, self.num_markers)],
                                             [y0 + i*yincr for i in range(0, self.num_markers)])
                self.current_object.set_marker(self.marker)
            else:
                self.current_object.set_data(self.object_crds_x, 
                                             self.object_crds_y)
                self.current_object.set_marker(None)
            r,theta = self.ReportDistanceAndAngle(self.object_crds_x, 
                                                  self.object_crds_y)
            self.distances.append(r)
            self.object_crds_x = []
            self.object_crds_y = []
            self.canvas.draw()


    def handle_key_press(self, event):
        if event.key == 'delete':
            self.DeleteLastPoint()

        if event.key == 'enter':
            #Print the average distance between points already selected (if any)
            self.ReportStats()
            # Clear all the arrays containing any existing coordinate data
            del self.object_crds_x [:] # x,y coordinates of the most recently
            del self.object_crds_y [:] # selected point (not yet part of a line)
            del self.objects_crds_x [:] # a list of coords for pairs of points
            del self.objects_crds_y [:] # already belonging to lines
            # Clear the array of distances (between these pairs of points)
            del self.distances [:]
            # Get rid of all the lines on the screen (between these points)
            i = len(self.object_ax_lines)-1
            while i >= 0:
                #j = len(self.ax.lines)-1
                #while (j>=0) and (self.ax.lines[j]!=self.object_ax_lines[i]):
                #    j -= 1
                #assert(j >= 0):
                j = self.ax.lines.index(self.object_ax_lines[i])
                del self.ax.lines[j]
                del self.object_ax_lines[i]
                i -= 1
            self.canvas.draw()

        elif event.key in ('d','D'):
            sys.stderr.write("Divide each line segment into how many intervals?: ")
            self.num_markers = int(sys.stdin.readline())+1
            if self.num_markers > 1:
                self.marker='o'
            else:
                self.marker='.'
            sys.stderr.write("(This only effects the line segments that you draw in the future.)\n")


    def handle_key_release(self, event):
        # do nothing for now
        pass


    def DeleteLastPoint(self):
        assert(len(self.object_crds_x) == len(self.object_crds_y))
        if len(self.object_crds_x) == 1:
            del self.object_crds_x[-1]
            del self.object_crds_y[-1]
            assert(len(self.object_ax_lines) > 0)
            #i = 0
            #while ((i < len(self.ax.lines)) and 
            #       (self.ax.lines[i] != self.object_ax_lines[-1])
            #    i += 1
            #if i < len(self.ax.lines):
            #    del self.ax.lines[i]
            assert(self.object_ax_lines[-1] in self.ax.lines)
            i = self.ax.lines.index(self.object_ax_lines[-1])
            del self.ax.lines[i]
            del self.object_ax_lines[-1]
            if len(self.object_ax_lines) > 0:
                self.current_object = self.object_ax_lines[-1]
        #elif len(self.object_crds_x) == 2:
        #    del self.object_crds_x[-1]
        #    del self.object_crds_y[-1]
        #    self.current_object.set_data(self.object_crds_x, 
        #                                 self.object_crds_y)
        else:
            if len(self.distances) > 0:
                assert(len(self.distances) == 
                       len(self.objects_crds_x) == 
                       len(self.objects_crds_y))
                del self.distances[-1]
                sys.stderr.write("  deleted previous distance measurement\n")
                self.object_crds_x = self.objects_crds_x[-1][0:1]
                self.object_crds_y = self.objects_crds_y[-1][0:1]
                del self.objects_crds_x[-1]
                del self.objects_crds_y[-1]
                self.current_object.set_data(self.object_crds_x, 
                                             self.object_crds_y)
                self.current_object.set_marker('.')
        self.canvas.draw()


    def ReportDistance(self, x_control, y_control):
        r = sqrt((x_control[1]-x_control[0])**2 + (y_control[1]-y_control[0])**2)
        #sys.stdout.write('distance = '+str(r)+'\n')
        sys.stdout.write('distance = '+str(r))
        return r


    def ReportDistanceAndAngle(self, x_control, y_control):
        r = sqrt((x_control[1]-x_control[0])**2 + (y_control[1]-y_control[0])**2)
        theta = atan2(y_control[1] - y_control[0], x_control[1] - x_control[0])
        theta *= 180.0 / pi
        sys.stdout.write('distance = '+str(r)+
                         '  ( angle_from_x = '+str(theta)+' )\n')
        return r,theta


    def ReportStats(self):
        if len(self.distances) == 0:
            return
        r_ave = 0.0
        for r in self.distances:
            r_ave += r
        r_ave /= len(self.distances)
        r_stddev = 0.0
        sys.stdout.write('average')
        if len(self.distances) > 1:
            sys.stdout.write('  stddev   num_measurements')
            for r in self.distances:
                r_stddev += (r - r_ave)*(r - r_ave)
            r_stddev = sqrt(r_stddev / (len(self.distances)-1))
            sys.stdout.write('\n'+str(r_ave)+' '+str(r_stddev)+' '+
                             str(len(self.distances))+'\n')
        else:
            sys.stdout.write('\n'+str(r_ave)+'\n')


    def disconnect(self):
        if self.mouse_cid:
            self.canvas.mpl_disconnect(self.mouse_cid)
        if self.key_press_cid:
            self.canvas.mpl_disconnect(self.key_press_cid)
        if self.key_release_cid:
            self.canvas.mpl_disconnect(self.key_release_cid)


    def Finalize(self):
        if len(self.object_crds_x) > 0:
            del self.object_crds_x[:]
            del self.object_crds_y[:]
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
        self.ReportStats()
        self.disconnect()


    def __del__(self):
        self.Finalize()









class TriangleBuilder:

    def __init__(self, ax,
                 set_x_clickable = [],
                 set_y_clickable = []):

        self.x_clickable = set_x_clickable
        self.y_clickable = set_y_clickable

        self.object_crds_x = [] # x,y coordinates of the most recently
        self.object_crds_y = [] # selected point (not yet part of a line)
        self.objects_crds_x = [] # a list of coords for pairs of points
        self.objects_crds_y = [] # that already belong to lines
        self.angles = [] # A list of angles between these triplets of points

        # The remaining data members are specific to matplotlib:
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.object_ax_lines = []
        self.color = (0.0, 1.0, 0.0)
        self.linewidth = 2.0
        #self.connect()

        self.mouse_cid    = None
        self.key_press_cid = None
        self.key_release_cid = None


    def connect(self):
        self.mouse_cid    = self.canvas.mpl_connect('button_press_event',
                                                    self.handle_mouse_press)
        self.key_press_cid = self.canvas.mpl_connect('key_press_event',
                                                     self.handle_key_press)
        self.key_release_cid = self.canvas.mpl_connect('key_release_event',
                                                       self.handle_key_release)

    def handle_mouse_press(self, event):
        if event.inaxes!=self.ax:
            return
        if event.button != 2:
            return

        if len(self.x_clickable) == 0:
            x = event.xdata
            y = event.ydata
        else:
            i_min_dist = FindClosest2D(event.xdata, event.ydata, self.x_clickable, self.y_clickable)
            x = self.x_clickable[i_min_dist]
            y = self.y_clickable[i_min_dist]

        # Disregard whenever the user clicks twice in the same location:
        if ((len(self.object_crds_x) != 0) and
            (x == self.object_crds_x[-1]) and
            (y == self.object_crds_y[-1])):
            return

        self.object_crds_x.append(x)
        self.object_crds_y.append(y)
        assert(len(self.object_crds_x) == len(self.object_crds_y))
        if len(self.object_crds_x) == 1:
            self.current_object, = self.ax.plot(self.object_crds_x,
                                                self.object_crds_y,
                                                marker='.',
                                                color=self.color, 
                                                linewidth=self.linewidth)
            self.current_object.set_zorder(30)
            self.object_ax_lines.append(self.current_object)
            self.canvas.draw()

        if len(self.object_crds_x) == 2:
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
            self.canvas.draw()

        if len(self.object_crds_x) == 3:
            self.objects_crds_x.append(self.object_crds_x)
            self.objects_crds_y.append(self.object_crds_y)
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
            self.current_object.set_marker(None)
            a = self.ReportAngle(self.object_crds_x, 
                                 self.object_crds_y)
            self.angles.append(a)
            self.object_crds_x = []
            self.object_crds_y = []
            self.canvas.draw()

    def handle_key_press(self, event):
        if event.key == 'delete':
            self.DeleteLastPoint()
        if event.key == 'enter':
            #Print the average distance between points already selected (if any)
            self.ReportStats()
            # Clear all the arrays containing any existing coordinate data
            del self.object_crds_x [:] # x,y coordinates of the most recently
            del self.object_crds_y [:] # selected point (not yet part of a line)
            del self.objects_crds_x [:] # a list of coords for pairs of points
            del self.objects_crds_y [:] # already belonging to lines
            # Clear the array of distances (between these pairs of points)
            del self.angles [:]
            # Get rid of all the lines on the screen (between these points)
            i = len(self.object_ax_lines)-1
            while i >= 0:
                #j = len(self.ax.lines)-1
                #while (j>=0) and (self.ax.lines[j]!=self.object_ax_lines[i]):
                #    j -= 1
                #assert(j >= 0):
                j = self.ax.lines.index(self.object_ax_lines[i])
                del self.ax.lines[j]
                del self.object_ax_lines[i]
                i -= 1
            self.canvas.draw()

    def handle_key_release(self, event):
        # do nothing for now
        pass

    def DeleteLastPoint(self):
        assert(len(self.object_crds_x) == len(self.object_crds_y))
        if len(self.object_crds_x) == 1:
            del self.object_crds_x[-1]
            del self.object_crds_y[-1]
            assert(len(self.object_ax_lines) > 0)
            #i = 0
            #while ((i < len(self.ax.lines)) and 
            #       (self.ax.lines[i] != self.object_ax_lines[-1])
            #    i += 1
            #if i < len(self.ax.lines):
            #    del self.ax.lines[i]
            assert(self.object_ax_lines[-1] in self.ax.lines)
            i = self.ax.lines.index(self.object_ax_lines[-1])
            del self.ax.lines[i]
            del self.object_ax_lines[-1]
            if len(self.object_ax_lines) > 0:
                self.current_object = self.object_ax_lines[-1]
        elif len(self.object_crds_x) == 2:
            del self.object_crds_x[-1]
            del self.object_crds_y[-1]
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
        else:
            if len(self.angles) > 0:
                assert(len(self.angles) == 
                       len(self.objects_crds_x) == 
                       len(self.objects_crds_y))
                del self.angles[-1]
                sys.stderr.write("  deleted previous angle measurement\n")
                self.object_crds_x = self.objects_crds_x[-1][0:2]
                self.object_crds_y = self.objects_crds_y[-1][0:2]
                del self.objects_crds_x[-1]
                del self.objects_crds_y[-1]
                self.current_object.set_data(self.object_crds_x, 
                                             self.object_crds_y)
                self.current_object.set_marker('.')
        self.canvas.draw()

    def ReportAngle(self, 
                    xcrds, 
                    ycrds):
        assert(len(xcrds) == len(ycrds) == 3)
        ra = (xcrds[0]-xcrds[1], ycrds[0]-ycrds[1])
        rb = (xcrds[2]-xcrds[1], ycrds[2]-ycrds[1])
        rc = (xcrds[2]-xcrds[0], ycrds[2]-ycrds[0])
        theta012 = abs(CalcAngle2D(ra, rb))
        theta201 = abs(CalcAngle2D(ra, rc))
        #direction02 = [xcrds[2]-xcrds[0], ycrds[2]-ycrds[0]]
        #NormalizeInPlace(direction02)
        #self.tri_height.set_data([xcrds[1], xcrds[0]+direction02[0]*length01*cos(theta201)],
        #                         [ycrds[1], ycrds[0]+direction02[1]*length01*cos(theta201)])
        length01 = Distance((xcrds[0],ycrds[0]), (xcrds[1],ycrds[1]))
        tri_height = length01 * sin(theta201)
        length02 = Distance((xcrds[0],ycrds[0]), (xcrds[2],ycrds[2]))

        sys.stdout.write('angle = '+str(theta012 * 180.0 / pi)+
                         '  ( height = '+str(tri_height)+
                         '  base = '+str(length02)+' )\n')
        return theta012

    def ReportStats(self):
        if len(self.angles) == 0:
            return
        theta_ave = 0.0
        for theta in self.angles:
            theta_ave += abs(theta)
        theta_ave /= len(self.angles)
        theta_stddev = 0.0
        sys.stdout.write('average')
        if len(self.angles) > 1:
            sys.stdout.write('  stddev   num_measurements')
            for theta in self.angles:
                theta_stddev += (abs(theta)-theta_ave) * (abs(theta)-theta_ave)
            theta_stddev = sqrt(theta_stddev / (len(self.angles)-1))
            sys.stdout.write('\n'+str(theta_ave*180.0/pi)+' '+str(theta_stddev*180/pi)+' '+
                             str(len(self.angles))+'\n')
        else:
            sys.stdout.write('\n'+str(theta_ave)+'\n')

    def disconnect(self):
        if self.mouse_cid:
            self.canvas.mpl_disconnect(self.mouse_cid)
        if self.key_press_cid:
            self.canvas.mpl_disconnect(self.key_press_cid)
        if self.key_release_cid:
            self.canvas.mpl_disconnect(self.key_release_cid)

    def Finalize(self):
        if len(self.object_crds_x) > 0:
            del self.object_crds_x[:]
            del self.object_crds_y[:]
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
        self.ReportStats()
        self.disconnect()

    def __del__(self):
        self.Finalize()






class TriangleHeights:

    def __init__(self, ax,
                 set_x_clickable = [],
                 set_y_clickable = []):

        self.x_clickable = set_x_clickable
        self.y_clickable = set_y_clickable

        self.object_crds_x = [] # x,y coordinates of the most recently
        self.object_crds_y = [] # selected point (not yet part of a line)
        self.objects_crds_x = [] # a list of coords for pairs of points
        self.objects_crds_y = [] # that already belong to lines
        self.heights = [] # A list of heights between these triplets of points

        # The remaining data members are specific to matplotlib:
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.object_ax_lines = []
        self.color = (0.3, 0.5, 0.0)
        self.linewidth = 2.0
        #self.connect()

        self.mouse_cid    = None
        self.key_press_cid = None
        self.key_release_cid = None


    def connect(self):
        self.mouse_cid    = self.canvas.mpl_connect('button_press_event',
                                                    self.handle_mouse_press)
        self.key_press_cid = self.canvas.mpl_connect('key_press_event',
                                                     self.handle_key_press)
        self.key_release_cid = self.canvas.mpl_connect('key_release_event',
                                                       self.handle_key_release)

    def handle_mouse_press(self, event):
        if event.inaxes!=self.ax:
            return
        if event.button != 2:
            return

        if len(self.x_clickable) == 0:
            x = event.xdata
            y = event.ydata
        else:
            i_min_dist = FindClosest2D(event.xdata, event.ydata, self.x_clickable, self.y_clickable)
            x = self.x_clickable[i_min_dist]
            y = self.y_clickable[i_min_dist]

        # Disregard whenever the user clicks twice in the same location:
        if ((len(self.object_crds_x) != 0) and
            (x == self.object_crds_x[-1]) and
            (y == self.object_crds_y[-1])):
            return

        self.object_crds_x.append(x)
        self.object_crds_y.append(y)
        assert(len(self.object_crds_x) == len(self.object_crds_y))
        if len(self.object_crds_x) == 1:
            self.current_object, = self.ax.plot(self.object_crds_x,
                                                self.object_crds_y,
                                                marker='.',
                                                color=self.color, 
                                                linewidth=self.linewidth)
            self.current_object.set_zorder(30)
            self.object_ax_lines.append(self.current_object)
            self.canvas.draw()

        if len(self.object_crds_x) == 2:
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
            self.canvas.draw()

        if len(self.object_crds_x) == 3:
            self.objects_crds_x.append(self.object_crds_x)
            self.objects_crds_y.append(self.object_crds_y)
            self.current_object.set_marker(None)
            xcrds = self.object_crds_x
            ycrds = self.object_crds_y
            r10 = (xcrds[0]-xcrds[1], ycrds[0]-ycrds[1])
            r01 = (xcrds[1]-xcrds[0], ycrds[1]-ycrds[0])
            r12 = (xcrds[2]-xcrds[1], ycrds[2]-ycrds[1])
            r02 = (xcrds[2]-xcrds[0], ycrds[2]-ycrds[0])
            theta012 = abs(CalcAngle2D(r10, r12))
            theta201 = abs(CalcAngle2D(r01, r02))
            direction02 = [xcrds[2]-xcrds[0], ycrds[2]-ycrds[0]]
            NormalizeInPlace(direction02)
            length01 = Distance((xcrds[0],ycrds[0]), (xcrds[1],ycrds[1]))
            length02 = Distance((xcrds[0],ycrds[0]), (xcrds[2],ycrds[2]))
            self.current_object.set_data([xcrds[0], 
                                          xcrds[1], 
                                          xcrds[2], 
                                          xcrds[0], 
                                          xcrds[2], 
                                          xcrds[1], 
                                          xcrds[0]+direction02[0]*length01*cos(theta201)],
                                          [ycrds[0],
                                           ycrds[1], 
                                           ycrds[2], 
                                           ycrds[0], 
                                           ycrds[2], 
                                           ycrds[1], 
                                           ycrds[0]+direction02[1]*length01*cos(theta201)])
            self.current_object.set_marker(None)
            self.current_object.set_linewidth(self.linewidth)
            tri_height = length01 * sin(theta201)

            sys.stdout.write('height = '+str(tri_height)+
                             '  ( base = '+str(length02)+
                             '   angle = '+str(theta012 * 180.0 / pi)+' )\n')

            self.heights.append(tri_height)
            self.object_crds_x = []
            self.object_crds_y = []
            self.canvas.draw()

    def handle_key_press(self, event):
        if event.key == 'delete':
            self.DeleteLastPoint()
        if event.key == 'enter':
            #Print the average distance between points already selected (if any)
            self.ReportStats()
            # Clear all the arrays containing any existing coordinate data
            del self.object_crds_x [:] # x,y coordinates of the most recently
            del self.object_crds_y [:] # selected point (not yet part of a line)
            del self.objects_crds_x [:] # a list of coords for pairs of points
            del self.objects_crds_y [:] # already belonging to lines
            # Clear the array of distances (between these pairs of points)
            del self.heights [:]
            # Get rid of all the lines on the screen (between these points)
            i = len(self.object_ax_lines)-1
            while i >= 0:
                #j = len(self.ax.lines)-1
                #while (j>=0) and (self.ax.lines[j]!=self.object_ax_lines[i]):
                #    j -= 1
                #assert(j >= 0):
                j = self.ax.lines.index(self.object_ax_lines[i])
                del self.ax.lines[j]
                del self.object_ax_lines[i]
                i -= 1
            self.canvas.draw()

    def handle_key_release(self, event):
        # do nothing for now
        pass

    def DeleteLastPoint(self):
        assert(len(self.object_crds_x) == len(self.object_crds_y))
        if len(self.object_crds_x) == 1:
            del self.object_crds_x[-1]
            del self.object_crds_y[-1]
            assert(len(self.object_ax_lines) > 0)
            #i = 0
            #while ((i < len(self.ax.lines)) and 
            #       (self.ax.lines[i] != self.object_ax_lines[-1])
            #    i += 1
            #if i < len(self.ax.lines):
            #    del self.ax.lines[i]
            assert(self.object_ax_lines[-1] in self.ax.lines)
            i = self.ax.lines.index(self.object_ax_lines[-1])
            del self.ax.lines[i]
            del self.object_ax_lines[-1]
            if len(self.object_ax_lines) > 0:
                self.current_object = self.object_ax_lines[-1]
        elif len(self.object_crds_x) == 2:
            del self.object_crds_x[-1]
            del self.object_crds_y[-1]
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
        else:
            if len(self.heights) > 0:
                assert(len(self.heights) == 
                       len(self.objects_crds_x) == 
                       len(self.objects_crds_y))
                del self.heights[-1]
                self.object_crds_x = self.objects_crds_x[-1][0:2]
                self.object_crds_y = self.objects_crds_y[-1][0:2]
                sys.stderr.write("  deleted previous height measurement\n")
                del self.objects_crds_x[-1]
                del self.objects_crds_y[-1]
                self.current_object.set_data(self.object_crds_x, 
                                             self.object_crds_y)
                self.current_object.set_marker('o')
        self.canvas.draw()

    def ReportStats(self):
        if len(self.heights) == 0:
            return
        height_ave = 0.0
        for height in self.heights:
            height_ave += abs(height)
        height_ave /= len(self.heights)
        height_stddev = 0.0
        sys.stdout.write('average')
        if len(self.heights) > 1:
            sys.stdout.write('  stddev   num_measurements')
            for height in self.heights:
                height_stddev += (abs(height)-height_ave) * (abs(height)-height_ave)
            height_stddev = sqrt(height_stddev / (len(self.heights)-1))
            sys.stdout.write('\n'+str(height_ave)+' '+str(height_stddev)+' '+
                             str(len(self.heights))+'\n')
        else:
            sys.stdout.write('\n'+str(height_ave)+'\n')

    def disconnect(self):
        if self.mouse_cid:
            self.canvas.mpl_disconnect(self.mouse_cid)
        if self.key_press_cid:
            self.canvas.mpl_disconnect(self.key_press_cid)
        if self.key_release_cid:
            self.canvas.mpl_disconnect(self.key_release_cid)

    def Finalize(self):
        if len(self.object_crds_x) > 0:
            del self.object_crds_x[:]
            del self.object_crds_y[:]
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
        self.ReportStats()
        self.disconnect()

    def __del__(self):
        self.Finalize()







class ThreePointCircleBuilder:
    def __init__(self, ax,
                 set_x_clickable = [],
                 set_y_clickable = []):

        self.object_crds_x = [] # x,y coordinates of the most recently
        self.object_crds_y = [] # selected point (not yet part of a line)
        self.objects_crds_x = [] # a list of coords for pairs of points
        self.objects_crds_y = [] # that already belong to lines
        self.radii = []   # A list of radii of circles
        self.centers = [] # A list of x,y pairs for the center of these circles
        self.x_clickable = set_x_clickable
        self.y_clickable = set_y_clickable
        self.mouse_position_x = -1
        self.mouse_position_y = -1
        self.copy_radius_offset = 0.0
        # The remaining data members are specific to matplotlib:
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.object_ax_lines = []
        self.circles = []
        self.color = (0.6, 0.3, 0.0)
        self.linewidth = 3.0
        #self.connect()
        self.mouse_cid    = None
        self.key_press_cid = None
        self.key_release_cid = None
        self.mouse_move_cid = None


    def connect(self):
        self.mouse_cid    = self.canvas.mpl_connect('button_press_event',
                                                    self.handle_mouse_press)
        self.key_press_cid = self.canvas.mpl_connect('key_press_event',
                                                     self.handle_key_press)
        self.key_release_cid = self.canvas.mpl_connect('key_release_event',
                                                       self.handle_key_release)
        self.mouse_move_cid = self.canvas.mpl_connect('motion_notify_event',
                                                      self.handle_mouse_move)

    def handle_mouse_press(self, event):
        if event.inaxes!=self.ax:
            return
        if event.button != 2:
            return

        if len(self.x_clickable) == 0:
            x = event.xdata
            y = event.ydata
        else:
            i_min_dist = FindClosest2D(event.xdata, event.ydata, 
                                       self.x_clickable, self.y_clickable)
            x = self.x_clickable[i_min_dist]
            y = self.y_clickable[i_min_dist]

        # Disregard whenever the user clicks twice in the same location:
        if ((len(self.object_crds_x) != 0) and
            (x == self.object_crds_x[-1]) and
            (y == self.object_crds_y[-1])):
            return

        self.object_crds_x.append(x)
        self.object_crds_y.append(y)
        assert(len(self.object_crds_x) == len(self.object_crds_y))
        if len(self.object_crds_x) == 1:
            self.current_object, = self.ax.plot(self.object_crds_x,
                                                self.object_crds_y,
                                                '.',
                                                color=self.color)
                                                #linewidth=self.linewidth)
            self.current_object.set_zorder(30)
            self.object_ax_lines.append(self.current_object)
            self.canvas.draw()

        if len(self.object_crds_x) == 2:
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
            self.canvas.draw()

        if len(self.object_crds_x) == 3:
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
            self.AddCircle(self.object_crds_x, self.object_crds_y)


    def handle_key_press(self, event):
        if event.key == 'delete':
            self.DeleteLastPoint()
        if event.key == 'enter':
            #Print the average distance between points already selected (if any)
            self.ReportStats()
            # Clear all the arrays containing any existing coordinate data
            del self.object_crds_x [:] # x,y coordinates of the most recently
            del self.object_crds_y [:] # selected point (not yet part of a line)
            del self.objects_crds_x [:] # a list of coords for pairs of points
            del self.objects_crds_y [:] # already belonging to lines
            # Clear the array of distances (between these pairs of points)
            del self.radii [:]
            del self.centers [:]
            # Get rid of all the lines on the screen (between these points)
            i = len(self.object_ax_lines)-1
            while i >= 0:
                #j = len(self.ax.lines)-1
                #while (j>=0) and (self.ax.lines[j]!=self.object_ax_lines[i]):
                #    j -= 1
                #assert(j >= 0):
                j = self.ax.lines.index(self.object_ax_lines[i])
                del self.ax.lines[j]
                del self.object_ax_lines[i]
                if len(self.circles) > i:
                    self.circles[i].remove()
                    del self.circles[i]
                i -= 1
            self.canvas.draw()

        if event.key == 'c':
            # Make a copy of the current circle, with the central point
            # located at the current mouse location
            assert(len(self.objects_crds_x) == len(self.radii))
            if len(self.objects_crds_x) > 0:
                r = self.radii[-1] + self.copy_radius_offset
                copy_scale = r / self.radii[-1]
                assert(len(self.objects_crds_x[-1]) == 3)
                self.object_crds_x = [ (self.mouse_position_x +
                                        copy_scale *
                                        (self.objects_crds_x[-1][0] - 
                                         self.objects_crds_x[-1][1])),
                                       self.mouse_position_x,
                                       (self.mouse_position_x +
                                        copy_scale *
                                        (self.objects_crds_x[-1][2] - 
                                         self.objects_crds_x[-1][1])) ]
                self.object_crds_y = [ (self.mouse_position_y +
                                        copy_scale *
                                        (self.objects_crds_y[-1][0] - 
                                         self.objects_crds_y[-1][1])),
                                       self.mouse_position_y,
                                       (self.mouse_position_y +
                                        copy_scale *
                                        (self.objects_crds_y[-1][2] - 
                                         self.objects_crds_y[-1][1])) ]
            #self.current_object, = self.ax.plot(self.object_crds_x,
            #                                    self.object_crds_y,
            #                                    '.',
            #                                    color=self.color)
            self.AddCircle(self.object_crds_x, self.object_crds_y)

        if event.key == 'a':
            sys.stderr.write('What number do you want to add to the radius? (copied circles only)?\n'
                             '  Enter a number or a numeric expression: ')
            s = sys.stdin.readline()
            self.copy_radius_offset = eval(s)


    def handle_mouse_move(self, event):
        self.mouse_position_x = event.xdata
        self.mouse_position_y = event.ydata


    def handle_key_release(self, event):
        # do nothing for now
        pass


    def AddCircle(self, set_object_crds_x, set_object_crds_y):
        self.object_crds_x = set_object_crds_x
        self.object_crds_y = set_object_crds_y
        self.objects_crds_x.append(self.object_crds_x)
        self.objects_crds_y.append(self.object_crds_y)

        # optional: draw the 3 points as dots
        if self.current_object == None:
            self.current_object, = self.ax.plot(self.object_crds_x,
                                                self.object_crds_y,
                                                '.',
                                                color=self.color)
                                                #linewidth=self.linewidth)
            self.current_object.set_zorder(30)
            self.object_ax_lines.append(self.current_object)
        else:
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
        self.current_object = None

        #self.current_object.set_data(self.object_crds_x, 
        #                             self.object_crds_y)
        radius = 0.0
        center_x = 0.0
        center_y = 0.0
        try:
            radius, center_x, center_y = \
                    CircleFrom3Points2D((self.object_crds_x[0],
                                         self.object_crds_y[0]),
                                        (self.object_crds_x[1],
                                         self.object_crds_y[1]),
                                        (self.object_crds_x[2],
                                         self.object_crds_y[2]))
            sys.stdout.write('radius = '+str(radius)+
                             '  ( xcenter,ycenter = '
                             +str(center_x)+' , '+str(center_y)+' )\n')
            sys.stdout.flush()
        except np.linalg.LinAlgError as err:
            sys.stdout.write('radius = infinity  (points are colinear)\n')

        self.radii.append(radius)
        self.centers.append((center_x, center_y))

        circle = patches.Circle((center_x, center_y), 
                                radius, 
                                color=self.color, 
                                linewidth=self.linewidth,
                                fill = False)
        #circle.set_zorder(60)
        self.circles.append(circle)
        self.ax.add_patch(circle)

        self.object_crds_x = []
        self.object_crds_y = []
        self.canvas.draw()


    def DeleteLastPoint(self):
        assert(len(self.object_crds_x) == len(self.object_crds_y))
        if len(self.object_crds_x) == 1:
            del self.object_crds_x[-1]
            del self.object_crds_y[-1]
            self._DeleteAxLines_()
            
        elif len(self.object_crds_x) == 2:
            del self.object_crds_x[-1]
            del self.object_crds_y[-1]
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
        else:
            assert(len(self.object_crds_x) == 0)
            if len(self.radii) > 0:
                assert(len(self.radii) == 
                       len(self.centers) == 
                       len(self.objects_crds_x) == 
                       len(self.objects_crds_y))
                del self.radii[-1]
                del self.centers[-1]
                self.object_crds_x = self.objects_crds_x[-1][0:2]
                self.object_crds_y = self.objects_crds_y[-1][0:2]
                sys.stderr.write("  deleted previous radius measurement\n")
                del self.objects_crds_x[-1]
                del self.objects_crds_y[-1]
                assert(self.current_object == None)
                assert(len(self.object_ax_lines) > 0)
                self.current_object = self.object_ax_lines[-1]
                self.current_object.set_data(self.object_crds_x, 
                                             self.object_crds_y)
                self.circles[-1].remove()
                del self.circles[-1]
        self.canvas.draw()


    def _DeleteAxLines_(self):
        assert(len(self.object_ax_lines) > 0)
        #i = 0
        #while ((i < len(self.ax.lines)) and 
        #       (self.ax.lines[i] != self.object_ax_lines[-1])
        #    i += 1
        #if i < len(self.ax.lines):
        #    del self.ax.lines[i]
        assert(self.object_ax_lines[-1] in self.ax.lines)
        i = self.ax.lines.index(self.object_ax_lines[-1])
        del self.ax.lines[i]
        del self.object_ax_lines[-1]
        #if len(self.object_ax_lines) > 0:
        #    self.current_object = self.object_ax_lines[-1]
        #else:
        #    self.current_object = None
        self.current_object = None


    def ReportStats(self):
        if len(self.radii) == 0:
            return
        r_ave = 0.0
        for r in self.radii:
            r_ave += abs(r)
        r_ave /= len(self.radii)
        r_stddev = 0.0
        sys.stdout.write('average')
        if len(self.radii) > 1:
            sys.stdout.write('  stddev   num_measurements')
            for r in self.radii:
                r_stddev += (abs(r) - r_ave)*(abs(r) - r_ave)
            r_stddev = sqrt(r_stddev / (len(self.radii)-1))
            sys.stdout.write('\n'+str(r_ave)+' '+str(r_stddev)+' '+
                             str(len(self.radii))+'\n')
        else:
            sys.stdout.write('\n'+str(r_ave)+'\n')


    def disconnect(self):
        if self.mouse_cid:
            self.canvas.mpl_disconnect(self.mouse_cid)
        if self.key_press_cid:
            self.canvas.mpl_disconnect(self.key_press_cid)
        if self.key_release_cid:
            self.canvas.mpl_disconnect(self.key_release_cid)


    def Finalize(self):
        if len(self.object_crds_x) > 0:
            del self.object_crds_x[:]
            del self.object_crds_y[:]
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
        self.ReportStats()
        self.disconnect()


    def __del__(self):
        self.Finalize()







class SplineHighlighter:
    def __init__(self, ax,
                 set_x_interp = [],
                 set_y_interp = [],
                 set_t_interp = [],
                 set_x_clickable = [],
                 set_y_clickable = []):

        self.object_crds_x = [] # x,y coordinates of the most recently
        self.object_crds_y = [] # selected point (not yet part of a line)
        self.object_crds_t = [] # t (spline paramter) of selected point 
        self.objects_crds_x = [] # a list of coords for pairs of points
        self.objects_crds_y = [] # that already belong to lines
        self.objects_crds_t = [] # (and t parameter values)
        self.lengths = []   # A list of contour lengths
        self.x_interp = set_x_interp
        self.y_interp = set_y_interp
        self.t_interp = set_t_interp
        self.x_clickable = set_x_clickable
        self.y_clickable = set_y_clickable
        # The remaining data members are specific to matplotlib:
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.object_ax_lines = []
        self.color = (0.0, 0.28, 0.0)
        self.linewidth = 4.0
        #self.connect()

        self.mouse_cid    = None
        self.key_press_cid = None
        self.key_release_cid = None


    def connect(self):
        self.mouse_cid    = self.canvas.mpl_connect('button_press_event',
                                                    self.handle_mouse_press)
        self.key_press_cid = self.canvas.mpl_connect('key_press_event',
                                                     self.handle_key_press)
        self.key_release_cid = self.canvas.mpl_connect('key_release_event',
                                                       self.handle_key_release)

    def handle_mouse_press(self, event):
        if event.inaxes!=self.ax:
            return
        if event.button != 2:
            return
            
        if len(self.x_interp) == 0:
            return
        else:
            if len(self.x_clickable) == 0:
                x = event.xdata
                y = event.ydata
            else:
                i_min_dist = FindClosest2D(event.xdata, event.ydata, 
                                           self.x_clickable, self.y_clickable)
                x = self.x_clickable[i_min_dist]
                y = self.y_clickable[i_min_dist]
                
            i_min_dist = FindClosest2D(x, y,
                                       self.x_interp, self.y_interp)
            x = self.x_interp[i_min_dist]
            y = self.y_interp[i_min_dist]
            t = self.t_interp[i_min_dist]

        self.object_crds_x.append(x)
        self.object_crds_y.append(y)
        self.object_crds_t.append(t)

        # The next 6 lines are for debugging only
        if not (len(self.object_crds_x) ==
                len(self.object_crds_y) ==
                len(self.object_crds_t)):
            print(self.object_crds_x)
            print(self.object_crds_y)
            print(self.object_crds_t)

        assert(len(self.object_crds_x) ==
               len(self.object_crds_y) ==
               len(self.object_crds_t))

        if len(self.object_crds_x) == 1:
            self.current_object, = self.ax.plot(self.object_crds_x,
                                                self.object_crds_y,
                                                marker='o',
                                                color=self.color,
                                                linewidth=self.linewidth)
            self.object_ax_lines.append(self.current_object)
            self.canvas.draw()

        elif len(self.object_crds_x) == 2:
            self.objects_crds_x.append(self.object_crds_x)
            self.objects_crds_y.append(self.object_crds_y)
            self.objects_crds_t.append(self.object_crds_t)
            coords_curve_x = []
            coords_curve_y = []
            length = -1.0
            i = 0
            while i < len(self.t_interp):
                if BetweenInclusive(self.object_crds_t[0],
                                    self.t_interp[i],
                                    self.object_crds_t[1]):
                    coords_curve_x.append(self.x_interp[i])
                    coords_curve_y.append(self.y_interp[i])
                    if length < 0.0:
                        length = 0.0
                    else:
                        length += Distance((coords_curve_x[-2],coords_curve_y[-2]),
                                           (coords_curve_x[-1],coords_curve_y[-1]))
                i += 1

            self.current_object.set_data(coords_curve_x,
                                         coords_curve_y)
            self.current_object.set_marker(None)
            self.lengths.append(length)
            sys.stdout.write('contour_length = '+str(length)+'\n')
            self.object_crds_x = []
            self.object_crds_y = []
            self.object_crds_t = []
            self.canvas.draw()

    def handle_key_press(self, event):
        if event.key == 'delete':
            self.DeleteLastPoint()
        if event.key == 'enter':
            #Print the average distance between points already selected (if any)
            self.ReportStats()
            # Clear all the arrays containing any existing coordinate data
            del self.object_crds_x [:] # x,y coordinates of the most recently
            del self.object_crds_y [:] # selected point (not yet part of a line)
            del self.object_crds_t [:]
            del self.objects_crds_x [:] # a list of coords for pairs of points
            del self.objects_crds_y [:] # already belonging to lines
            del self.objects_crds_t [:]
            # Clear the array of distances (between these pairs of points)
            del self.lengths [:]
            # Get rid of all the lines on the screen (between these points)
            i = len(self.object_ax_lines)-1
            while i >= 0:
                #j = len(self.ax.lines)-1
                #while (j>=0) and (self.ax.lines[j]!=self.object_ax_lines[i]):
                #    j -= 1
                #assert(j >= 0):
                j = self.ax.lines.index(self.object_ax_lines[i])
                del self.ax.lines[j]
                del self.object_ax_lines[i]
                i -= 1
            self.canvas.draw()

    def handle_key_release(self, event):
        # do nothing for now
        pass

    def DeleteLastPoint(self):
        assert(len(self.object_crds_x) == len(self.object_crds_y))
        if len(self.object_crds_x) == 1:
            del self.object_crds_x[-1]
            del self.object_crds_y[-1]
            del self.object_crds_t[-1]
            assert(len(self.object_ax_lines) > 0)
            #i = 0
            #while ((i < len(self.ax.lines)) and 
            #       (self.ax.lines[i] != self.object_ax_lines[-1])
            #    i += 1
            #if i < len(self.ax.lines):
            #    del self.ax.lines[i]
            assert(self.object_ax_lines[-1] in self.ax.lines)
            i = self.ax.lines.index(self.object_ax_lines[-1])
            del self.ax.lines[i]
            del self.object_ax_lines[-1]
            if len(self.object_ax_lines) > 0:
                self.current_object = self.object_ax_lines[-1]
        #elif len(self.object_crds_x) == 2:
        #    del self.object_crds_x[-1]
        #    del self.object_crds_y[-1]
        #    self.current_object.set_data(self.object_crds_x, 
        #                                 self.object_crds_y)
        else:
            if len(self.lengths) > 0:
                assert(len(self.lengths) == 
                       len(self.objects_crds_x) == 
                       len(self.objects_crds_y) ==
                       len(self.objects_crds_t))
                del self.lengths[-1]
                self.object_crds_x = self.objects_crds_x[-1][0:1]
                self.object_crds_y = self.objects_crds_y[-1][0:1]
                self.object_crds_t = self.objects_crds_t[-1][0:1]
                del self.objects_crds_x[-1]
                del self.objects_crds_y[-1]
                del self.objects_crds_t[-1]
                sys.stderr.write("  deleted previous curve-length measurement\n")
                self.current_object.set_data(self.object_crds_x, 
                                             self.object_crds_y)
                self.current_object.set_marker('o')
        self.canvas.draw()

    def ReportStats(self):
        if len(self.lengths) == 0:
            return
        l_ave = 0.0
        for l in self.lengths:
            l_ave += l
        l_ave /= len(self.lengths)
        l_stddev = 0.0
        sys.stdout.write('average')
        if len(self.lengths) > 1:
            sys.stdout.write('  stddev   num_measurements')
            for l in self.lengths:
                l_stddev += (l - l_ave)*(l - l_ave)
            l_stddev = sqrt(l_stddev / (len(self.lengths)-1))
            sys.stdout.write('\n'+str(l_ave)+' '+str(l_stddev)+' '+
                             str(len(self.lengths))+'\n')
        else:
            sys.stdout.write('\n'+str(l_ave)+'\n')

    def disconnect(self):
        if self.mouse_cid:
            self.canvas.mpl_disconnect(self.mouse_cid)
        if self.key_press_cid:
            self.canvas.mpl_disconnect(self.key_press_cid)
        if self.key_release_cid:
            self.canvas.mpl_disconnect(self.key_release_cid)

    def Finalize(self):
        if len(self.object_crds_x) > 0:
            del self.object_crds_x[:]
            del self.object_crds_y[:]
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
        self.ReportStats()
        self.disconnect()

    def __del__(self):
        self.Finalize()







class SplineTangentCircles:

    def __init__(self, 
                 ax,
                 set_x_interp = [],
                 set_y_interp = [],
                 set_t_interp = [],
                 set_x_control = [],
                 set_y_control = [],
                 set_t_control = [],
                 set_c3a = [],
                 set_c3b = [],
                 set_c1a = [],
                 set_c1b = [],
                 set_x_clickable = [],
                 set_y_clickable = []):

        self.object_crds_x = [] # x,y coordinates of the most recently
        self.object_crds_y = [] # selected point (not yet part of a line)
        self.object_crds_t = [] # t (spline paramter) of selected point 
        self.objects_crds_x = [] # a list of coords for pairs of points
        self.objects_crds_y = [] # that already belong to lines
        self.objects_crds_t = [] # (and t parameter values)
        self.x_interp = set_x_interp
        self.y_interp = set_y_interp
        self.t_interp = set_t_interp
        self.x_control = set_x_control
        self.y_control = set_y_control
        self.t_control = set_t_control
        self.c3a = set_c3a
        self.c3b = set_c3b
        self.c1a = set_c1a
        self.c1b = set_c1b
        self.x_clickable = set_x_clickable
        self.y_clickable = set_y_clickable
        self.radii = []
        self.centers = []
        self.circles = []
        # The remaining data members are specific to matplotlib:
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.object_ax_lines = []
        self.color = (0.9, 0.4, 0.0)
        self.linewidth = 4.0
        #self.connect()

        self.mouse_cid    = None
        self.key_press_cid = None
        self.key_release_cid = None


    def connect(self):
        self.mouse_cid    = self.canvas.mpl_connect('button_press_event',
                                                    self.handle_mouse_press)
        self.key_press_cid = self.canvas.mpl_connect('key_press_event',
                                                     self.handle_key_press)
        self.key_release_cid = self.canvas.mpl_connect('key_release_event',
                                                       self.handle_key_release)

    def handle_mouse_press(self, event):
        if event.inaxes!=self.ax:
            return
        if event.button != 2:
            return
            
        if len(self.x_interp) == 0:
            return

        if len(self.x_clickable) == 0:
            x = event.xdata
            y = event.ydata
        else:
            i_min_dist = FindClosest2D(event.xdata, event.ydata, 
                                       self.x_clickable, self.y_clickable)
            x = self.x_clickable[i_min_dist]
            y = self.y_clickable[i_min_dist]
                
        i_min_dist = FindClosest2D(x, y,
                                   self.x_interp, self.y_interp)
        x = self.x_interp[i_min_dist]
        y = self.y_interp[i_min_dist]
        t = self.t_interp[i_min_dist]

        tangent_vect = SplineInterpEvalD1(t,
                                          self.c3a, self.c3b, 
                                          self.c1a, self.c1b, 
                                          self.t_control)
        NormalizeInPlace(tangent_vect)

        # normal_vect is the tangent_vect rotated by 90 degrees
        normal_vect = (-tangent_vect[1], tangent_vect[0])

        radius_curvature = 0.0
        x0 = 0.0
        y0 = 0.0
        try:
            Kcurve = SplineInterpCurvature2D(t,
                                             self.c3a, self.c3b, 
                                             self.c1a, self.c1b, 
                                             self.t_control)
        except ZeroDivisionError as err:
            sys.stdout.write('radius = 0\n')
        try:
            radius_curvature = 1.0/Kcurve
            sys.stdout.write('radius = '+str(abs(radius_curvature))+
                             '  ( intersection x,y = '+str(x)+' , '+str(y)+' )\n')
            sys.stdout.flush()

            # Now draw a circle tangent to the curve
            x0 = self.x_interp[i_min_dist] + normal_vect[0]*radius_curvature
            y0 = self.y_interp[i_min_dist] + normal_vect[1]*radius_curvature
            #self.curvature_circle = plt.Circle((x0,y0),
            #                                   abs(radius_curvature),
            #                                   color='g', fill=False)

        except ZeroDivisionError as err:
            sys.stdout.write('radius = infinity (inflection point)\n')
            sys.stdout.flush()

        self.radii.append(radius_curvature)
        self.centers.append((x0, y0))
        circle = patches.Circle((x0, y0), 
                                radius_curvature, 
                                color=self.color, 
                                linewidth=self.linewidth,
                                fill = False)
        #circle.set_zorder(70)
        self.circles.append(circle)
        self.ax.add_patch(circle)
        # Optional: draw a point (marker) where the user clicked
        self.current_object, = self.ax.plot([x],
                                            [y],
                                            marker='.',
                                            color=self.color,
                                            linewidth=self.linewidth)
        self.current_object.set_zorder(30)
        self.object_ax_lines.append(self.current_object)

        self.canvas.draw()


    def handle_key_press(self, event):
        if event.key == 'delete':
            self.DeleteLastPoint()
        if event.key == 'enter':
            #Print the average distance between points already selected (if any)
            self.ReportStats()
            # Clear all the arrays containing any existing coordinate data
            del self.object_crds_x [:] # x,y coordinates of the most recently
            del self.object_crds_y [:] # selected point (not yet part of a line)
            del self.objects_crds_x [:] # a list of coords for pairs of points
            del self.objects_crds_y [:] # already belonging to lines
            # Clear the array of distances (between these pairs of points)
            del self.radii [:]
            del self.centers [:]
            # Get rid of all the lines on the screen (between these points)
            i = len(self.object_ax_lines)-1
            while i >= 0:
                #j = len(self.ax.lines)-1
                #while (j>=0) and (self.ax.lines[j]!=self.object_ax_lines[i]):
                #    j -= 1
                #assert(j >= 0):
                j = self.ax.lines.index(self.object_ax_lines[i])
                del self.ax.lines[j]
                del self.object_ax_lines[i]
                self.circles[i].remove()
                del self.circles[i]
                i -= 1
            self.canvas.draw()

    def handle_key_release(self, event):
        # do nothing for now
        pass

    def DeleteLastPoint(self):
        assert(len(self.object_ax_lines) ==
               len(self.radii) == 
               len(self.centers) == 
               len(self.circles))

        if len(self.radii) == 0:
            return

        del self.radii[-1]
        del self.centers[-1]
        self.circles[-1].remove()
        del self.circles[-1]
        sys.stderr.write("  deleted previous radius measurement\n")

        #i = 0
        #while ((i < len(self.ax.lines)) and 
        #       (self.ax.lines[i] != self.object_ax_lines[-1])
        #    i += 1
        #if i < len(self.ax.lines):
        #    del self.ax.lines[i]
        assert(self.object_ax_lines[-1] in self.ax.lines)
        i = self.ax.lines.index(self.object_ax_lines[-1])
        del self.ax.lines[i]
        del self.object_ax_lines[-1]
        if len(self.object_ax_lines) > 0:
            self.current_object = self.object_ax_lines[-1]
        self.canvas.draw()


    def ReportStats(self):
        if len(self.radii) == 0:
            return
        r_ave = 0.0
        for r in self.radii:
            r_ave += abs(r)
        r_ave /= len(self.radii)
        r_stddev = 0.0
        sys.stdout.write('average')
        if len(self.radii) > 1:
            sys.stdout.write('  stddev   num_measurements')
            for r in self.radii:
                r_stddev += (abs(r) - r_ave)*(abs(r) - r_ave)
            r_stddev = sqrt(r_stddev / (len(self.radii)-1))
            sys.stdout.write('\n'+str(r_ave)+' '+str(r_stddev)+' '+
                             str(len(self.radii))+'\n')
        else:
            sys.stdout.write('\n'+str(r_ave)+'\n')

    def disconnect(self):
        if self.mouse_cid:
            self.canvas.mpl_disconnect(self.mouse_cid)
        if self.key_press_cid:
            self.canvas.mpl_disconnect(self.key_press_cid)
        if self.key_release_cid:
            self.canvas.mpl_disconnect(self.key_release_cid)

    def Finalize(self):
        if len(self.object_crds_x) > 0:
            del self.object_crds_x[:]
            del self.object_crds_y[:]
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
        self.ReportStats()
        self.disconnect()

    def __del__(self):
        self.Finalize()








class SplineTangentLines:

    def __init__(self, 
                 ax,
                 set_image_physical_size,
                 set_x_interp = [],
                 set_y_interp = [],
                 set_t_interp = [],
                 set_x_control = [],
                 set_y_control = [],
                 set_t_control = [],
                 set_c3a = [],
                 set_c3b = [],
                 set_c1a = [],
                 set_c1b = [],
                 set_x_clickable = [],
                 set_y_clickable = []):

        self.object_crds_x = [] # x,y coordinates of the most recently
        self.object_crds_y = [] # selected point (not yet part of a line)
        self.object_crds_t = [] # t (spline paramter) of selected point 
        self.objects_crds_x = [] # a list of coords for pairs of points
        self.objects_crds_y = [] # that already belong to lines
        self.objects_crds_t = [] # (and t parameter values)
        self.x_interp = set_x_interp
        self.y_interp = set_y_interp
        self.t_interp = set_t_interp
        self.x_control = set_x_control
        self.y_control = set_y_control
        self.t_control = set_t_control
        self.c3a = set_c3a
        self.c3b = set_c3b
        self.c1a = set_c1a
        self.c1b = set_c1b
        self.x_clickable = set_x_clickable
        self.y_clickable = set_y_clickable
        self.image_physical_size = set_image_physical_size
        self.angles = []
        self.centers = []
        # The remaining data members are specific to matplotlib:
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.object_ax_lines = []
        self.color = (0.69, 0.0, 0.69)
        self.linewidth = 1.5
        #self.connect()

        self.mouse_cid    = None
        self.key_press_cid = None
        self.key_release_cid = None


    def connect(self):
        self.mouse_cid    = self.canvas.mpl_connect('button_press_event',
                                                    self.handle_mouse_press)
        self.key_press_cid = self.canvas.mpl_connect('key_press_event',
                                                     self.handle_key_press)
        self.key_release_cid = self.canvas.mpl_connect('key_release_event',
                                                       self.handle_key_release)

    def handle_mouse_press(self, event):
        if event.inaxes!=self.ax:
            return
        if event.button != 2:
            return
            
        if len(self.x_interp) == 0:
            return

        if len(self.x_clickable) == 0:
            x = event.xdata
            y = event.ydata
        else:
            i_min_dist = FindClosest2D(event.xdata, event.ydata, 
                                       self.x_clickable, self.y_clickable)
            x = self.x_clickable[i_min_dist]
            y = self.y_clickable[i_min_dist]
                
        i_min_dist = FindClosest2D(x, y,
                                   self.x_interp, self.y_interp)
        x = self.x_interp[i_min_dist]
        y = self.y_interp[i_min_dist]
        t = self.t_interp[i_min_dist]

        tangent_vect = SplineInterpEvalD1(t,
                                          self.c3a, self.c3b, 
                                          self.c1a, self.c1b, 
                                          self.t_control)
        NormalizeInPlace(tangent_vect)

        theta = atan2(tangent_vect[1], tangent_vect[0]) * 180.0/pi
        self.angles.append(theta)
        self.centers.append((x,y))

        sys.stdout.write('angle_from_x = '+str(theta)+
                         ' ( intersection x,y = '+str(x)+' , '+str(y)+' )\n')

        line_data_x = [self.x_interp[i_min_dist] - 
                       5.0*tangent_vect[0]*self.image_physical_size,
                       self.x_interp[i_min_dist],  # <-- visible dot
                       self.x_interp[i_min_dist] + 
                       5.0*tangent_vect[0]*self.image_physical_size ]
        line_data_y = [self.y_interp[i_min_dist] - 
                       5.0*tangent_vect[1]*self.image_physical_size,
                       self.y_interp[i_min_dist],  # <-- visible dot
                       self.y_interp[i_min_dist] + 
                       5.0*tangent_vect[1]*self.image_physical_size ]

        self.current_object, = self.ax.plot(line_data_x,
                                            line_data_y,
                                            marker='.',
                                            color=self.color,
                                            linewidth=self.linewidth)

        self.current_object.set_zorder(80)

        self.object_ax_lines.append(self.current_object)


        self.canvas.draw()


    def handle_key_press(self, event):
        if event.key == 'delete':
            self.DeleteLastPoint()
        if event.key == 'enter':
            #Print the average distance between points already selected (if any)
            self.ReportStats()
            # Clear all the arrays containing any existing coordinate data
            del self.object_crds_x [:] # x,y coordinates of the most recently
            del self.object_crds_y [:] # selected point (not yet part of a line)
            del self.objects_crds_x [:] # a list of coords for pairs of points
            del self.objects_crds_y [:] # already belonging to lines
            # Clear the array of distances (between these pairs of points)
            del self.angles [:]
            del self.centers [:]
            # Get rid of all the lines on the screen (between these points)
            i = len(self.object_ax_lines)-1
            while i >= 0:
                #j = len(self.ax.lines)-1
                #while (j>=0) and (self.ax.lines[j]!=self.object_ax_lines[i]):
                #    j -= 1
                #assert(j >= 0):
                j = self.ax.lines.index(self.object_ax_lines[i])
                del self.ax.lines[j]
                del self.object_ax_lines[i]
                i -= 1
            self.canvas.draw()

    def handle_key_release(self, event):
        # do nothing for now
        pass

    def DeleteLastPoint(self):
        assert(len(self.object_ax_lines) ==
               len(self.angles) == 
               len(self.centers))

        if len(self.angles) == 0:
            return

        del self.angles[-1]
        del self.centers[-1]
        sys.stderr.write("  deleted previous tangent-angle measurement\n")

        #i = 0
        #while ((i < len(self.ax.lines)) and 
        #       (self.ax.lines[i] != self.object_ax_lines[-1])
        #    i += 1
        #if i < len(self.ax.lines):
        #    del self.ax.lines[i]
        assert(self.object_ax_lines[-1] in self.ax.lines)
        i = self.ax.lines.index(self.object_ax_lines[-1])
        del self.ax.lines[i]
        del self.object_ax_lines[-1]
        if len(self.object_ax_lines) > 0:
            self.current_object = self.object_ax_lines[-1]
        self.canvas.draw()


    def ReportStats(self):
        if len(self.angles) == 0:
            return
        theta_ave = 0.0
        for theta in self.angles:
            theta_ave += theta
        theta_ave /= len(self.angles)
        theta_stddev = 0.0
        sys.stdout.write('average')
        if len(self.angles) > 1:
            sys.stdout.write('  stddev   num_measurements')
            for theta in self.angles:
                theta_stddev += (theta - theta_ave)*(theta - theta_ave)
            theta_stddev = sqrt(theta_stddev / (len(self.angles)-1))
            sys.stdout.write('\n'+str(theta_ave)+' '+str(theta_stddev)+' '+
                             str(len(self.angles))+'\n')
        else:
            sys.stdout.write('\n'+str(theta_ave)+'\n')

    def disconnect(self):
        if self.mouse_cid:
            self.canvas.mpl_disconnect(self.mouse_cid)
        if self.key_press_cid:
            self.canvas.mpl_disconnect(self.key_press_cid)
        if self.key_release_cid:
            self.canvas.mpl_disconnect(self.key_release_cid)

    def Finalize(self):
        if len(self.object_crds_x) > 0:
            del self.object_crds_x[:]
            del self.object_crds_y[:]
            self.current_object.set_data(self.object_crds_x, 
                                         self.object_crds_y)
        self.ReportStats()
        self.disconnect()

    def __del__(self):
        self.Finalize()









#http://matplotlib.org/users/event_handling.html

class SplineBuilder:
    def __init__(self, ax, 
                 set_nsub = 100, 
                 set_alpha = 0.5,
                 set_show_maxima = False,
                 set_show_inflections = False,
                 set_min_bump_height = 0.0,
                 set_bump_width = 0.0):

        self.nsub = set_nsub # The number of line segments which approximate
                             # each interval between control points in the
                             # spline interpolation.
        self.alpha_exponent = set_alpha
        self.bump_width = set_bump_width
        self.min_bump_height = set_min_bump_height

        self.x_control = [] # a list of the x-coordinates clicked on by the user
        self.y_control = [] # a list of the y-coordinates clicked on by the user
        self.t_control = [] # t_control[i] = "time" for ith control point.t[0]=0
                            #(This parameter is needed for spline interpolation.
                            # It is calculated from x_control and y_control.)

        # SPLINE COEFFICIENTS
        # These coefficients appear in the equation for the cubic polynomial
        # which approximates the curve in each interval between control points.
        # They are calculated from: x_control, y_control, and t_control
        self.c3a = []
        self.c3b = []
        self.c1a = []
        self.c1b = []

        # Various lists of coordinates corresponding to the spline 
        # evaluated at different places along the curve:
        self.x_interp = [] # spline evaluated many times along the curve
        self.y_interp = [] # (at even intervals between control points)
        self.t_interp = []
        self.show_maxima = set_show_maxima
        self.x_maxima = [] # location of all the places on the curve
        self.y_maxima = [] # where curvature is maximized
        self.t_maxima = []
        self.k_max_curvatures = []
        self.curve_maxima_directions = []
        self.show_inflections = set_show_inflections
        self.x_inflections = [] # location of all of the inflection points
        self.y_inflections = [] # (locations on the spline where curvature = 0)
        self.t_inflections = []


        # The remaining data members are specific to matplotlib:
        self.ax = ax
        self.canvas = ax.figure.canvas
        # The curve that appears in matplotlib
        self.curve, = ax.plot([], [], color=(0,0,1), linewidth=2.0)

        # The matplotlib object of the dots representing inflection points:
        self.inflections, = ax.plot([], [], color=(0.6,0,0.6), marker='o', linestyle='None')


        # The matplotlib object of the dots representing curvature maxima:

        #self.maxima, = ax.plot([], [], 
        #                       color=(1,0,0), marker='o', linestyle='None')
        # I want to be able to change the color of the maxima depending on their
        # curvature, so I decided to use "scatter()" instead of "plot()".
        #self.maxima = ax.scatter([], [], c=[], cmap='hsv', vmin=0, vmax=1, s=25)
        self.maxima = ax.scatter([], [], c=[], cmap='spectral', vmin=0, vmax=1, s=30)
        self.maxima.set_zorder(20)
        self.maxima_colors = []

        # The next 4 lines were originally for testing/debugging,
        # however removing them causes mysterious matplotlib glitches...
        clrs=np.array([0.2*i+0.1 for i in range(0, 5)])
        xy=np.random.rand(2, 5)
        self.maxima.set_offsets(xy)
        self.maxima.set_array(clrs)

        self.mouse_cid    = None
        self.key_press_cid = None
        self.key_release_cid = None




    def connect(self):
        self.mouse_cid    = self.canvas.mpl_connect('button_press_event',
                                                    self.handle_mouse_press)
        self.key_press_cid = self.canvas.mpl_connect('key_press_event',
                                                     self.handle_key_press)
        self.key_release_cid = self.canvas.mpl_connect('key_release_event',
                                                       self.handle_key_release)

    def handle_mouse_press(self, event):
        if event.inaxes!=self.ax:
            return
        if event.button != 2:
            return

        elif ((len(self.x_control) == 0) or 
              (event.xdata != self.x_control[-1]) or
              (event.ydata != self.y_control[-1])):

            #sys.stderr.write('added point ('+str(event.x)+','+
            #                 str(event.y)+')\n')
            self.x_control.append(event.xdata)
            self.y_control.append(event.ydata)

            self.CalcSplineParams()
            self.CalcInterp()
            self.CalcInflectionPoints()
            self.CalcCurvatureMaxima()
            self.redraw()



    def redraw(self):
        assert(len(self.x_control) == len(self.y_control))
        #self.ax.figure.clf()

        if len(self.x_control) == 0:
            self.curve.set_data([], [])
            self.inflections.set_data([], [])
            #self.maxima.set_data([],[])
            self.maxima.set_offsets([])
            self.canvas.draw()
        elif len(self.x_control) == 1:
            self.curve.set_data(self.x_control, self.y_control)
            self.curve.set_marker('.')
            self.inflections.set_data([], [])
            #self.maxima.set_data([],[])
            self.maxima.set_offsets([])
            self.canvas.draw()
        elif len(self.x_control) >= 2:
            #    (copy the coordinates into the matplotlib object)
            self.curve.set_data(self.x_interp, self.y_interp)
            self.curve.set_marker(None)

            # (copy the coordinates into the matplotlib object)
            #self.maxima.set_data(self.x_maxima, self.y_maxima)
            self.maxima.set_offsets(np.array([self.x_maxima, 
                                              self.y_maxima]).transpose())
            self.maxima.set_array(np.array(self.maxima_colors))

            # (copy the coordinates into the matplotlib object)
            self.inflections.set_data(self.x_inflections, self.y_inflections)

            # Draw the curve:
            self.canvas.draw()


    def handle_key_press(self, event):
        if event.key == 'delete':
            self.DeleteLastPoint()
        elif event.key.lower() == 'd':
            sys.stderr.write("Discard shallow curvature maxima criteria:\n"
                             "  Minimum Height (or Depth) of a bump: ")
            height = sys.stdin.readline()
            self.min_bump_height = float(height)
            sys.stderr.write("  The bump must rise above its surroundings by this distance.\n"
                             "  When considering the surroundings, how far away from the\n"
                             "  center of the bump should I search?: ")
            width = sys.stdin.readline()
            self.bump_width = float(width) * 2.0
            self.DiscardShallowMaxima()


    def DeleteLastPoint(self):
        if len(self.x_control) > 0:
            del self.x_control[-1]
            del self.y_control[-1]
            # If an inflection point is located on this segment, 
            # then delete it also
            while ((len(self.t_inflections) >= 1) and 
                   (len(self.t_control) >= 2) and 
                   (self.t_inflections[-1] > self.t_control[-2])):
                assert(self.t_inflections[-1] <= self.t_control[-1])
                del self.t_inflections[-1]
                #del self.x_inflections[-1]
                #del self.y_inflections[-1]
                # x_inflections and y_inflections are numpy arrays, so we can't use 'del'
                self.x_inflections = self.x_inflections[:-1]
                self.y_inflections = self.y_inflections[:-1]
            while ((len(self.t_maxima) >= 1) and 
                   (len(self.t_control) >= 2) and 
                   (self.t_maxima[-1] > self.t_control[-2])):
                assert(self.t_maxima[-1] <= self.t_control[-1])
                del self.t_maxima[-1]
                self.x_maxima = self.x_maxima[:-1]
                self.y_maxima = self.y_maxima[:-1]
                # x_maxima and y_maxima are numpy arrays, so we can't use 'del'
                #del self.x_maxima[-1]
                #del self.y_maxima[-1]
            # The next 3 lines are unnecessary because these arrays
            # are rebuilt from scratch each time the user clicks the mouse:
            #del self.t_control[-1]
            #del self.x_interp[-1*self.nsub:]
            #del self.y_interp[-1*self.nsub:]
            self.x_interp = self.x_interp[-1*self.nsub:]
            self.y_interp = self.y_interp[-1*self.nsub:]

        self.CalcSplineParams()
        self.CalcInterp()
        self.CalcInflectionPoints()
        self.CalcCurvatureMaxima()
        self.redraw()


    def handle_key_release(self, event):
        pass


    def CalcSplineParams(self):
        if len(self.x_control) < 2:
            return

        #copy the coordinates into an Nx2 dimensional array (x)
        x = np.array([self.x_control, self.y_control]).transpose()

        # Figure out the spline parameters: c3a, c3b, c1a, c1b:
        self.c3a, self.c3b, self.c1a, self.c1b, self.t_control = \
            CalcNaturalCubicSplineCoeffs(x, self.alpha_exponent)


    def CalcInterp(self):
        if len(self.x_control) < 2:
            return
        # Recalculate the points along the curve at regular intervals:
        sp_points = SplineInterpEvalNtimesPerPoint(self.nsub, 
                                                   self.c3a, self.c3b,
                                                   self.c1a, self.c1b, 
                                                   self.t_control)
        self.x_interp, self.y_interp = sp_points[0].transpose()
        self.t_interp = sp_points[1]




    def SaveControlPoints(self, out_file):
        assert(len(self.x_control) == len(self.y_control))
        N = len(self.x_control)
        if N == 0:
            sys.stderr.write("  Nothing to save.  (Draw a spline curve before pressing this button.)\n")
            return
        for i in range(0, N):
            out_file.write(str(self.x_control[i])+" "+
                           str(self.y_control[i])+"\n")

    def SaveInterp(self, out_file):
        assert(len(self.x_interp) == len(self.y_interp))
        N = len(self.x_interp)
        if N == 0:
            sys.stderr.write("  Nothing to save.  (Draw a spline curve before pressing this button.)\n")
            return
        for i in range(0, N):
            out_file.write(str(self.x_interp[i])+" "+
                           str(self.y_interp[i])+"\n")

    def LoadControlPoints(self, in_file):
        self.x_control = []
        self.y_control = []
        for line in in_file:
            tokens = line.split()
            self.x_control.append(float(tokens[0]))
            self.y_control.append(float(tokens[1]))
        self.CalcSplineParams()
        self.CalcInterp()
        self.CalcInflectionPoints()
        self.CalcCurvatureMaxima()
        self.redraw()



    def CalcInflectionPoints(self):
        if self.show_inflections and len(self.t_control) > 0:
            self.t_inflections = FindInflectionPoints(self.nsub, 
                                                      self.c3a, self.c3b,
                                                      self.c1a, self.c1b, 
                                                      self.t_control)
            self.x_inflections, self.y_inflections = \
                SplineInterpEvalMany(self.t_inflections,
                                     self.c3a, self.c3b, self.c1a, self.c1b, 
                                     self.t_control).transpose()
        else:
            del self.x_inflections[:]
            del self.y_inflections[:]
            del self.t_inflections[:]



    def CalcCurvatureMaxima(self):
        if self.show_maxima and len(self.t_control) > 0:
            self.t_maxima, self.k_max_curvatures, self.curve_maxima_directions = \
                             FindCurvatureMaxima(self.nsub, 
                                                 self.c3a, self.c3b,
                                                 self.c1a, self.c1b, 
                                                 self.t_control)

            self.x_maxima, self.y_maxima = \
                SplineInterpEvalMany(self.t_maxima,
                                     self.c3a, self.c3b, self.c1a, self.c1b, 
                                     self.t_control).transpose()

            # now color the maxima according to positive or negative:
            self.maxima_colors = [0.18
                                  for i in range(0, len(self.k_max_curvatures))]
            for i in range(0, len(self.k_max_curvatures)):
                if self.k_max_curvatures[i] >= 0:
                    self.maxima_colors[i] = 0.86
        else:
            self.x_maxima = np.array([])
            self.y_maxima = np.array([])
            self.t_maxima = np.array([])
            self.k_max_curvatures = np.array([])
            self.curve_maxima_directions = np.array([])
            self.maxima_colors = np.array([])




    def ShowMaxima(self, set_show_maxima):
        self.show_maxima = set_show_maxima
        if self.show_maxima:
            self.CalcCurvatureMaxima()
        else:
            self.x_maxima = np.array([])
            self.y_maxima = np.array([])
            self.t_maxima = np.array([])
            self.k_max_curvatures = np.array([])
            self.curve_maxima_directions = np.array([]) 
            self.maxima_colors = np.array([])
        self.redraw()



    def ShowInflections(self, set_show_inflections):
        self.show_inflections = set_show_inflections
        if self.show_inflections:
            self.CalcInflectionPoints()
        else:
            self.x_inflections = []
            self.y_inflections = []
            self.t_inflections = []
        self.redraw()



    def DiscardShallowMaxima(self):

        self.CalcCurvatureMaxima()

        is_discarded = IdentifyShallowPoints(self.min_bump_height,
                                             self.bump_width,
                                             self.k_max_curvatures,
                                             self.curve_maxima_directions,
                                             self.x_interp,
                                             self.y_interp,
                                             self.t_interp,
                                             self.x_maxima,
                                             self.y_maxima,
                                             self.t_maxima)

        assert(len(self.x_maxima) == len(self.y_maxima) == 
               len(self.t_maxima) == len(is_discarded))

        # OLD CODE:
        #i = len(is_discarded)-1
        #while i >= 0:
        #    if is_discarded[i]:
        #        del self.x_maxima[i]
        #        del self.y_maxima[i]
        #        del self.t_maxima[i]
        #        del self.k_max_curvatures[i]
        #        del self.curve_maxima_directions[i]
        #        del self.maxima_colors[i]
        #    i -= 1
        #
        # Unfortunately after converting these to numpy arrays, you can no longer use 'del'
        # Use the following code instead:

        x_maxima_old = [x for x in self.x_maxima]
        y_maxima_old = [y for y in self.y_maxima]
        t_maxima_old = [t for t in self.t_maxima]
        k_max_curvatures_old = [k for k in self.k_max_curvatures]
        curve_maxima_directions_old = [n for n in self.curve_maxima_directions]
        maxima_colors_old = [c for c in self.maxima_colors]

        i = 0
        j = 0
        while i < len(is_discarded):
            if not is_discarded[i]:
                self.x_maxima[j] = x_maxima_old[i]
                self.y_maxima[j] = y_maxima_old[i]
                self.t_maxima[j] = t_maxima_old[i]
                self.k_max_curvatures[j] = k_max_curvatures_old[i]
                self.curve_maxima_directions[j] = curve_maxima_directions_old[i]
                self.maxima_colors[j] = maxima_colors_old[i]
                j += 1
            i += 1
        self.x_maxima = self.x_maxima[:j]
        self.y_maxima = self.y_maxima[:j]
        self.t_maxima = self.t_maxima[:j]
        self.k_max_curvatures = self.k_max_curvatures[:j]
        self.curve_maxima_directions = self.curve_maxima_directions[:j]
        self.maxima_colors = self.maxima_colors[:j]

        self.redraw()



    def DiscardShallowIntervals(self):

        self.CalcCurvatureMaxima()

        is_discarded = IdentifyShallowIntervals(self.min_bump_height,
                                                self.bump_width,
                                                self.k_max_curvatures,
                                                self.curve_maxima_directions,
                                                self.x_interp,
                                                self.y_interp,
                                                self.t_interp,
                                                self.x_maxima,
                                                self.y_maxima,
                                                self.t_maxima)

        assert(len(self.x_maxima) == len(self.y_maxima) == 
               len(self.t_maxima) == len(is_discarded))

        # OLD CODE:
        #i = len(is_discarded)-1
        #while i >= 0:
        #    if is_discarded[i]:
        #        del self.x_maxima[i]
        #        del self.y_maxima[i]
        #        del self.t_maxima[i]
        #        del self.k_max_curvatures[i]
        #        del self.curve_maxima_directions[i]
        #        del self.maxima_colors[i]
        #    i -= 1
        #
        # Unfortunately after converting these to numpy arrays, you can no longer use 'del'
        # Use the following code instead:

        x_maxima_old = [x for x in self.x_maxima]
        y_maxima_old = [y for y in self.y_maxima]
        t_maxima_old = [t for t in self.t_maxima]
        k_max_curvatures_old = [k for k in self.k_max_curvatures]
        curve_maxima_directions_old = [n for n in self.curve_maxima_directions]
        maxima_colors_old = [c for c in self.maxima_colors]

        i = 0
        j = 0
        while i < len(is_discarded):
            if not is_discarded[i]:
                self.x_maxima[j] = x_maxima_old[i]
                self.y_maxima[j] = y_maxima_old[i]
                self.t_maxima[j] = t_maxima_old[i]
                self.k_max_curvatures[j] = k_max_curvatures_old[i]
                self.curve_maxima_directions[j] = curve_maxima_directions_old[i]
                self.maxima_colors[j] = maxima_colors_old[i]
                j += 1
            i += 1
        self.x_maxima = self.x_maxima[:j]
        self.y_maxima = self.y_maxima[:j]
        self.t_maxima = self.t_maxima[:j]
        self.k_max_curvatures = self.k_max_curvatures[:j]
        self.curve_maxima_directions = self.curve_maxima_directions[:j]
        self.maxima_colors = self.maxima_colors[:j]

        self.redraw()



    def disconnect(self):
        if self.mouse_cid:
            self.canvas.mpl_disconnect(self.mouse_cid)
        if self.key_press_cid:
            self.canvas.mpl_disconnect(self.key_press_cid)
        if self.key_release_cid:
            self.canvas.mpl_disconnect(self.key_release_cid)

    def Finalize(self):
        self.disconnect()

    def __del__(self):
        self.Finalize()













class WorldObjects:

    def SelectMode(self, label):
        # Finish up any remaining tasks with the currently selected object (if any)
        if self.mode in self.objects_by_name:
            #print('disconnecting \"'+self.mode+'\"')
            self.objects_by_name[self.mode].Finalize()

        self.mode = label

        # Then, figure out which object the user selected, and "connect" it.
        # This is necessary to allow it to respond to mousclicks, keypresses,...
        if label in self.objects_by_name:
            self.objects_by_name[label].connect()

        plt.draw()
        self.CopySplineCoords()
        self.LoadClickablePoints()

    def ShowGrid(self, label):
        self.ax.grid(True)

    def HideGrid(self, label):
        self.ax.grid(False)

    def HandleCheckButton(self, label):
        assert(label in self.checkb_status)
        self.checkb_status[label] = not self.checkb_status[label]

        #if label=='show grid':
        #    self.ax.grid(self.checkb_status[label])
        # COMMENTING OUT
        # Something about using self.ax.grid() here is glitchy.
        # I guess I'll try doing this with a regular button instead...

        if label=='show inflection points':
            self.spline.ShowInflections(self.checkb_status[label])
            self.CopySplineCoords()
            self.LoadClickablePoints()
        elif label=='show curvature maxima':
            self.spline.ShowMaxima(self.checkb_status[label])
            self.CopySplineCoords()
            self.LoadClickablePoints()
        else:
            # The 'snap to dots' check button option should be 
            # mutually exclusive with any other check-button choices.
            # (Currently this only includes 'snap to dots', 
            #  although I probably will add other choices later.)
            #if (self.checkb_status['snap to dots'] and
            #    self.checkb_status['snap to curve']):
            #    if label == 'snap to dots':
            #        self.check_snap.set_active['snap to curve']
            #        self.checkb_status['snap to curve'] = False
            #    elif label == 'snap to curve':
            #        self.check_snap.set_active['snap to dots']
            #        self.checkb_status['snap to dots'] = False

            self.LoadClickablePoints()

    def LoadClickablePoints(self):
        del self.x_clickable[:]
        del self.y_clickable[:]
        del self.t_clickable[:]
        if self.checkb_status['snap to dots']:
            #self.x_clickable = (self.spline.x_maxima + 
            #                   self.spline.x_inflections)
            #self.y_clickable = (self.spline.y_maxima + 
            #                   self.spline.y_inflections)
            #self.t_clickable = (self.spline.t_maxima + 
            #                   self.spline.t_inflections)
            for i in range(0, len(self.spline.x_maxima)):
                self.x_clickable.append(self.spline.x_maxima[i])
                self.y_clickable.append(self.spline.y_maxima[i])
                self.t_clickable.append(self.spline.t_maxima[i])
            for i in range(0, len(self.spline.x_inflections)):
                self.x_clickable.append(self.spline.x_inflections[i])
                self.y_clickable.append(self.spline.y_inflections[i])
                self.t_clickable.append(self.spline.t_inflections[i])
        elif self.checkb_status['snap to curve']:
            for i in range(0, len(self.spline.x_interp)):
                self.x_clickable.append(self.spline.x_interp[i])
                self.y_clickable.append(self.spline.y_interp[i])
                self.t_clickable.append(self.spline.t_interp[i])

    def CopySplineCoords(self):

        self.spline_highlighter.x_interp = self.spline.x_interp
        self.spline_highlighter.y_interp = self.spline.y_interp
        self.spline_highlighter.t_interp = self.spline.t_interp

        self.spline_tangent_circles.x_interp = self.spline.x_interp
        self.spline_tangent_circles.y_interp = self.spline.y_interp
        self.spline_tangent_circles.t_interp = self.spline.t_interp
        self.spline_tangent_circles.x_control = self.spline.x_control
        self.spline_tangent_circles.y_control = self.spline.y_control
        self.spline_tangent_circles.t_control = self.spline.t_control
        self.spline_tangent_circles.c3a = self.spline.c3a
        self.spline_tangent_circles.c3b = self.spline.c3b
        self.spline_tangent_circles.c1a = self.spline.c1a
        self.spline_tangent_circles.c1b = self.spline.c1b

        self.spline_tangent_lines.x_interp = self.spline.x_interp
        self.spline_tangent_lines.y_interp = self.spline.y_interp
        self.spline_tangent_lines.t_interp = self.spline.t_interp
        self.spline_tangent_lines.x_control = self.spline.x_control
        self.spline_tangent_lines.y_control = self.spline.y_control
        self.spline_tangent_lines.t_control = self.spline.t_control
        self.spline_tangent_lines.c3a = self.spline.c3a
        self.spline_tangent_lines.c3b = self.spline.c3b
        self.spline_tangent_lines.c1a = self.spline.c1a
        self.spline_tangent_lines.c1b = self.spline.c1b


    def DeleteLastPoint(self, label):
        if self.mode in self.objects_by_name:
            self.objects_by_name[self.mode].DeleteLastPoint()

    def DiscardShallowIntervals(self, label):
        self.spline.DiscardShallowIntervals()

    def DiscardShallowMaxima(self, label):
        self.spline.DiscardShallowMaxima()

    def SaveSpline(self, label):
        sys.stderr.write("Enter a filename for writing: ")
        filename = sys.stdin.readline().strip()
        out_file = None
        try:
            in_file = open(filename, 'r')
            in_file.close()
            sys.stderr.write("File exists. Overwrite? (Y/N): ")
            response = sys.stdin.readline().strip()
            if (len(response) > 0) and (response[0].lower() == 'y'):
                out_file = open(filename, 'w')
        except IOError:
            out_file = open(filename, 'w')
        if out_file:
            self.spline.SaveControlPoints(out_file)
            out_file.close()

    def SaveInterp(self, label):
        sys.stderr.write("Enter a filename for writing: ")
        filename = sys.stdin.readline().strip()
        out_file = None
        try:
            in_file = open(filename, 'r')
            in_file.close()
            sys.stderr.write("File exists. Overwrite? (Y/N): ")
            response = sys.stdin.readline().strip()
            if (len(response) > 0) and (response[0].lower() == 'y'):
                out_file = open(filename, 'w')
        except IOError:
            out_file = open(filename, 'w')
        if out_file:
            self.spline.SaveInterp(out_file)
            out_file.close()

    def LoadSpline(self, label):
        sys.stderr.write("Enter a filename for reading: ")
        filename = sys.stdin.readline().strip()
        try:
            in_file = open(filename, 'r')
            self.spline.LoadControlPoints(in_file)
            in_file.close()
        except IOError: 
            raise InputError('Error: Unable to open file \"'+filename+'\"\n'
                             '       for reading.\n')

    def __init__(self, img, fig, ax,
                 set_nsub = 50,
                 pixel_width = 1.0,
                 alpha_exponent = 0.5,
                 set_min_bump_height = 0.0,
                 set_bump_width = 0.0):
                 #set_show_grid = False):

        self.img = img
        self.fig = fig
        self.ax = ax
        self.nsub = set_nsub
        self.pixel_width = pixel_width
        self.alpha_exponent = alpha_exponent
        self.im=ax.imshow(self.img, interpolation = 'bicubic', 
                          extent=[0, self.img.shape[1]*pixel_width, 
                                  0, self.img.shape[0]*pixel_width])
        self.image_physical_size = (self.img.shape[0]+self.img.shape[1])*pixel_width
        self.ax.set_title('Click with the MIDDLE mouse-button to draw/measure')

        plt.subplots_adjust(left=0.47)

        self.ax.autoscale(enable=False)


        self.ax_chk = plt.axes([0.02, 0.14, 0.36, 0.23])
        self.checkb_status = {'show curvature maxima':False,
                              'show inflection points':False,
                              'snap to curve':False, 
                              'snap to dots':False}
                              #'show grid':set_show_grid}
        selfcheck_buttons = CheckButtons(self.ax_chk, 
                                         ('show curvature maxima',
                                          'show inflection points',
                                          'snap to curve',
                                          'snap to dots'),
                                         #'show grid'), 
                                         (self.checkb_status['show curvature maxima'],
                                          self.checkb_status['show inflection points'],
                                          self.checkb_status['snap to curve'],
                                          self.checkb_status['snap to dots']))
                                          #self.checkb_status['show grid']))

        selfcheck_buttons.on_clicked(self.HandleCheckButton)
        self.x_clickable = []
        self.y_clickable = []
        self.t_clickable = []


        self.ax_rad = plt.axes([0.02, 0.38, 0.36, 0.47])
        self.radio = RadioButtons(self.ax_rad, 
                                  ('measure distances', 
                                   'measure angles',
                                   'measure triangle heights',
                                   'measure 3point curvature', 
                                   'draw spline', 
                                   'measure spline curvature', 
                                   'measure spline tangents', 
                                   'measure curve lengths'))
        self.mode = 'measure distances'
        self.radio.on_clicked(self.SelectMode)
        self.spline = SplineBuilder(self.ax, 
                                    self.nsub, 
                                    self.alpha_exponent,
                                    self.checkb_status['show curvature maxima'],
                                    self.checkb_status['show inflection points'],
                                    set_min_bump_height,
                                    set_bump_width)

        self.spline_highlighter = SplineHighlighter(self.ax, 
                                                    self.spline.x_interp,
                                                    self.spline.y_interp,
                                                    self.spline.t_interp,
                                                    self.x_clickable,
                                                    self.y_clickable)
        self.spline_tangent_circles = \
            SplineTangentCircles(self.ax,
                                 self.spline.x_interp,
                                 self.spline.y_interp,
                                 self.spline.t_interp,
                                 self.spline.x_control,
                                 self.spline.y_control,
                                 self.spline.t_control,
                                 self.spline.c3a,
                                 self.spline.c3b,
                                 self.spline.c1a,
                                 self.spline.c1b,
                                 self.x_clickable,
                                 self.y_clickable)

        self.spline_tangent_lines = \
            SplineTangentLines(self.ax,
                               self.image_physical_size,
                               self.spline.x_interp,
                               self.spline.y_interp,
                               self.spline.t_interp,
                               self.spline.x_control,
                               self.spline.y_control,
                               self.spline.t_control,
                               self.spline.c3a,
                               self.spline.c3b,
                               self.spline.c1a,
                               self.spline.c1b,
                               self.x_clickable,
                               self.y_clickable)


        self.objects_by_name = {'draw spline': self.spline,
                                'measure distances':
                                LineSegmentBuilder(self.ax,
                                                   self.x_clickable,
                                                   self.y_clickable),
                                'measure angles':
                                TriangleBuilder(self.ax,
                                                self.x_clickable,
                                                self.y_clickable),
                                'measure triangle heights':
                                TriangleHeights(self.ax,
                                                self.x_clickable,
                                                self.y_clickable),
                                'measure 3point curvature':
                                ThreePointCircleBuilder(self.ax,
                                                        self.x_clickable,
                                                        self.y_clickable),
                                'measure spline curvature': self.spline_tangent_circles,
                                'measure spline tangents': self.spline_tangent_lines,
                                'measure curve lengths':self.spline_highlighter}


        self.SelectMode(self.mode)

        self.ax_del = plt.axes([0.02, 0.02, 0.07, 0.10])
        self.b_del = Button(self.ax_del, 'delete\nlast\npoint')
        self.b_del.on_clicked(self.DeleteLastPoint)

        self.ax_savespl = plt.axes([0.10, 0.02, 0.07, 0.10])
        self.b_savespl = Button(self.ax_savespl, 'save\nspline')
        self.b_savespl.on_clicked(self.SaveSpline)

        self.ax_loadspl = plt.axes([0.18, 0.02, 0.07, 0.10])
        self.b_loadspl = Button(self.ax_loadspl, 'load\nspline')
        self.b_loadspl.on_clicked(self.LoadSpline)

        if (set_min_bump_height > 0):
            self.ax_discard_maxima = plt.axes([0.26, 0.02, 0.09, 0.10])
            self.b_discard_maxima = Button(self.ax_discard_maxima,
                                           'discard\nshallow\nmaxima')
            self.b_discard_maxima.on_clicked(self.DiscardShallowMaxima)

            #self.ax_discard_intervals = plt.axes([0.38, 0.02, 0.14, 0.10])
            #self.b_discard_intervals = Button(self.ax_discard_intervals,
            #                                  'DISCARD\nSHALLOW\nINTERVALS')
            #self.b_discard_intervals.on_clicked(self.DiscardShallowIntervals)

        self.ax_show_grid = plt.axes([0.36, 0.02, 0.05, 0.10])
        self.b_show_grid = Button(self.ax_show_grid,'show\ngrid')
        self.b_show_grid.on_clicked(self.ShowGrid)

        self.ax_hide_grid = plt.axes([0.42, 0.02, 0.05, 0.10])
        self.b_hide_grid = Button(self.ax_hide_grid,'hide\ngrid')
        self.b_hide_grid.on_clicked(self.HideGrid)

        self.ax_saveallpoints = plt.axes([0.48, 0.02, 0.09, 0.10])
        self.b_saveallpoints = Button(self.ax_saveallpoints, 'save\nall points')
        self.b_saveallpoints.on_clicked(self.SaveInterp)

        self.ax.grid(True)

        plt.show()


    def Finalize(self):
        if self.mode in self.objects_by_name:
            self.objects_by_name[label].Finalize()





def main():

    usage_example = "    "+g_program_name+" image_file [-p pixel_width] [-discard h w] [-nspline n]"

    # The resolution of each spline is determined by "nsub"
    # Each interval between mouse-clicks is divided into "nsub" intervals:
    nsub = 50
    # The width of each pixel (in nm, for example).  By default, set to 1
    pixel_width = 1.0
    # The alpha_exponent is a parameter that effects spline interpolation.
    # A value of 0.5 corresponds to centripital Catmull-Rom spline.
    # A value of 1.0 corresponds to chordal Catmull-Rom spline.
    alpha_exponent = 0.5

    # In order to not be discarded, a "bump" (curvature maxima), must have
    # stray from its straight-line approximation by a distance of at least
    # "min_bump_height".  This must happen somewhere along the curve within
    # an interval of length "bump_width" centered at the curvature maxima.
    min_bump_height = 4.0
    bump_width = 450.0


    #show_grid = True
    #grid_h_spacing = 0.0
    #grid_h_offset = 0.0
    #grid_v_spacing = 0.0
    #grid_v_offset = 0.0


    #######  Main Code Below: #######

    sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+' ')
    if sys.version < '3':
        sys.stderr.write(' (python version < 3)\n')
    else:
        sys.stderr.write('\n')

    try:

        argv = [arg for arg in sys.argv]

        # Loop over the remaining arguments not processed yet.
        # These arguments are specific to the lttree.py program
        # and are not understood by ttree.py:
        i = 1
        while i < len(argv):
            #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')
            if ((argv[i].lower() == '-?') or
                (argv[i].lower() == '--?') or
                (argv[i].lower() == '-help') or
                (argv[i].lower() == '-help')):
                if i+1 >= len(argv):
                    sys.stdout.write(usage_example+'\n')
                    sys.exit(0)

            elif argv[i].lower() == '-p':
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number.\n')
                pixel_width = float(argv[i+1])
                sys.stderr.write('   pixel_width = '+str(pixel_width)+'\n')
                del(argv[i:i+2])

            elif argv[i].lower() == '-alpha':
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number.\n')
                alpha_exponent = float(argv[i+1])
                sys.stderr.write('   alpha = '+str(alpha_exponent)+'\n')
                del(argv[i:i+2])

            elif argv[i].lower() == '-nspline':
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number.\n')
                nsub = int(argv[i+1])
                sys.stderr.write('   spline resolution = '+str(nsub)+'\n')
                del(argv[i:i+2])

            elif argv[i].lower() == '-discard':
                if i+2 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by two numbers.\n')
                min_bump_height = float(argv[i+1])
                bump_width = float(argv[i+2]) * 2
                del(argv[i:i+3])


            #elif argv[i].lower() == '-grid':
            #    show_grid = True
            #    del(argv[i:i+1])
            #elif argv[i].lower() == '-nogrid':
            #    show_grid = False
            #    del(argv[i:i+1])
            #
            #elif argv[i].lower() == '-hgrid':
            #    if i+2 >= len(argv):
            #        raise InputError('Error: '+argv[i]+' flag should be followed by two numbers.\n')
            #    grid_h_spacing = float(argv[i+1])
            #    grid_h_offset  = float(argv[i+2])
            #    sys.stderr.write('   minimum_bump_height = '+str(min_bump_height)+'   interval_of_width = '+str(bump_width)+'\n')
            #    del(argv[i:i+3])
            #
            #elif argv[i].lower() == '-vgrid':
            #    if i+2 >= len(argv):
            #        raise InputError('Error: '+argv[i]+' flag should be followed by two numbers.\n')
            #    grid_v_spacing = float(argv[i+1])
            #    grid_v_offset  = float(argv[i+2])
            #    sys.stderr.write('   minimum_bump_height = '+str(min_bump_height)+'   interval_of_width = '+str(bump_width)+'\n')
            #    del(argv[i:i+3])

            elif argv[i][0] == '-':
                raise InputError('Error('+g_program_name+'):\n'
                                 'Unrecogized command line argument \"'+argv[i]+'\"\n')
            else:
                i += 1

        if len(argv) == 1:
            raise InputError("Error: Expected an image file name.\n\n"+
                             "Usage: \n\n"+
                             "       "+usage_example+"\n")

        if len(argv) > 2:
            raise InputError("Error('+g_program_name+'):\n"+
                             "       More arguments than expected.\n\n"+
                             "Usage: \n\n"+
                             "       "+usage_example+"\n")


        filename_in = argv[1]

        if min_bump_height > 0.0:
            sys.stderr.write(' criteria for discarding shallow bumps:\n')
            sys.stderr.write('   minimum_bump_height = '+str(min_bump_height)+'   search_radius = '+str(bump_width/2)+'\n')


        # ------------ Done parsing argument list ----------

        # Load a color image in grayscale
        err_msg = "Error: Invalid or unsupported image file: \""+filename_in+"\"\n" + \
                  "       (Either convert the image to PNG format, or install the \"cv2\" python module.)\n"

        try:
            import cv2
            img = cv2.imread(filename_in)
            b,g,r = cv2.split(img)
            img = cv2.merge([r,g,b])
        except ImportError:
            try:
                import matplotlib.image as mpimg
                img = mpimg.imread(filename_in)
            except IOError:
                raise InputError(err_msg)

        fig, ax = plt.subplots()

        # "world_objects" contains all of the visible objects in the plot,
        # as well as event handlers for responding to mouse-clicks

        world_objects = WorldObjects(img, fig, ax,
                                     nsub, pixel_width, alpha_exponent,
                                     min_bump_height, bump_width)
                                     #show_grid)




    except (ValueError, InputError) as err:
        sys.stderr.write('\n'+str(err)+'\n')
        sys.exit(-1)




if __name__ == "__main__":
    main()
