#! /usr/bin/python

import ctypes as ct
import numpy as np
from scipy.interpolate import splrep, splev, interp1d
import os
import time
import sys
lib = ct.cdll.LoadLibrary(os.path.dirname(__file__)+'/lib/_tools.so')

DOC_PATH = '../Fig/'


mod = np.mod
sin = np.sin
cos = np.cos
pi = np.pi
PI2 = 2.*pi


#===


# for a given array 'x' of voltages, and a 'threshhold', returns in an array 'ti', the indices j of x where crossings take place i.e.,
#          where x[j] < threshhold and x[j+1] > threshhold
lib.crossings.argtypes = [ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_double]
def crossings(x, threshold):
    x = np.array(x)
    ti = x.size*np.ones((1+x.size/2), float)	 # this is definitively the maximum number of possible crossings

    lib.crossings(x.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_uint(x.size),
        ti.ctypes.data_as(ct.POINTER(ct.c_double)),
        ct.c_double(threshold))

    return ti[:ti.searchsorted(x.size)]


def find_crossings(x, trigger):
    ti = [crossings(x[j], trigger) for j in xrange(x.shape[0])]
    return ti


def unwrapped_phase(recurrences, times):
        f = interp1d(recurrences, PI2*np.arange(recurrences.size))
        return f(times)


def unwrapped_phase_differences(recurrences, reference_index=None):

        start = np.max([t_j[0] for t_j in recurrences]) # last of the first recurrence events,
        end = np.min([t_j[-1] for t_j in recurrences])  # first of the last recurrence events:        in between, all phases can be interpolated

        if reference_index == None:
                times = np.sort(np.concatenate(recurrences))
                times = times[times.searchsorted(start):times.searchsorted(end)]    # these are all times in between

        else:
            t_j = np.asarray(recurrences[reference_index])
            times = t_j[t_j.searchsorted(start):t_j.searchsorted(end)]

        phases = np.array([unwrapped_phase(t_j, times) for t_j in recurrences])

        phase_differences = np.array([phases[i]-phases[0] for i in xrange(1, len(recurrences), 1)])
        for i in xrange(phase_differences.shape[0]):
                phase_differences[i] -= phase_differences[i, 0]  # differences should start diffusing at zero

        return times, phase_differences



"""
ti has indices of crossings for all the cells, row by row
first row - j1,j2,j3,j4 ... where j's refer to the indices of V[0][] where V[0][j] < threshhold and V[0][j+1] > threshhold
"""
def compute_phase_difference(ti):
    #compute for each of the cells, the indices to insert into its crossings, the crossings of the first cell, to maintain a sorted order
    idx = [ti[i].searchsorted(ti[0]) for i in xrange(1, len(ti), 1)]	 # idx contains closest crossings

    #compute inverses of j2-j1,j3-j2,j4-j3 ...
    #dti0_inv = 1./(ti[0][1:]-ti[0][:-1]) 	# period for normalization

    d, Ti = [], []
    for i in xrange(ti[0].size-1):	 # for each return time of Neuron 1

 #       try:
            d.append([])
            Ti.append(ti[0][i])

            #original wrong????
            #d[-1].append((ti[j][idx[j-1][i]]-ti[0][i])*dti0_inv[i])
            #TODO edited right????

            for j in xrange(1, len(ti), 1):
                indx = idx[j-1][i]
                if indx == 0 | indx == len(ti[j]) :
                    d.pop(); Ti.pop();
                    break;
                else:
                    d[-1].append((ti[0][i]-ti[j][indx-1])*1./(ti[j][indx]-ti[j][indx-1]));

#        except:
#            d.pop()
#            Ti.pop()

    return np.asarray(Ti), np.asarray(d)



def phase_difference(V, V_trigger=-40.): # V should be in millivolt!
    ti = find_crossings(V, V_trigger)
    return compute_phase_difference(ti)



#===


class splineLS1D():
    def __init__(self, s=0., isphase=False, verbose=0):
        self.isphase = isphase
        self.s = s
        self.tck = None
        self.verbose = verbose

    def set(self, what, value):
        print "# Setting", what, "to", value
        if what == 's':
            self.s = float(value)
        if what == 'isphase':
            self.isphase = value
            if value:
                self.xmin, self.xmax = 0., 2.*np.pi
        else:
            print "#", what, "cannot be set to", value
            exit(-1)

    def get(self, what):
        if what == 'params':
            return self.tck

    def makeModel(self, f, x, w=None):
        #print "splineLS1D: computing smoothing spline to the data ..."
        f, x = np.asarray(f), np.asarray(x)

        if self.isphase:	# we have to append x=0, and x=2*pi, in order to have the algorithm compute 2*pi periodic values
            self.xmin, self.xmax=0., 2.*np.pi
            x = np.mod(x, 2.*np.pi)
            fmin, fmax = f[x.argmin()], f[x.argmax()]
            xmin, xmax = x.min(), x.max()
            f, x = list(f), list(x)
            if self.xmin != xmin: f.append(fmin-xmin/(xmin-(xmax-2.*np.pi))*(fmin-fmax)), x.append(self.xmin)
            if self.xmax != xmax: f.append(fmax+(2.*np.pi-xmax)/(xmin-(xmax-2.*np.pi))*(fmin-fmax)), x.append(self.xmax)
            f, x = np.array(f), np.array(x)

            if np.iterable(w):
                w = list(w)
                w.extend([1.])
                w = np.asarray(w)

        else:
            self.xmin = x.min(); self.xmax = x.max()

        index = np.argsort(x)
        xs, fs = x[index], f[index]
        self.tck = splrep(xs, fs, w=w, s=self.s, per=self.isphase)

        if self.verbose:
            print "computation done"
            figure()
            plot(x, f, 'k,')
            plot(xs, self(xs), 'k--')

            if self.isphase:
                ph.xticks(fontsize=20)

            if self.verbose == 1:
                show()

    def __call__(self, x):

        if self.isphase:
            return splev(np.mod(x, 2.*np.pi), self.tck)

        if self.tck:
            XMIN = x>self.xmin; XMAX = x<self.xmax
            return splev(x*XMIN*XMAX+XMAX*(1.-XMIN)*self.xmin+XMIN*(1.-XMAX)*self.xmax, self.tck)

        else:
            return np.zeros((x.size), float)

    def df(self, x):
        if self.isphase:
            return spalde(np.mod(x, 2.*np.pi), self.tck)[:, 1]

        return np.array(spalde(x, self.tck))[:, 1]

    def saveCoefs(self, filename=None, file=None, close=True):
        p = self.get('params')
        if not file: file = open(filename, 'w')
        #print "# ", self.tck
        file.write(repr(self.xmin)+','+repr(self.xmax)+','+str(self.tck[0].size)+',\n')
        file.write(repr(self.tck[0][0]))
        for i in range(1, self.tck[0].size):
            file.write(','+repr(self.tck[0][i]))
        file.write('\n')
        file.write(repr(self.tck[1][0]))
        for i in range(1, self.tck[1].size):
            file.write(','+repr(self.tck[1][i]))
        if close: file.close()
        else:
            file.write('\n')
            return file

    def loadCoefs(self, filename=None, file=None, close=True):
        self.pindex = []
        p = []
        if not file: file = open(filename, 'r')
        tmp = file.readline().split(',')
        self.xmin = float(tmp[0]); self.xmax = float(tmp[1]); Nc = int(tmp[2])
        c, t = [], []
        tmp = file.readline().split(',')
        for i in range(Nc):
            c.append(float(tmp[i]))
        tmp = file.readline().split(',')
        for i in range(Nc):
            t.append(float(tmp[i]))
        self.tck = (c, t, 3)
        if close: file.close()
        else: return file

############################################################
########################## PLOTTING ########################
############################################################


def adjustForPlotting(x, y, ratio, threshold):	# ratio = xscale/yscale

    x, y = np.asarray(x), np.asarray(y)

    dx, dy = x[1:]-x[:-1], ratio*(y[1:]-y[:-1])

    xnew, ynew = [x[0]], [y[0]]
    for i in xrange(1, x.size, 1):
        if np.sqrt((x[i]-xnew[-1])**2 + (ratio*(y[i]-ynew[-1]))**2) > threshold:
            xnew.append(x[i])
            ynew.append(y[i])

    return np.asarray(xnew), np.asarray(ynew)




def plot_phase_2D(phase_1, phase_2, **kwargs):
    from pylab import plot, subplot

    try:
        PI = kwargs.pop('PI')

    except:
        PI = np.pi

    try:
        ax = kwargs.pop('axes')

    except:
        ax = subplot(111)

    j0 = 0

    for j in xrange(1, phase_1.size):

        if abs(phase_1[j]-phase_1[j-1]) < PI and abs(phase_2[j]-phase_2[j-1]) < PI:
            continue

        else:

            try:
                ax.plot(phase_1[j0:j], phase_2[j0:j], '-', **kwargs)

            except:
                pass

            j0 = j

    try:
        ax.plot(phase_1[j0:], phase_2[j0:], '-', **kwargs)

    except:
        pass


def plot_phase_3D(phase_1, phase_2, phase_3, ax, **kwargs):
    from pylab import plot, subplot

    try:
        PI = kwargs.pop('PI')

    except:
        PI = np.pi

    #assert isinstance(Axes3D)

    j0 = 0

    for j in xrange(1, phase_1.size):

        if abs(phase_1[j]-phase_1[j-1]) < PI and abs(phase_2[j]-phase_2[j-1]) < PI and abs(phase_3[j]-phase_3[j-1]) < PI:
            continue

        else:

            try:
                ax.plot(phase_1[j0:j], phase_2[j0:j], phase_3[j0:j], '-', **kwargs)

            except:
                pass

            j0 = j

    try:
        ax.plot(phase_1[j0:], phase_2[j0:], phase_3[j0:], '-', **kwargs)

    except:
        pass



rb = 0.4
gb = 0.4
bb = 0.4
sh1 = pi/4.+0.5
sh2 = -pi/4.-0.5


def clmap_bin(th1, th2):
    r = sin((th1+sh1)/2.)**6*sin((th2+sh2)/2.)**6+cos(th2/2.)**6 + 0.1*(sin(th1/2.)**6*sin(th2/2.)**6)
    g = cos(th1/2.)**6+0.3*(sin(th1/2.)**6*sin(th2/2.)**6) + 0.1*(sin(th1/2.)**6*sin(th2/2.)**6)
    #b = sin((th1+sh1)/2.)**6*sin((th2+sh2)/2.)**6 + sin(th1/2.)**6*sin(th2/2.)**6
    r = cos(.5*(th1))*cos(0.5*(th2-pi))
    g = cos(.5*(th1-pi))*cos(0.5*(th2))
    b = cos(.5*(th1-pi))*cos(0.5*(th2-pi))
    k = cos(.5*(th1-4.*pi/3.))*cos(0.5*(th2-2.*pi/3.))
    p = cos(.5*(th1-2.*pi/3.))*cos(0.5*(th2-4.*pi/3.))
    x = np.zeros((3), float)
    ncolor = argmax(array([r, g, b, k**2, p**2])**2)
    if ncolor < 3:
        x[ncolor] = 1.
        return (x[0], x[1], x[2], 1.)
    else:
        if ncolor == 3:
            return (0., 0., 0., 1.)
        if ncolor == 4:
            return (1., 0., 1., 1.)


def clmap_bin_river(th1, th2):
    r = sin((th1+sh1)/2.)**6*sin((th2+sh2)/2.)**6+cos(th2/2.)**6 + 0.1*(sin(th1/2.)**6*sin(th2/2.)**6)
    g = cos(th1/2.)**6+0.3*(sin(th1/2.)**6*sin(th2/2.)**6) + 0.1*(sin(th1/2.)**6*sin(th2/2.)**6)
    b = sin((th1+sh1)/2.)**6*sin((th2+sh2)/2.)**6 + sin(th1/2.)**6*sin(th2/2.)**6
    x = np.zeros((3), float)
    n = argmax([r, g, b])
    if n == 0: x[0] = 1.
    else: x = x+0.6
    return (x[0], x[1], x[2], 1.)


def clmap_bin_crossriver(th1, th2):
    r = sin((th1+sh1)/2.)**6*sin((th2+sh2)/2.)**6+cos(th2/2.)**6 + 0.1*(sin(th1/2.)**6*sin(th2/2.)**6)
    g = cos(th1/2.)**6+0.3*(sin(th1/2.)**6*sin(th2/2.)**6) + 0.1*(sin(th1/2.)**6*sin(th2/2.)**6)
    b = sin((th1+sh1)/2.)**6*sin((th2+sh2)/2.)**6 + sin(th1/2.)**6*sin(th2/2.)**6
    x = np.zeros((3), float)
    n = argmax([r, g, b])
    if n == 1: x[1] = 1.
    else: x = x+0.6
    return (x[0], x[1], x[2], 1.)


def clmap(th1, th2):
    r = sin((th1+sh1)/2.)**6*sin((th2+sh2)/2.)**6+cos(th2/2.)**6 + 0.1*(sin(th1/2.)**6*sin(th2/2.)**6)
    g = cos(th1/2.)**6+0.3*(sin(th1/2.)**6*sin(th2/2.)**6) + 0.1*(sin(th1/2.)**6*sin(th2/2.)**6)
    b = sin((th1+sh1)/2.)**6*sin((th2+sh2)/2.)**6 + sin(th1/2.)**6*sin(th2/2.)**6
    return (r/(1.+rb), g/(1.+gb), b/(1.+bb), 1.)


def clmap_grey(th1, th2):
    if th2 > np.pi/4. and th2 < 7.*np.pi/4.:
        return (0.2, 0.2, 0.2, 1.)

    else:
        return clmap(th1, th2)


def clmap2(th1, th2):
    th1, th2 = mod(7.*th1, 2.*pi), mod(7.*th2, 2.*pi)
    r = sin((th1+sh1)/2.)**6*sin((th2+sh2)/2.)**6+cos(th2/2.)**6 + 0.1*(sin(th1/2.)**6*sin(th2/2.)**6)
    g = cos(th1/2.)**6+0.3*(sin(th1/2.)**6*sin(th2/2.)**6) + 0.1*(sin(th1/2.)**6*sin(th2/2.)**6)
    b = sin((th1+sh1)/2.)**6*sin((th2+sh2)/2.)**6 + sin(th1/2.)**6*sin(th2/2.)**6
    return (r/(1.+rb), g/(1.+gb), b/(1.+bb), 1.)


def torus_distance(t_1, t_2): # phases range (0, 1)
    delta = abs(t_1-t_2)
    return np.sqrt((min(delta[0], 1.-delta[0])**2+min(delta[1], 1.-delta[1])**2).sum())


torus_patterns = [np.array([0.5, 0.5]), np.array([0.5, 0.]), np.array([0., 0.5]), np.array([1./3., 2./3.]), np.array([2./3., 1./3.])]
torus_colors = [(0., 0., 1.), (1., 0., 0.), (0., 0.5, 0.), (0., 0., 0.), (0.5, 0.5, 0.5)]   # blue, green, red, black, grey
def clmap_patterns(phase_1, phase_2): # phases in range (0, 1)
    phase = np.array([phase_1, phase_2])
    # select which pattern
    return torus_colors[
            np.argmin(
                [torus_distance(phase, torus_patterns[i])
                    for i in xrange(len(torus_patterns))])]



#==============================================

import pylab as pl
import matplotlib.patches as mpatches

xblue, yblue, xgreen, ygreen, xred, yred, width, radius = 0.23, 0.65, 0.77, 0.65, 0.47, 0.25, 0.03, 0.1
def three_cells_alt(coupling_strength, ax):

    c = coupling_strength+0.00001
    patches = []
    art = mpatches.Circle(np.array([xblue, yblue]), radius, fc='b')
    patches.append(art)
    art = mpatches.Circle(np.array([xgreen, ygreen]), radius, fc='g')
    patches.append(art)
    art = mpatches.Circle(np.array([xred, yred]), radius, fc='r')
    patches.append(art)

    sh = 0.02
    shr = 0.25
    alpha = 1.3

    # 1 blue -> green
    xbluegreen = [xblue+alpha*radius, xgreen-alpha*radius]
    ybluegreen = np.array([yblue, ygreen])+sh
    ax.add_line(pl.Line2D(xbluegreen, ybluegreen, lw=5., c='k'))
    ax.add_patch(mpatches.Circle(np.array([xbluegreen[1], ybluegreen[1]]), radius/4., fc='k'))
    ax.text((xblue+xgreen)/2.-0.05, (yblue+ygreen)/2.+0.05, repr(c[2])[:6], fontsize=15)

    # 2 green -> blue
    xgreenblue = xbluegreen
    ygreenblue = ybluegreen-2.*sh
    ax.add_line(pl.Line2D(xgreenblue, ygreenblue, lw=5., c='k'))

    ax.add_patch(mpatches.Circle(np.array([xgreenblue[0], ygreenblue[0]]), radius/4., fc='k'))

    ax.text((xblue+xgreen)/2.-0.05, (yblue+ygreen)/2.-0.1, repr(c[0])[:6], fontsize=15)

    # 3 green -> red
    xredgreen = np.array([xred+0., xgreen+0.])*0.45+0.375
    yredgreen = np.array([yred+0., ygreen+0.])*0.45+0.2185
    ax.add_line(pl.Line2D(xredgreen, yredgreen, lw=5., c='k'))
    ax.add_patch(mpatches.Circle(np.array([xredgreen[0], yredgreen[0]]), radius/4., fc='k'))
    ax.text((xgreen+xred)/2.+0.015, (ygreen+yred)/2.-0.015, repr(c[5])[:6], fontsize=15, rotation=45)

    # 4 red -> green
    xredgreen = np.array([xred+0., xgreen+0.])*0.45+0.33
    yredgreen = np.array([yred+0., ygreen+0.])*0.45+0.255
    ax.add_line(pl.Line2D(xredgreen, yredgreen, lw=5., c='k'))
    ax.add_patch(mpatches.Circle(np.array([xredgreen[1], yredgreen[1]]), radius/4., fc='k'))
    ax.text((xgreen+xred)/2.-0.1, (ygreen+yred)/2.+0.07, repr(c[3])[:6], fontsize=15, rotation=45)

    # 5 red -> blue
    xbluered = np.array([xblue+0., xred+0.])*0.45+0.17
    ybluered = np.array([yblue+0., yred+0.])*0.45+0.23
    ax.add_line(pl.Line2D(xbluered, ybluered, lw=5., c='k'))
    ax.add_patch(mpatches.Circle(np.array([xbluered[0], ybluered[0]]), radius/4., fc='k'))
    ax.text((xblue+xred)/2.-0.1, (yblue+yred)/2.-0.0, repr(c[1])[:6], fontsize=15, rotation=-55)

    # 6 blue -> red
    xbluered = np.array([xblue+0., xred+0.])*0.45+0.2
    ybluered = np.array([yblue+0., yred+0.])*0.45+0.265
    ax.add_line(pl.Line2D(xbluered, ybluered, lw=5., c='k'))
    ax.add_patch(mpatches.Circle(np.array([xbluered[1], ybluered[1]]), radius/4., fc='k'))
    ax.text((xblue+xred)/2.-0.02, (yblue+yred)/2.+0.1, repr(c[4])[:6], fontsize=15, rotation=-55)

    ax.add_patch(patches[0])
    ax.add_patch(patches[1])
    ax.add_patch(patches[2])

    ax.text(xblue-0.04, yblue-0.05, r'$1$', fontsize=40)
    ax.text(xgreen-0.04, ygreen-0.05, r'$2$', fontsize=40)
    ax.text(xred-0.04, yred-0.05, r'$3$', fontsize=40)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()













if __name__ == '__main__':

    K = 0.01*np.ones((6), float)
    K[2] = 0.02 		# (1) -> (2)
    K[5] = 0.02		# (2) -> (3)
    K[1] = 0.02		# (3) -> (1)

    fig = pl.figure()
    ax = fig.add_axes([0., 0., 1., 1.])
    three_cells_alt(K, ax=ax)
    fig.show()
    time.sleep(1)









