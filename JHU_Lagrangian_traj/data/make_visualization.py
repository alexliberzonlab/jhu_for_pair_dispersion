#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  3 11:32:09 2018

@author: ron
"""


from flowtracks.io import Scene
from flowtracks.trajectory import Trajectory
from mayavi import mlab
import numpy as np


S = Scene('jhu_traj_t4_n150.h5')

t = []

for tr in S.iter_trajectories():
    t.append(tr)
    
tt = t[400:800:4]



def get_a(tr):
    return np.linalg.norm(tr.accel()[0,:])
a = [get_a(tr) for tr in tt]
A = np.amax(a)


def get_x(tr):
    return np.linalg.norm(tr.pos()[0,:] - tr.pos()[-1,:])
x = [get_x(tr) for tr in tt]
X = np.amax(x)

for tr in tt[:50]:
    c = ( float(get_x(tr)/X ), 0.2, 1.0 - float(get_x(tr)/X ))
    mlab.plot3d(tr.pos()[:,0], tr.pos()[:,1], tr.pos()[:,2], color = c,
                tube_radius = 0.01)