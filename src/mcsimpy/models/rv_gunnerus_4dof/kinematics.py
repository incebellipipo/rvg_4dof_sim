#!/usr/bin/env python3

# MCSim_python
# Copyright (C) 2022, NTNU - Norges teknisk-naturvitenskapelige universitet
# This file is part of MCSim_python.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#
# Created By: M. Marley
# Created Date: 2022-02-04


# Log:
# - Tested:  2022-02-04 M.Marley input/output relations for individual functions
#

"""
Library of kinematics functions
"""

import numpy as np

# =============================================================================
# Functions
# =============================================================================
def rad2pipi(angle):
    """Maps angle to [-pi,pi]
       Input
           angle: in radians
       Output
           angle: in [-pi,pi]
    """
    # Created: 2022-02-04 M.Marley
    # Tested: 2022-02-04 M.Marley
    angle = np.arctan2(np.sin(angle),np.cos(angle))
    return angle

def Rotx(angle):
    """Rotation matrix about x-axis
       Input
           angle: in radians
       Output
           Rx: 3x3 matrix
    """
    # Created: 2022-02-04 M.Marley
    # Tested: 2022-02-04 M.Marley
    Rx = np.array([[1, 0, 0],
          [0, np.cos(angle), -np.sin(angle)],
          [0, np.sin(angle), np.cos(angle)]])
    return Rx

def Roty(angle):
    """Rotation matrix about y-axis
       Input
           angle: in radians
       Output
           Ry: 3x3 matrix
    """
    # Created: 2022-02-04 M.Marley
    # Tested: 2022-02-04 M.Marley
    Ry = np.array([[np.cos(angle),0, np.sin(angle)],
          [0,1,0],
          [-np.sin(angle), 0, np.cos(angle)]])
    return Ry

def Rotz(angle):
    """Rotation matrix about z-axis
       Input
           angle: in radians
       Output
           Rz: 3x3 matrix
    """
    # Created: 2022-02-04 M.Marley
    # Tested: 2022-02-04 M.Marley
    Rz = np.array([[np.cos(angle), -np.sin(angle), 0],
          [np.sin(angle), np.cos(angle), 0],
          [0, 0, 1]])
    return Rz


def Rotzyx(angles):
    """Rotation matrix about zyx-axis
       Input
           angles: euler angles in radians
       Output
           R: 3x3 matrix
    """
    # Created: 2022-05-04 M.Marley
    # Tested:
    Rzyx = Rotz(angles[2])@Roty(angles[1])@Rotx(angles[0])
    return Rzyx


def Smat(x):
    """Skew-symmetric cross-product operator matrix
       Input
           x: 3x1 vector
       Output
           S: 3x3 matrix
    """
    # Created: 2022-02-04 M.Marley
    # Tested: 2022-02-04 M.Marley
    S = np.array([[0, -x[2], x[1]],
          [x[2], 0, -x[0]],
          [-x[1], x[0], 0]])
    return S

def dot_eta3(psi,nu):
    """Derivative of 3DOF pose vector p=[x,y,psi]
       Input
           psi: heading
           nu: 3DOF body-fixed velocities
       Output
           d_eta3: 3x1 vector
    """
    # Created: 2022-02-04 M.Marley
    # Tested: 2022-02-04 M.Marley
    d_eta3 = Rotz(psi)@nu
    return d_eta3

def dot_eta6(eta,nu):
    """Derivative of 6DOF position and orientation vector
       Input
           eta: position and orientation
           nu: 6DOF body-fixed velocities
       Output
           d_eta6: 6x1 vector
    """
    # Created: 2022-05-02 M.Marley
    # Tested: 2022-05-02 M.Marley

    Rzyx = Rotzyx(eta[3:6])

    phi = eta[3]
    theta = eta[4]
#    psi = eta[5]
    T = np.array([[1,np.sin(phi)*np.tan(theta),np.cos(phi)*np.tan(theta)],
                  [0, np.cos(phi), -np.sin(phi)],
                  [0, np.sin(phi)/np.cos(theta), np.cos(phi)/np.cos(theta)]])

    d_eta6 = np.zeros(len(eta))

    d_eta6[0:3] = Rzyx@nu[0:3]
    d_eta6[3:6] = T@nu[3:6]

    return d_eta6

def Rot_planar(angle):
    """Planar rotation matrix R: \R -> SO(2)
       Input
           angle: in radians
       Output
           R: 2x2 matrix
    """
    # Created: 2022-11-15 M.Marley
    # Tested: 2022-11-15 M.Marley
    R = np.array([[np.cos(angle), -np.sin(angle)],
          [np.sin(angle), np.cos(angle)]])
    return R

def dot_eta4(psi,nu):
    """4DOF kinematics
       Input
           angle: in radians
       Output
           d_eta4: derivatives
    """
    # Created: 2022-11-15 M.Marley
    # Tested: 2022-11-15 M.Marley

    R = np.array([[np.cos(psi), -np.sin(psi),0,0],
          [np.sin(psi), np.cos(psi),0,0],
          [0,0,1,0],[0,0,0,1]])

    d_eta4 = R@nu

    return d_eta4
