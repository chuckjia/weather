#!/usr/bin/env python3
#
# Physical model
#
# Author: Chuck Jia
# Date: November 7, 2019
#


from boundaries import PhysicalBoundary
from godunov import PhysicalFlux
from initialcond import PhysicalInitialCondition as PIC
from mesh import Mesh
from parameters import Parameters
import physics as phys
import numpy as np
from numpy import (
    cos,
    exp,
    pi,
    sin,
)
from sources import PhysicalSources
from velocities import PhysicalUW


num_space_steps = 100
Nt = 1000
Dt = 0.5

x0 = 0.0
xf = 9e4
zf = 1.6e4


def z0_fcn(x):
    # height = 2500.0
    # steepness = 6000.0
    # x_center = 0.5 * (xf + x0)
    #
    # return height * exp(- ((x - x_center) / steepness) ** 2)

    return 0


def z0_der_fcn(x):
    return 0


param = Parameters(
    x0=x0,
    xf=xf,
    z0_fcn=z0_fcn,
    z0_der_fcn=z0_der_fcn,
    zf=zf,
    Nx=num_space_steps,
    Nz=num_space_steps,
    tf=Dt * Nt,
    Nt=Nt,
    Mesh=Mesh,
    UW=PhysicalUW,
    Flux=PhysicalFlux,
    Boundary=PhysicalBoundary,
    Sources=PhysicalSources,
    init_fcn_list=PIC().get_init_fcn_list(),
    rho_o_fcn=phys.rho_o_fcn,
    time_method="rk4",
    num_msg=200,
    num_csv=200,
    l2_error_period=-1,
    ave_period_list=[10, 10, -1, -1],
)
