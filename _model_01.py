#!/usr/bin/env python3

#
# Test model: 1D along x-direction
#
# Author: Chuck Jia
# Date: November 7, 2019
#


from boundaries import ZeroDirichletBoundary
from godunov import PhysicalFlux
from mesh import Mesh
from parameters import Parameters
import numpy as np
from numpy import (
    cos,
    exp,
    pi,
    sin,
)
from sources import Sources
from velocities import UW


num_space_steps = 200
Nt = 100

Nx, Nz = num_space_steps, num_space_steps
tf = 1.0
Dt = tf / Nt

x0 = 0.0
xf = 5e4
zf = 1000.0


def z0_fcn(x):
    return 200


def z0_der_fcn(x):
    return 0


class UWT(UW):
    def __init__(self, mesh):
        self.c_x = 2 * pi / mesh.param.xf
        self.c_z = pi / 1e3

        super().__init__(mesh)
        self.u_ct_m, self.w_ct_m = self.uw_fcn(self.mesh.x_ct_m, self.mesh.z_ct_m)

    def uw_fcn(self, x, z):
        cx_x = self.c_x * x
        cz_z = self.c_z * z

        u = 7.5 + cos(cx_x) * cos(cz_z)
        w = cos(cx_x) * sin(cz_z)

        return u, w

    def uw_der_fcn(self, x, z):
        c1 = self.c_x
        c2 = self.c_z
        c1_x = self.c_x * x
        c2_z = self.c_z * z

        du_dx = -c1 * sin(c1_x) * cos(c2_z)
        dw_dz = c2 * cos(c1_x) * cos(c2_z)

        return du_dx, dw_dz


class SourcesT(Sources):
    two_pi = 2.0 * pi

    def __init__(self, mesh, uw):
        super().__init__(mesh, uw)

        self.c_x = 2 * self.two_pi / (self.xf - self.x0)
        self.c_z = pi / 200.0

        self.exact_solution_fcn_list = (self.exact_solution,) * 4

    def exact_solution(self, x, z, t):
        return cos(self.two_pi * t) * sin(self.c_x * x) * sin(self.c_z * z)

    def source_fcn_pointwise(self, x, z, t):
        twopi_t = self.two_pi * t
        c_x, c_z = self.c_x, self.c_z
        cx_x, cz_z = self.c_x * x, self.c_z * z

        cos_twopi_t = cos(twopi_t)

        sol = self.exact_solution(x, z, t)
        dsol_dx = cos_twopi_t * c_x * cos(cx_x) * sin(cz_z)
        dsol_dz = cos_twopi_t * c_z * sin(cx_x) * cos(cz_z)

        u, w = self.uw.u_ct_m, self.uw.w_ct_m
        du_dx, dw_dz = self.uw.uw_der_fcn(x, z)

        source_val = -self.two_pi * sin(twopi_t) * sin(cx_x) * sin(cz_z) + \
            sol * (du_dx + dw_dz) + u * dsol_dx + w * dsol_dz

        source_tuple = (source_val,) * 4

        return (source_tuple, 0)

    def source_fcn(self, t):
        return self.source_fcn_pointwise(x=self.x_ct_m, z=self.z_ct_m, t=t)


param = Parameters(
    x0=x0,
    xf=xf,
    z0_fcn=z0_fcn,
    z0_der_fcn=z0_der_fcn,
    zf=zf,
    Nx=Nx,
    Nz=Nz,
    tf=tf,
    Nt=Nt,
    Mesh=Mesh,
    UW=UWT,
    Flux=PhysicalFlux,
    Boundary=ZeroDirichletBoundary,
    Sources=SourcesT,
    init_fcn_list=None,
    rho_o_fcn=None,
    time_method="rk4",
    num_msg=20,
    num_csv=1,
    l2_error_period=-1,
    ave_period=-1,
)


if __name__ == "__main__":
    mesh = Mesh(param=param)
    st = SourcesT(mesh=mesh)
