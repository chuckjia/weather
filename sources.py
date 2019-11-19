import numpy as np
from numpy import (
    exp,
    log,
    pi,
)
import physics as phys
import utilities as util

class Sources():
    def __init__(self, mesh, uw):
        self.param = mesh.param
        self.mesh = mesh
        self.uw = uw

        # Aliases
        self.get_solutions = mesh.get_solution_matrices
        self.rho_o_m = mesh.rho_o_m
        self.x0, self.xf, self.z0_fcn, self.zf = \
            self.param.x0, self.param.xf, self.param.z0_fcn, self.param.zf
        self.x_ct_m, self.z_ct_m = mesh.x_ct_m, mesh.z_ct_m

        # Exact solution for tests
        self.exact_solution_fcn_list = None

    def source_fcn(self, t):
        return None

    def create_ghost_mask(self):
        """A mask in which ghosts cells are 0 and real cells are 1. It can be
           used to mask the ghost cells duing calculation of the sources."""

        param = self.param

        mask = np.ones((param.Nx + 2, param.Nz + 2))
        mask[0, :] = 0
        mask[-1, :] = 0
        mask[:, 0] = 0
        mask[:, -1] = 0

        return mask


class PhysicalSources(Sources):
    def __init__(self, mesh, uw):
        super().__init__(mesh, uw)
        self.ghost_mask = self.create_ghost_mask()

    def alpha_fcn(self, T):
        """
        A continuous function on temperature.

        Args:
            T : A float or a matrix. The temperature, or the matrix of
                temperatures in the entire mesh

        Returns:
            alpha, one_minus_alpha : The alpha and (1 - alpha) value, or the
                                     corresponding matrix for the entire mesh
        """

        T_w = 273.15
        T_i = 263.15

        alpha = (T >= T_i) * (T <= T_w) * (T - T_i) / (T_w - T_i) + (T > T_w)
        one_minus_alpha = 1 - alpha

        return alpha * self.ghost_mask, one_minus_alpha * self.ghost_mask

    e_oo = phys.e_oo
    L_v = phys.L_v
    L_s = phys.L_s
    R_v = phys.R_v

    def e_s_fcn(self, T):
        """Returns pair e_sw, e_si, the saturated water vapor pressures over
           water and over ice, resp. Defined in Eq.201910201631."""

        T_oo = 273.16
        common_part = (1.0 / T_oo - 1.0 / T) / self.R_v

        e_sw = self.e_oo * exp(self.L_v * common_part)
        e_si = self.e_oo * exp(self.L_s * common_part)

        return e_sw, e_si

    epsilon = 0.622

    def q_vs_fcn(self, T, z, p_e, e_s_pair):
        """The saturated water vapor mixing ratios."""

        e_sw, e_si = e_s_pair
        epsilon = self.epsilon

        q_vsw = epsilon * e_sw / (p_e - e_sw)
        q_vsi = epsilon * e_si / (p_e - e_si)

        return q_vsw, q_vsi

    c_p = phys.c_p
    Delta_t = phys.Delta_t

    c1_con = Delta_t * L_v * L_v / (c_p * R_v)
    c2_con = Delta_t * L_s * L_s / (c_p * R_v)

    def con_source(self, T, r_qv, rho_o, q_vs_pair, alpha_pair):
        """CON: Cloud bulk condensation rate from water vapor."""

        alpha, one_minus_alpha = alpha_pair
        q_vsw, q_vsi = q_vs_pair

        T_squared = T ** 2
        con_w = (r_qv - rho_o * q_vsw) / (
            self.Delta_t + self.c1_con / T_squared * q_vsw
        )
        con_i = (r_qv - rho_o * q_vsi) / (
            self.Delta_t + self.c2_con / T_squared * q_vsi
        )

        return alpha * con_w + one_minus_alpha * con_i

    N_d = phys.N_d
    D_d = phys.D_d
    c_aut = 2.16e3 * N_d / D_d

    def aut_source(self, r_qc, T, z, rho_o, alpha_pair):
        alpha, one_minus_alpha = alpha_pair

        psi = 1e3 * alpha * r_qc
        aut_r = (psi ** 3) / (3e5 * psi + self.c_aut)

        tau_a = -800.0 * exp(-((T + 15.0) ** 2)) + 1000.0
        aut_s = one_minus_alpha * r_qc / tau_a

        return aut_r + aut_s

    N0_D_bar = 1e7  # Fixed parameter for the Marshall-Palmer size distribution
    rho_w = phys.rho_w

    c1_D_r_bar = ((pi * 1e10) ** (1.0 / 12.0)) / (
        (pi / 6.0 * rho_w * N0_D_bar) ** (1.0 / 3.0)
    )
    c2_D_s_bar = ((40.0 / N0_D_bar) ** 0.5) * (5e5 ** (1.0 / 3.0))

    def D_bar_fcn(self, r_qp, alpha_pair):
        alpha, one_minus_alpha = alpha_pair

        D_r_bar = self.c1_D_r_bar * (alpha * r_qp) ** 0.25
        D_s_bar = self.c2_D_s_bar * ((one_minus_alpha * r_qp) ** (1.0 / 6.0))

        return D_r_bar, D_s_bar

    c1_acc_r = 26 * pi
    c2_acc_s = 0.06 * pi

    def acc_source(self, r_qc, D_bar_pair, alpha_pair):
        alpha, one_minus_alpha = alpha_pair
        D_r_bar, D_s_bar = D_bar_pair

        acc_r = self.c1_acc_r * alpha * (D_r_bar ** 2.5) * r_qc
        acc_s = self.c2_acc_s * one_minus_alpha * (D_s_bar ** 2.25) * r_qc

        return acc_r + acc_s

    c1_dep_r = 2 * pi
    c2_dep_s = 4.0 / 3.0 * pi

    def G_fcn(self, T_e, e_s_pair):
        e_sw, e_si = e_s_pair
        e_s = e_sw * (T_e >= 273.15) + e_si * (T_e < 273.15)
        G = (1 / 2.2e7) / (T_e / e_s + 100.0 / T_e)

        return G

    def dep_source(self, qv, T_e, D_bar_pair, q_vs_pair, e_s_pair):
        D_r_bar, D_s_bar = D_bar_pair
        q_vsw, q_vsi = q_vs_pair

        F_r = 0.78 + 6.88367634335026082226e2 * (D_r_bar ** 0.75)
        F_s = 0.65 + 1.74413302244983611899e2 * (D_s_bar ** 0.625)

        G = self.G_fcn(T_e=T_e, e_s_pair=e_s_pair)

        dep_r = self.c1_dep_r * D_r_bar * (qv / q_vsw - 1) * F_r * G
        dep_s = self.c2_dep_s * D_s_bar * (qv / q_vsi - 1) * F_s * G

        return dep_r + dep_s

    def V_Ts_fcn(self, r_qp, one_minus_alpha):
        return 1.70816352295439766 * (one_minus_alpha * r_qp) ** (1.0 / 12.0)

    c_source = L_v / c_p

    def source_fcn_pointwise(self, r_theta, r_qv, r_qc, r_qp, x, z, rho_o):
        """Compute and return the sources. This function is vectorized. It can
           be applied both pointwise and on numpy matrices."""

        theta_e = phys.theta_e_fcn(z=z)
        p_e = phys.p_e_fcn(z=z)

        theta = r_theta / rho_o
        qv = r_qv / rho_o

        T_e = phys.theta_to_T(theta=theta_e, p_e=p_e)
        T = phys.theta_to_T(theta=theta, p_e=p_e)
        alpha_pair = self.alpha_fcn(T=T)

        e_s_pair = self.e_s_fcn(T=T)
        q_vs_pair = self.q_vs_fcn(T=T, z=z, p_e=p_e, e_s_pair=e_s_pair)
        D_bar_pair = self.D_bar_fcn(r_qp=r_qp, alpha_pair=alpha_pair)

        con = self.con_source(
            T=T,
            r_qv=r_qv,
            rho_o=rho_o,
            q_vs_pair=q_vs_pair,
            alpha_pair=alpha_pair,
        )
        acc = self.acc_source(
            r_qc=r_qc, D_bar_pair=D_bar_pair, alpha_pair=alpha_pair,
        )
        aut = self.aut_source(
            r_qc=r_qc, T=T, z=z, rho_o=rho_o, alpha_pair=alpha_pair,
        )
        dep = self.dep_source(
            qv=qv,
            T_e=T_e,
            D_bar_pair=D_bar_pair,
            q_vs_pair=q_vs_pair,
            e_s_pair=e_s_pair,
        )

        V_Ts = self.V_Ts_fcn(r_qp=r_qp, one_minus_alpha=alpha_pair[1])

        source_tuple = (
            theta_e / T_e * self.c_source * (con + dep),
            -con - dep,
            con - acc - aut,
            acc + aut + dep,
        )

        return source_tuple, V_Ts

    def source_fcn(self, t):
        r_theta_m, r_qv_m, r_qc_m, r_qp_m = self.get_solutions()

        return self.source_fcn_pointwise(
            r_theta=r_theta_m,
            r_qv=r_qv_m,
            r_qc=r_qc_m,
            r_qp=r_qp_m,
            x=self.x_ct_m,
            z=self.z_ct_m,
            rho_o=self.rho_o_m,
        )
