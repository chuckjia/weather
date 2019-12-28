from numpy import exp
import numpy as np
import physics as phys
import utilities as util


class InitialCondition:
    def init_theta_fcn(self, x, z):
        return

    def init_qv_fcn(self, x, z):
        return

    def init_qc_fcn(self, x, z):
        return

    def init_qp_fcn(self, x, z):
        return

    def get_init_fcn_list(self):
        return (
            self.init_theta_fcn,
            self.init_qv_fcn,
            self.init_qc_fcn,
            self.init_qp_fcn,
        )



class PhysicalInitialCondition(InitialCondition):
    T0 = phys.T0

    def init_theta_fcn(self, x, z):
        Delta_T = 50
        p_e = phys.p_e_fcn(z)
        init_T = self.T0 - (1 - p_e / phys.p0) * Delta_T
        init_theta = phys.T_to_theta(init_T, p_e)

        return init_theta

    def init_qv_fcn(self, x, z):
        theta = self.init_theta_fcn(x, z)
        p_e = phys.p_e_fcn(z)
        T = phys.theta_to_T(theta, p_e)
        saturation = 0.5

        return saturation * 3.801664 * exp(17.67 * (T - 273.15) / (T - 29.65))

    def init_qc_fcn(self, x, z):
        return np.zeros(x.shape)

    def init_qp_fcn(self, x, z):
        return np.zeros(x.shape)
