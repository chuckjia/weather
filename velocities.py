import numpy as np
from numpy import (
    sin,
    cos,
    maximum as max,
    minimum as min,
)
import time


class UW():
    def __init__(self, mesh):
        self.param = mesh.param
        self.mesh = mesh
        self.save_matrices_to_csv = mesh.save_matrices_to_csv

        x_bl_m, x_br_m, x_tl_m, x_tr_m = mesh.get_x_matrices()
        z_bl_m, z_br_m, z_tl_m, z_tr_m = mesh.get_z_matrices()
        x_ct_m, z_ct_m = mesh.get_center_points_matrices()

        self.x_l_side_m = x_bl_m
        self.x_r_side_m = x_br_m
        self.x_b_side_m = x_ct_m
        self.x_t_side_m = x_ct_m

        self.z_l_side_m = 0.5 * (z_bl_m + z_tl_m)
        self.z_r_side_m = 0.5 * (z_br_m + z_tr_m)
        self.z_b_side_m = 0.5 * (z_bl_m + z_br_m)
        self.z_t_side_m = 0.5 * (z_tl_m + z_tr_m)

        self.set_uw_matrices(uw_fcn=self.uw_fcn)

    def uw_fcn(self, x, z):
        return None, None

    def set_uw_matrices(self, uw_fcn):
        uw_list = []
        for x_side_m, z_side_m in self.get_xz_side_matrix_pairs():
            uw_list.extend(uw_fcn(x_side_m, z_side_m))

        self.u_l_m, self.w_l_m, \
        self.u_r_m, self.w_r_m, \
        self.u_b_m, self.w_b_m, \
        self.u_t_m, self.w_t_m = uw_list

    def get_x_side_matrices(self):
        return self.x_l_side_m, self.x_r_side_m, self.x_b_side_m, self.x_t_side_m

    def get_z_side_matrices(self):
        return self.z_l_side_m, self.z_r_side_m, self.z_b_side_m, self.z_t_side_m

    def get_xz_side_matrix_pairs(self):
        return (
            (self.x_l_side_m, self.z_l_side_m),
            (self.x_r_side_m, self.z_r_side_m),
            (self.x_b_side_m, self.z_b_side_m),
            (self.x_t_side_m, self.z_t_side_m),
        )

    def get_u_matrices(self):
        return self.u_l_m, self.u_r_m, self.u_b_m, self.u_t_m

    def get_w_matrices(self):
        return self.w_l_m, self.w_r_m, self.w_b_m, self.w_t_m

    def projection(self):
        return

    def save_uw_to_csv(self):
        filename_list = [
            "u_l_m", "u_r_m", "u_b_m", "u_t_m", "w_l_m", "w_r_m", "w_b_m", "w_t_m",
        ]
        matrix_list = [
            self.u_l_m, self.u_r_m, self.u_b_m, self.u_t_m,
            self.w_l_m, self.w_r_m, self.w_b_m, self.w_t_m,
        ]
        self.save_matrices_to_csv(
            matrix_list=matrix_list,
            filename_list=filename_list,
            results_folder_name="results",
        )


class PhysicalUW(UW):
    def __init__(self, mesh):
        start_time = time.time()

        super().__init__(mesh=mesh)
        self.set_uw_matrices(uw_fcn=self.uw_fcn)

        print(
            "Initialized velocities. "
            "Time used: %1.2fs." % (time.time() - start_time)
        )

    def uw_fcn(self, x, z):
        # X = 10e3
        # Z = 15e3
        # A = 4.8e4
        # S = 2.5e-2
        # x_c = 30e3
        # pi = np.pi
        # rho_oo = 1.225
        #
        # def rho_o_fcn(z):
        #     return -1e-5 * rho_oo * z + rho_oo
        #
        # drho_o_dz = -1e-5 * rho_oo
        # x_hat = max(-X, min(X, x - x_c))
        # z_hat = min(z, Z)
        # c_x = pi / X
        # c_z = pi / Z
        # rho_o = rho_o_fcn(z)
        #
        # # Compute u
        # neg_dpsi_dz_large_z = S * z
        # c_z_times_z = c_z * z
        # neg_dpsi_dz_small_z = (-A / rho_oo) * sin(c_x * x_hat) * (
        #     drho_o_dz * sin(c_z_times_z) + c_z * rho_o * cos(c_z_times_z)
        # ) + neg_dpsi_dz_large_z
        #
        # u_m = (
        #     neg_dpsi_dz_small_z * (z < Z) + neg_dpsi_dz_large_z * (z >= Z)
        # ) / rho_o
        #
        # # Compute w
        # neg_dpsi_dx = (-A / rho_oo * c_x) * sin(c_z * z_hat) * cos(c_x * (x - x_c))
        # w_m = neg_dpsi_dx * (x > (x_c - X)) * (x < (x_c + X)) / rho_o

        # Right flowing u and w

        def normalize(x, y, norm=1):
            vec_size = (x ** 2 + y ** 2) ** 0.5 / norm
            return x / vec_size, y / vec_size

        u_m = np.ones(x.shape)
        w_m = np.zeros(x.shape)

        for i in range(1, x.shape[0] - 1):
            # c = cos(2 * pi * x[i, 0])
            for j in range(1, x.shape[1] - 1):
                x_vec = x[i + 1, j] - x[i, j]
                z_vec = z[i + 1, j] - z[i, j]
                u_m[i, j], w_m[i, j] = normalize(x=x_vec, y=z_vec, norm=0.5)

        return u_m, w_m


if __name__ == "__main__":
    pass
