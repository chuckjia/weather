import numpy as np
import os
import physics as phys
import time
import utilities as util

class Mesh():
    def __init__(self, param):
        start_time = time.time()

        self.param = param
        self.get_params = param.get_space_params

        x0, xf, Nx, Nz = param.x0, param.xf, param.Nx, param.Nz
        self.ave_period_list = self.param.ave_period_list
        self.Dx = (xf - x0) / Nx
        self.grid_shape = (Nx + 2, Nz + 2)

        self.initialize_grid_coord_matrices()
        self.initialize_cell_center_coord_matrices()
        self.initialize_solution_matrices()
        self.initialize_cell_side_matrices()
        self.initialize_cell_vol_matrix()
        self.initialize_rho_o()

        self.solution_setter_tuple = (
            self.set_theta, self.set_qv, self.set_qc, self.set_qp
        )

        print(
            "Initialized mesh grid. "
            "Time used: %1.2fs." % (time.time() - start_time)
        )

    def get_x_z_matrices(self):
        return (
            self.x_bl_m, self.x_br_m, self.x_tl_m, self.x_tr_m,
            self.z_bl_m, self.z_br_m, self.z_tl_m, self.z_tr_m,
        )

    def get_x_matrices(self):
        return self.x_bl_m, self.x_br_m, self.x_tl_m, self.x_tr_m

    def get_z_matrices(self):
        return self.z_bl_m, self.z_br_m, self.z_tl_m, self.z_tr_m

    def get_center_points_matrices(self):
        return self.x_ct_m, self.z_ct_m

    def initialize_solution_matrices(self):
        grid_shape = self.grid_shape
        self.theta_m = np.zeros(grid_shape)
        self.qv_m = np.zeros(grid_shape)
        self.qc_m = np.zeros(grid_shape)
        self.qp_m = np.zeros(grid_shape)

    def get_solution_matrices(self):
        return self.theta_m, self.qv_m, self.qc_m, self.qp_m

    def set_theta(self, theta_m):
        self.theta_m = theta_m

    def set_qv(self, qv_m):
        self.qv_m = qv_m

    def set_qc(self, qc_m):
        self.qc_m = qc_m

    def set_qp(self, qp_m):
        self.qp_m = qp_m

    def set_solutions(self, theta_m, qv_m, qc_m, qp_m):
        self.theta_m = theta_m
        self.qv_m = qv_m
        self.qc_m = qc_m
        self.qp_m = qp_m

    def initialize_grid_coord_matrices(self):
        x0, xf, z0_fcn, zf, Nx, Nz = self.get_params()

        def get_z_lists_for_x(x):
            z0 = z0_fcn(x)
            Dz = (zf - z0) / Nz
            return (
                [z0 + Dz * (j - 1) for j in range(Nz + 2)],  # Cell bottom
                [z0 + Dz * j for j in range(Nz + 2)],  # Cell top
            )

        # Aliases for convenience
        Dx = (xf - x0) / Nx
        grid_shape = (Nx + 2, Nz + 2)

        # 2D arrays for x-coords: the trailing "_m" means it is a matrix
        self.x_bl_m = np.zeros(grid_shape)
        self.x_br_m = np.zeros(grid_shape)
        self.x_tl_m = self.x_bl_m
        self.x_tr_m = self.x_br_m

        # 2D arrays for z-coords
        self.z_bl_m = np.zeros(grid_shape)
        self.z_br_m = np.zeros(grid_shape)
        self.z_tl_m = np.zeros(grid_shape)
        self.z_tr_m = np.zeros(grid_shape)

        # Compute x coords at cell corners
        for i in range(Nx + 2):
            x = x0 + Dx * (i - 1)

            self.x_bl_m[i, :] = x
            self.x_br_m[i, :] = x + Dx

        # Compute z coords at cell corners
        for i in range(Nx + 2):
            x_left = x0 + Dx * (i - 1)
            self.z_bl_m[i, :], self.z_tl_m[i, :] = get_z_lists_for_x(x_left)

            x_right = x0 + Dx * i
            self.z_br_m[i, :], self.z_tr_m[i, :] = get_z_lists_for_x(x_right)

    def initialize_cell_center_coord_matrices(self):
        x1_m, x2_m, x4_m, x3_m, z1_m, z2_m, z4_m, z3_m = self.get_x_z_matrices()

        x_124_m = (x1_m + x2_m + x4_m) / 3.0
        x_123_m = (x1_m + x2_m + x3_m) / 3.0
        x_234_m = (x2_m + x3_m + x4_m) / 3.0
        x_134_m = (x1_m + x3_m + x4_m) / 3.0

        z_124_m = (z1_m + z2_m + z4_m) / 3.0
        z_123_m = (z1_m + z2_m + z3_m) / 3.0
        z_234_m = (z2_m + z3_m + z4_m) / 3.0
        z_134_m = (z1_m + z3_m + z4_m) / 3.0

        k_l_m = (z_134_m - z_123_m) / (x_134_m - x_123_m)
        k_r_m = (z_234_m - z_124_m) / (x_234_m - x_124_m)

        self.x_ct_m = (
            k_l_m * x_123_m - k_r_m * x_124_m - z_123_m + z_124_m
        ) / (k_l_m - k_r_m)

        self.z_ct_m = k_l_m * (self.x_ct_m - x_123_m) + z_123_m

    def initialize_cell_side_matrices(self):
        def compute_horizontal_sides(x_l_m, z_l_m, x_r_m, z_r_m):
            side_vec_x_m = x_r_m - x_l_m
            side_vec_z_m = z_r_m - z_l_m
            side_len_m = (side_vec_x_m ** 2 + side_vec_z_m ** 2) ** 0.5
            n_side_x_m = -side_vec_z_m / side_len_m
            n_side_z_m = side_vec_x_m / side_len_m

            return side_len_m, n_side_x_m, n_side_z_m

        self.top_side_len_m, self.n_top_x_m, self.n_top_z_m = \
            compute_horizontal_sides(
                self.x_tl_m, self.z_tl_m, self.x_tr_m, self.z_tr_m
            )
        self.bott_side_len_m, self.n_bott_x_m, self.n_bott_z_m = \
            compute_horizontal_sides(
                self.x_bl_m, self.z_bl_m, self.x_br_m, self.z_br_m
            )

        self.left_side_len_m = self.z_tl_m - self.z_bl_m
        self.right_side_len_m = self.z_tr_m - self.z_br_m

    def get_cell_top_side_len_matrix(self):
        return self.top_side_len_m

    def get_cell_left_side_len_matrix(self):
        return self.left_side_len_m

    def get_cell_right_side_len_matrix(self):
        return self.right_side_len_m

    def get_cell_top_normal_vector_matrices(self):
        return self.n_top_x_m, self.n_top_z_m

    def initialize_cell_vol_matrix(self):
        self.cell_vol_m = \
            (self.left_side_len_m + self.right_side_len_m) * (self.Dx * 0.5)

    def get_cell_vol_matrix(self):
        return self.cell_vol_m

    def apply_initial_conditions(self, init_fcn_list):
        solution_setters = self.solution_setter_tuple

        if len(init_fcn_list) != len(solution_setters):
            raise Exception("Wrong number of initial conditions!")

        num_init_fcn_args = util.get_function_arg_number(init_fcn_list[0])

        if num_init_fcn_args == 2:
            init_fcn_args = (self.x_ct_m, self.z_ct_m)
        elif num_init_fcn_args == 3:
            init_fcn_args = (self.x_ct_m, self.z_ct_m, 0)
        else:
            raise Exception(
                "Incorrect number of input arguments for initial condition "
                "functions! Only 2 arguments (x, z) or 3 arguments (x, z, t) "
                "are allowed."
            )

        for init_fcn in init_fcn_list:
            if util.get_function_arg_number(init_fcn) != num_init_fcn_args:
                raise Exception(
                    "Inconsistent number of input arguments in the initial "
                    "condition functions!"
                )

        rho_o_m = self.rho_o_m
        if rho_o_m is None:
            for set_sol, init_fcn in zip(solution_setters, init_fcn_list):
                set_sol(init_fcn(*init_fcn_args))
        else:
            for set_sol, init_fcn in zip(solution_setters, init_fcn_list):
                set_sol(rho_o_m * init_fcn(*init_fcn_args))

    def initialize_rho_o(self):
        if self.param.rho_o_fcn is None:
            self.rho_o_m = None
            return

        self.rho_o_m = self.param.rho_o_fcn(z=self.z_ct_m)

    def average_solutions(self, step_no, force=False):
        if self.ave_period_list is None:
            return

        for sol_m, set_sol, ave_period in zip(
            self.get_solution_matrices(),
            self.solution_setter_tuple,
            self.ave_period_list,
        ):
            if force or step_no % ave_period == 0:
                new_sol_m = 0.5 * (
                    sol_m + util.get_M_i_minus_one(M=sol_m, padding="mirror")
                )
                set_sol(new_sol_m)

    def save_matrices_to_csv(
        self, matrix_list, filename_list,
        results_folder_name="", filename_suffix="",
    ):
        for matrix, filename in zip(matrix_list, filename_list):
            filepath = os.path.join(
                results_folder_name, filename + filename_suffix + ".csv"
            )
            np.savetxt(filepath, matrix, delimiter=',')

    def save_coord_matrices_to_csv(self):
        filename_list = [
            "x_bl_m", "x_br_m", "x_tl_m", "x_tr_m",
            "z_bl_m", "z_br_m", "z_tl_m", "z_tr_m",
            "x_ct_m", "z_ct_m",
            "cell_vol_m",
        ]
        matrix_list = [
            self.x_bl_m, self.x_br_m, self.x_tl_m, self.x_tr_m,
            self.z_bl_m, self.z_br_m, self.z_tl_m, self.z_tr_m,
            self.x_ct_m, self.z_ct_m,
            self.cell_vol_m,
        ]

        self.save_matrices_to_csv(
            matrix_list=matrix_list,
            filename_list=filename_list,
            results_folder_name="results",
        )

    def save_solution_matrices_to_csv(self, t=None, theta_to_T=False):
        filename_list = ["theta", "qv", "qc", "qp"]

        rho_o_m = self.rho_o_m
        theta_m = self.theta_m
        qv_m = self.qv_m
        qc_m = self.qc_m
        qp_m = self.qp_m

        if rho_o_m is not None:
            theta_m = theta_m / rho_o_m
            qv_m = qv_m / rho_o_m
            qc_m = qc_m / rho_o_m
            qp_m = qp_m / rho_o_m

        if theta_to_T:
            p_e = phys.p_e_fcn(self.z_ct_m)
            theta_m = phys.theta_to_T(theta_m, p_e)

        matrix_list = [theta_m, qv_m, qc_m, qp_m]

        if t is None:
            suffix = ""
        else:
            suffix = "_" + str(t)

        self.save_matrices_to_csv(
            matrix_list=matrix_list,
            filename_list=filename_list,
            results_folder_name="solutions",
            filename_suffix=suffix,
        )


if __name__ == "__main__":
    pass
