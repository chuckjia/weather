import numpy as np
import utilities as util


class TimeStepping():
    def __init__(self, param):
        # Initialize objects using param
        self.param = param
        self.mesh, self.uw, self.flux, self.boundary, self.sources = \
            param.initialize_objects()

        # Aliasing functions for convenience
        self.apply_boundary_conditions = self.boundary.apply_boundary_conditions
        self.average_solutions = self.mesh.average_solutions
        self.compute_all_fluxes = self.flux.compute_all_fluxes
        self.get_solutions = self.mesh.get_solution_matrices
        self.save_solutions = self.mesh.save_solution_matrices_to_csv
        self.set_solutions = self.mesh.set_solutions
        self.source_fcn = self.sources.source_fcn

        # Aliasing attributes for convenience
        self.init_fcn_list = param.init_fcn_list
        self.solution_setter_tuple = self.mesh.solution_setter_tuple
        self.time_method = param.time_method

        # Aliasing parameters
        self.cell_vol_m = self.mesh.cell_vol_m
        self.exact_solution_fcn_list = self.sources.exact_solution_fcn_list
        self.solution_shape = (param.Nx + 2, param.Nz + 2)
        self.Nt = param.Nt
        self.x_ct_m, self.z_ct_m = self.mesh.x_ct_m, self.mesh.z_ct_m

        # Compute common parameters
        self.Dt = param.tf / param.Nt
        self.inv_cell_col_m = 1.0 / self.cell_vol_m

        # For showing progress and saving results
        self.csv_period = self.compute_step_period(param.num_csv)
        self.l2_error_period = param.l2_error_period
        self.psg_period = self.compute_step_period(param.num_msg)

    def compute_step_period(self, total_num):
        total_num = min(total_num, self.Nt)

        if total_num <= 0:
            step_period = self.Nt + 1
        else:
            step_period = self.Nt // total_num

        return step_period

    def show_progress(self, step_no, prefix="=> Progress: ", end="\r"):
        if step_no % self.psg_period == 0:
            t = round(step_no * self.Dt, 10)
            print(
                prefix + "%1.2f%%, t = " % (step_no / self.Nt * 100) +
                str(t) + "s  ",
                end=end,
            )

    def save_step_result_to_csv(self, step_no):
        if step_no % self.csv_period != 0:
            return

        t = self.Dt * step_no
        self.save_solutions(t=("%1.4fs" % t), theta_to_T=True)
        print("Saved solutions at t = " + str(round(t, 10)) + "s.    ")

    def compute_l2_error(self, step_no, force=False):
        def l2_norm(sol_m):
            return np.sum(self.cell_vol_m * (sol_m ** 2)) ** 0.5

        def l2_error(numer_sol_m, exact_sol_m):
            abs_l2_error = l2_norm(sol_m=numer_sol_m - exact_sol_m)
            exact_sol_norm = l2_norm(sol_m=exact_sol_m)
            relative_l2_error = abs_l2_error / exact_sol_norm

            return relative_l2_error, abs_l2_error, exact_sol_norm

        exact_fcn_list = self.exact_solution_fcn_list
        if exact_fcn_list is None:
            return

        if force or step_no % self.l2_error_period == 0:
            t = self.Dt * step_no
            xzt_tuple = (self.x_ct_m, self.z_ct_m, t)
            l2_error_norm_tuple = tuple(
                l2_error(
                    numer_sol_m=sol_m,
                    exact_sol_m=exact_fcn(*xzt_tuple)
                )
                for sol_m, exact_fcn in
                zip(self.get_solutions(), exact_fcn_list)
            )

            rel_errors, abs_errors, l2_norms = tuple(zip(*l2_error_norm_tuple))

            msg_prefix = "\nl2 errors at t = " + str(round(t, 10)) + "s:\n"
            indent = " " * 4
            header = indent + " " * 12 + \
                "{:^10} | {:^10} | {:^10} | {:^10}\n".format(
                    "theta", "qv", "qc", "qp"
                )
            msg = msg_prefix + header + \
                indent + "rel error : %1.4e | %1.4e | %1.4e | %1.4e\n" + \
                indent + "abs error : %1.4e | %1.4e | %1.4e | %1.4e\n" + \
                indent + "l2 norm   : %1.4e | %1.4e | %1.4e | %1.4e\n"
            print(msg % (rel_errors + abs_errors + l2_norms))

    def save_exact_solution_to_csv(self, t):
        exact_sol_matrix_list = [
            exact_sol_fcn(self.x_ct_m, self.z_ct_m, t)
            for exact_sol_fcn in self.exact_solution_fcn_list
        ]
        filename_list = ["exact_theta", "exact_qv", "exact_qc", "exact_qp"]
        self.mesh.save_matrices_to_csv(
            matrix_list=exact_sol_matrix_list,
            filename_list=filename_list,
            results_folder_name="solutions",
            filename_suffix=("_%1.4fs" % t),
        )

    def apply_initial_conditions(self):
        init_fcn_list = self.param.init_fcn_list

        if init_fcn_list is None:
            init_fcn_list = self.exact_solution_fcn_list

        if init_fcn_list is None:
            raise Exception("Did not assign correct initial conditions!")

        self.param.save_to_txt()
        self.param.save_to_csv()
        self.mesh.save_coord_matrices_to_csv()
        self.mesh.apply_initial_conditions(init_fcn_list=init_fcn_list)

    def compute_R(self, t):
        S_tuple, V_Ts = self.source_fcn(t=t)
        GG_t_tuple, GG_b_tuple, FF_r_tuple, FF_l_tuple = \
            self.compute_all_fluxes(V_Ts=V_Ts)

        R_tuple = tuple(
            self.inv_cell_col_m * (GG_b - GG_t + FF_l - FF_r) + S
            for GG_t, GG_b, FF_r, FF_l, S in zip(
                GG_t_tuple, GG_b_tuple, FF_r_tuple, FF_l_tuple, S_tuple,
            )
        )

        return R_tuple

    def forward_euler(self):
        Dt = self.Dt
        self.apply_initial_conditions()

        for step_no in range(self.Nt):
            t = Dt * step_no
            self.save_step_result_to_csv(step_no=step_no)
            self.show_progress(step_no=step_no, end="\n")

            solution_tuple = self.get_solutions()
            R_tuple = self.compute_R(t)
            new_solution_tuple = tuple(
                sol + Dt * R
                for sol, R in zip(solution_tuple, R_tuple)
            )

            self.apply_boundary_conditions()
            self.set_solutions(*new_solution_tuple)
            self.compute_l2_error(step_no=step_no + 1)

            # self.save_exact_solution_to_csv(t=t)

        print("\n[Final time]", end=" ")
        self.compute_l2_error(step_no=self.Nt, force=True)
        self.save_step_result_to_csv(step_no=self.Nt)

    def runge_kutta_intermediate_step(
        self, Delta_t, t, solution_copy_tuple, k_m_tuple,
    ):
        R_tuple = self.compute_R(t)

        for k_m, R in zip(k_m_tuple, R_tuple):
            k_m += R

        for sol_copy, R, set_sol in zip(
            solution_copy_tuple, R_tuple, self.solution_setter_tuple,
        ):
            set_sol(sol_copy + Delta_t * R)

        self.apply_boundary_conditions()

    def runge_kutta_final_step(self, Delta_t, t, solution_copy_tuple, k_m_tuple):
        R_tuple = self.compute_R(t)

        for sol_copy, R, k_m, set_sol in zip(
            solution_copy_tuple, R_tuple, k_m_tuple, self.solution_setter_tuple,
        ):
            set_sol(sol_copy + Delta_t * (k_m + R))

        self.apply_boundary_conditions()

    def runge_kutta_4(self):
        solution_shape, Dt = self.solution_shape, self.Dt
        half_Dt = 0.5 * Dt
        self.apply_initial_conditions()
        Delta_t_list = [half_Dt, half_Dt, Dt]

        for step_no in range(self.Nt):
            t = Dt * step_no
            self.show_progress(step_no=step_no, end="\n")
            self.save_step_result_to_csv(step_no=step_no)

            # Step 0
            solution_copy_tuple = tuple(
                np.copy(sol_m) for sol_m in self.get_solutions()
            )
            k_m_tuple = tuple(np.zeros(solution_shape) for _ in range(4))

            # Steps 1-3
            intermediate_t_list = [t + half_Dt, t + half_Dt, t + Dt]
            for Delta_t, interm_t in zip(Delta_t_list, intermediate_t_list):
                self.runge_kutta_intermediate_step(
                    Delta_t=Delta_t,
                    t=interm_t,
                    solution_copy_tuple=solution_copy_tuple,
                    k_m_tuple=k_m_tuple,
                )

            # Step 4
            self.runge_kutta_final_step(
                Delta_t=Dt / 6.0,
                t=t + Dt,
                solution_copy_tuple=solution_copy_tuple,
                k_m_tuple=k_m_tuple,
            )

            self.average_solutions(step_no=step_no + 1)
            self.compute_l2_error(step_no=step_no + 1)

        self.save_step_result_to_csv(step_no=self.Nt)
        print("\n[Final time]", end=" ")
        self.compute_l2_error(step_no=self.Nt, force=True)

    def timestep(self):
        time_method = self.param.time_method
        if time_method == "forward_euler":
            self.forward_euler()
        elif time_method == "rk4":
            self.runge_kutta_4()
        else:
            raise Exception("Unknown time stepping method!")


if __name__ == "__main__":
    pass
