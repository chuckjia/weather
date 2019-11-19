import utilities as util
import numpy as np
import testing
import time


class Flux():
    def __init__(self, mesh, uw):
        start_time = time.time()

        self.param = mesh.param
        self.mesh = mesh
        self.uw = uw

        self.get_solutions = mesh.get_solution_matrices

        # For top fluxes
        # top_flux_velocity is everything in (3.26) except u_check, i.e.
        # |\Gamma_{top}| \vec{n}_{top} \cdot (u_{top}, w_{top})
        self.top_flux_velocity = mesh.top_side_len_m * (
            mesh.n_top_x_m * uw.u_t_m + mesh.n_top_z_m * uw.w_t_m
        )
        self.top_forward_cells = self.top_flux_velocity >= 0
        self.top_back_cells = self.top_flux_velocity < 0

        # Matrices for right fluxes
        self.right_flux_velocity = mesh.right_side_len_m * uw.u_r_m
        self.right_forward_cells = self.right_flux_velocity >= 0
        self.right_back_cells = self.right_flux_velocity < 0

        print(
            "Initialized fluxes. "
            "Time used: %1.2fs." % (time.time() - start_time)
        )

    def top_flux_source(self, sol_m, forward_cells, back_cells):
        """The function u_check as in (3.27)."""

        off_sol_m = util.pad_zeros_right(M=sol_m[:, 1:])
        return forward_cells * sol_m + back_cells * off_sol_m

    def compute_top_fluxes(self):
        return tuple(
            self.top_flux_velocity * self.top_flux_source(
                sol_m=sol_m,
                forward_cells=self.top_forward_cells,
                back_cells=self.top_back_cells,
            )
            for sol_m in self.get_solutions()
        )

    def get_bott_fluxes(self, GG_list):
        return tuple(util.pad_zeros_left(M=GG[:, :-1]) for GG in GG_list)

    def right_flux_source(self, sol_m):
        """The function u_check as in (3.30)."""

        off_sol_m = util.pad_zeros_bottom(M=sol_m[1:, :])
        return self.right_forward_cells * sol_m + \
            self.right_back_cells * off_sol_m

    def compute_right_fluxes(self):
        return tuple(
            self.right_flux_velocity * self.right_flux_source(sol_m=sol_m)
            for sol_m in self.get_solutions()
        )

    def get_left_fluxes(self, FF_list):
        return tuple(util.pad_zeros_top(M=FF[:-1, :]) for FF in FF_list)


class PhysicalFlux(Flux):
    def __init__(self, mesh, uw):
        super().__init__(mesh, uw)

        self.GG_qp_factor = mesh.top_side_len_m * mesh.n_top_z_m

    def compute_top_fluxes(self, V_Ts):
        """Override the classical Godunov fluxes, as the velocities for qp is
           different."""

        solutions_list = self.get_solutions()

        top_flux_velocity_qp = self.top_flux_velocity - self.GG_qp_factor * V_Ts
        top_forward_cells_qp = top_flux_velocity_qp >= 0
        top_back_cells_qp = top_flux_velocity_qp < 0

        flux_velocity_list = [self.top_flux_velocity] * 3 + [top_flux_velocity_qp]
        forward_cells_list = [self.top_forward_cells] * 3 + [top_forward_cells_qp]
        back_cells_list = [self.top_back_cells] * 3 + [top_back_cells_qp]

        return tuple(
            flux_v * self.top_flux_source(
                sol_m=sol_m, forward_cells=fc, back_cells=bc,
            )
            for sol_m, flux_v, fc, bc in zip(
                solutions_list,
                flux_velocity_list,
                forward_cells_list,
                back_cells_list,
            )
        )

    def compute_all_fluxes(self, V_Ts):
        GG_t_tuple = self.compute_top_fluxes(V_Ts=V_Ts)
        GG_b_tuple = self.get_bott_fluxes(GG_list=GG_t_tuple)
        FF_r_tuple = self.compute_right_fluxes()
        FF_l_tuple = self.get_left_fluxes(FF_list=FF_r_tuple)

        return GG_t_tuple, GG_b_tuple, FF_r_tuple, FF_l_tuple


if __name__ == "__main__":
    pass
