class Boundary():
    def __init__(self, mesh):
        self.param = mesh.param
        self.mesh = mesh
        self.get_solutions = mesh.get_solution_matrices

    def apply_boundary_conditions(self):
        return

    def apply_left_zero_dirichlet_condition(self):
        for sol_m in self.get_solutions():
            sol_m[0, :] = - sol_m[1, :]

    def apply_left_dirichlet_condition(self, boundary_val_array_tuple):
        sol_list = self.get_solutions()
        for sol_m, boundary_vals in zip(sol_list, boundary_val_array_tuple):
            sol_m[0, :] = 2 * boundary_vals - sol_m[1, :]

    def apply_left_dirichlet_condition_alt(self, boundary_val_array_tuple):
        sol_list = self.get_solutions()
        for sol_m, boundary_vals in zip(sol_list, boundary_val_array_tuple):
            sol_m[0, :] = boundary_vals

    def apply_all_boundary_zero_dirichlet_condition(self):
        for sol_m in self.get_solutions():
            sol_m[0, :] = -sol_m[1, :]
            sol_m[-1, :] = -sol_m[-2, :]
            sol_m[:, 0] = -sol_m[:, 1]
            sol_m[:, -1] = -sol_m[:, -2]

    def apply_right_neumann_condition(self):
        for sol_m in self.get_solutions():
            sol_m[-1, :] = sol_m[-2, :]

    def apply_bott_neumann_condition(self):
        for sol_m in self.get_solutions():
            sol_m[:, 0] = sol_m[:, 1]


class ZeroDirichletBoundary(Boundary):
    def __init__(self, mesh):
        super().__init__(mesh=mesh)

    def apply_boundary_conditions(self):
        self.apply_all_boundary_zero_dirichlet_condition()


class PhysicalBoundary(Boundary):
    def __init__(self, mesh):
        super().__init__(mesh=mesh)
        mesh.apply_initial_conditions(init_fcn_list=self.param.init_fcn_list)
        self.set_boundary_values()

    def set_boundary_values(self):
        self.left_boundary_val_array_tuple = tuple(
            (sol_m[0, :] + sol_m[1, :]) / 2 for sol_m in self.get_solutions()
        )

    def smooth_bott_boundary(self):
        theta_m = self.get_solutions()[0]
        theta_m[:, 1] = theta_m[:, 2]

    def apply_boundary_conditions(self):
        # Left Dirichlet boundary condition
        self.apply_left_dirichlet_condition_alt(
            boundary_val_array_tuple=self.left_boundary_val_array_tuple,
        )

        # Right Neumann boundary condition
        self.apply_right_neumann_condition()

        self.apply_bott_neumann_condition()
        self.smooth_bott_boundary()


if __name__ == "__main__":
    pass
