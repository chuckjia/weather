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

    def apply_left_dirichlet_condition(self, boundary_val_array):
        for sol_m in self.get_solutions():
            sol_m[0, :] = 2 * boundary_val_array - sol_m[1, :]

    def apply_all_boundary_zero_dirichlet_condition(self):
        for sol_m in self.get_solutions():
            sol_m[0, :] = -sol_m[1, :]
            sol_m[-1, :] = -sol_m[-2, :]
            sol_m[:, 0] = -sol_m[:, 1]
            sol_m[:, -1] = -sol_m[:, -2]

    def apply_right_neumann_condition(self):
        for sol_m in self.get_solutions():
            sol_m[-1, :] = sol_m[-2, :]


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
        self.twice_left_boundary_val_array_tuple = tuple(
            sol_m[0, :] + sol_m[1, :] for sol_m in self.get_solutions()
        )

    def apply_boundary_conditions(self):
        # Left Dirichlet boundary condition
        for sol_m, twice_boundary_val in zip(
            self.get_solutions(), self.twice_left_boundary_val_array_tuple,
        ):
            sol_m[0, :] = twice_boundary_val - sol_m[1, :]

        # Right Neumann boundary condition
        self.apply_right_neumann_condition()


if __name__ == "__main__":
    pass
