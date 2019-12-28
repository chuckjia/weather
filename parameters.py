import utilities as util


class Parameters():
    def __init__(
        self, x0, xf, z0_fcn, z0_der_fcn, zf, Nx, Nz, tf, Nt,
        Mesh, UW, Flux, Boundary, Sources,
        init_fcn_list, rho_o_fcn=None, time_method="rk4",
        num_msg=-1, num_csv=-1, l2_error_period=-1, ave_period_list=-1,
    ):
        self.x0         = x0
        self.xf         = xf
        self.z0_fcn     = z0_fcn
        self.z0_der_fcn = z0_der_fcn
        self.zf         = zf

        self.Nx = Nx
        self.Nz = Nz

        self.tf = tf
        self.Nt = Nt

        self.Mesh     = Mesh
        self.UW       = UW
        self.Flux     = Flux
        self.Boundary = Boundary
        self.Sources  = Sources

        self.init_fcn_list = init_fcn_list
        self.rho_o_fcn     = rho_o_fcn
        self.time_method   = time_method

        self.num_msg         = num_msg
        self.num_csv         = num_csv
        self.l2_error_period = l2_error_period if l2_error_period >= 1 else (Nt + 2)

        if ave_period_list == -1:
            self.ave_period_list = None
        elif isinstance(ave_period_list, int):
            self.ave_period_list = (ave_period_list,) * 4
        elif util.is_sequence_type(ave_period_list) and len(ave_period_list) == 4:
            for ave_period in ave_period_list:
                if not isinstance(ave_period, int):
                    raise Exception("Each ave_period needs to be an integer!")
            self.ave_period_list = tuple(
                ave_period if ave_period >= 1 else (Nt + 2)
                for ave_period in ave_period_list
            )
        else:
            raise Exception("Wrong format for ave_period_list!")

        util.create_folders(("results", "solutions"))


    def get_space_params(self):
        return self.x0, self.xf, self.z0_fcn, self.zf, self.Nx, self.Nz

    def get_time_params(self):
        return self.tf, self.Nt

    def initialize_objects(self):
        mesh = self.Mesh(param=self)
        uw = self.UW(mesh=mesh)
        flux = self.Flux(mesh=mesh, uw=uw)
        boundary = self.Boundary(mesh=mesh)
        sources = self.Sources(mesh=mesh, uw=uw)

        return mesh, uw, flux, boundary, sources

    def __str__(self):
        sep = "===== " * 6
        indent = " " * 4
        info = indent + ("\n" + indent).join([
            "x0 = %1.2f m",
            "xf = %1.2f m",
            "z0 = [function]",
            "zf = %1.2f m",
            "Nx = %d",
            "Nz = %d",
            "tf = %1.4f s",
            "Nt = %d (# of time steps)",
            "",
            "Time method: %s",
            "",
            "# of progress messages = %d",
            "# of saved CSV results = %d",
            "l2 error period        = %d",
            "Averaging period       = %s"
        ]) % (
            self.x0,
            self.xf,
            self.zf,
            self.Nx,
            self.Nz,
            self.tf,
            self.Nt,
            self.time_method,
            self.num_msg,
            self.num_csv,
            self.l2_error_period,
            str(self.ave_period_list),
        )


        return sep + "\nParameters:\n" + info + "\n" + sep

    def save_to_csv(self):
        data = (
            (self.x0, "x0,%1.20e"),
            (self.xf, "xf,%1.20e"),
            (self.zf, "zf,%1.20e"),
            (self.Nx, "Nx,%d"),
            (self.Nz, "Nz,%d"),
            (self.tf, "tf,%1.20e"),
            (self.Nt, "Nt,%d"),
            (self.num_csv, "num_csv,%d"),
        )

        with open("results/param.csv", "w") as f:
            f.write(
                "\n".join([s % d for d, s in data])
            )

    def save_to_txt(self):
        with open("results/param.txt", "w") as f:
            f.write(str(self))


from boundaries import Boundary
from godunov import Flux
from mesh import Mesh
from velocities import UW
from sources import Sources


class DummyParameters(Parameters):
    def __init__(self):
        super().__init__(
            x0=0,
            xf=1,
            z0_fcn=lambda z : 0,
            z0_der_fcn=lambda z : 0,
            zf=1,
            Nx=2,
            Nz=2,
            tf=1,
            Nt=1,
            Mesh=Mesh,
            UW=UW,
            Flux=Flux,
            Boundary=Boundary,
            Sources=Sources,
            init_fcn_list=None,
            rho_o_fcn=None,
            time_method="rk4",
            num_msg=-1,
            num_csv=-1,
            l2_error_period=-1,
            ave_period_list=-1,
        )
