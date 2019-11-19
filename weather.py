from _model_physical import param
import testing as tttt
from testing import show_time_used
from timesteps import TimeStepping
import time


start_time = time.time()
ts = TimeStepping(param=param)
print(ts.param, "\n")
ts.timestep()
# theta_m = ts.mesh.get_solution_matrices()[0]
# tttt.print_matrix(theta_m[1:-1, 1:-1], accuracy=10, matrix_name="theta")
show_time_used(start_time=start_time)
