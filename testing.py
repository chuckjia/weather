import numpy as np
import time


def show_time_used(start_time):
    print("Time used: %fs." % (time.time() - start_time))


def show_progress(progress):
    print("Progress: %1.2f%%" % progress, end="\r")


def validate_grid_matrices(
    x_bl_m, x_br_m, x_tl_m, x_tr_m, z_bl_m, z_br_m, z_tl_m, z_tr_m,
):
    def check_matrices_are_equal(matrix_list, msg_prefix=""):
        msg_prefix = msg_prefix.strip()
        for i in range(len(matrix_list) - 1):
            for j in range(i + 1, len(matrix_list)):
                if not np.array_equal(matrix_list[i], matrix_list[j]):
                    msg = ": Matrices %d and %d are not the same." % (i, j)
                    print(msg_prefix + msg)

    x_bl_m = x_bl_m[1:, 1:]
    z_bl_m = z_bl_m[1:, 1:]

    x_br_m = x_br_m[:-1, 1:];
    z_br_m = z_br_m[:-1, 1:];

    x_tl_m = x_tl_m[1:, :-1]
    z_tl_m = z_tl_m[1:, :-1]

    x_tr_m = x_tr_m[:-1, :-1]
    z_tr_m = z_tr_m[:-1, :-1]

    x_matrix_list = [x_bl_m, x_br_m, x_tl_m, x_tr_m]
    z_matrix_list = [z_bl_m, z_br_m, z_tl_m, z_tr_m]
    check_matrices_are_equal(x_matrix_list, msg_prefix="For x")
    check_matrices_are_equal(z_matrix_list, msg_prefix="For z")

    print("Checking completed.")


def remove_ghost_cells(M):
    return M[1:-1, 1:-1]
