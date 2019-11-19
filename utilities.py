import numpy as np
import os


def pad_zeros_left(M):
    """Pad 1 zero column on the left."""

    zero_col = np.zeros((M.shape[0], 1))

    return np.concatenate((zero_col, M), axis=1)


def pad_zeros_right(M):
    """Pad 1 zero column on the right."""

    zero_col = np.zeros((M.shape[0], 1))

    return np.concatenate((M, zero_col), axis=1)


def pad_zeros_bottom(M):
    """Pad 1 zero row on the bottom."""

    zero_row = np.zeros((1, M.shape[1]))

    return np.concatenate((M, zero_row), axis=0)


def pad_zeros_top(M):
    """Pad 1 zero row on the top."""

    zero_row = np.zeros((1, M.shape[1]))

    return np.concatenate((zero_row, M), axis=0)


def pad_mirror_top(M):
    """Pad 1 row on the top. The padded row mirrors the original first row on
       the top."""

    padded_row = np.reshape(M[0, :], (1, M.shape[1]))

    return np.concatenate((padded_row, M), axis=0)


def get_M_i_minus_one(M, padding="zero"):
    """Return a matrix M_new in which each element M_new[i, j] = M[i - 1, j].
       The top row of M_new is NOT intended to be used. It is padded with zeros
       to keep the original matrix shape."""

    off_M = M[:-1, :]

    if padding == "zero":
        return pad_zeros_top(M=off_M)
    elif padding == "mirror":
        return pad_mirror_top(M=off_M)
    else:
        raise Exception(
            "Wrong padding type: Only 'zero' and 'mirror' are supported!"
        )


def get_M_i_plus_one(M):
    """Return a matrix M_new in which each element M_new[i, j] = M[i + 1, j].
       The bottom row is NOT intended to be used. It is padded with zeros
       to keep the original matrix shape."""

    return pad_zeros_bottom(M=M[1:, :])


def get_M_j_minus_one(M):
    """Return a matrix M_new in which each element M_new[i, j] = M[i, j - 1].
       The left row is NOT intended to be used. It is padded with zeros
       to keep the original matrix shape."""

    return pad_zeros_left(M=M[:, :-1])


def get_M_j_plus_one(M):
    """Return a matrix M_new in which each element M_new[i, j] = M[i, j + 1].
       The right row is NOT intended to be used. It is padded with zeros
       to keep the original matrix shape."""

    return pad_zeros_right(M=M[:, 1:])


def get_function_arg_number(f):
    arg_names = f.__code__.co_varnames
    narg = f.__code__.co_argcount

    if arg_names[0] == "self":
        narg -= 1

    return narg


def is_sequence_type(var):
    return isinstance(var, list) or isinstance(var, tuple)


def print_matrix(M, accuracy=2, format="e", matrix_name="matrix", indent=2):
    nrow, ncol = M.shape
    if format != "e" and format != "f":
        raise Exception("Incorrect format parameter!")

    row_format_str = " " * indent + ("% 2." + str(accuracy) + format + ", ") * ncol
    row_format_str = row_format_str[:-2] + ";"
    print(matrix_name + " = [")
    for row in M:
        print(row_format_str % tuple(row))
    print("]")


def create_folder(folder_name):
    if os.path.exists(folder_name):
        if os.path.isdir(folder_name):
            return
        else:
            raise Exception("A file with the name '%s' exists!" % folder_name)

    os.mkdir(folder_name)
    print("Folder created: %s" % folder_name)


def create_folders(folder_name_list):
    for folder_name in folder_name_list:
        create_folder(folder_name=folder_name)
