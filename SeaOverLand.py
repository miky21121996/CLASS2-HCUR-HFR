# Paolo Oliveri, 26/04/2016 
# Routine SeaOverLand
# Check if fill_value is correct for the data

import numpy as np
import numpy.ma as ma


def seaoverland(input_matrix, depth=1):  # depth is to select the number of consequential mask points to fill
    # depth loop
    for d in range(depth):
        if np.sum(input_matrix.mask) == 0:  # nothing to fill
            return input_matrix
        else:
            # Create a m x n x 8 3D matrix in which, third dimension fixed, the other dimensions
            #  contains values that are shifted in one of the 8 possible direction compared to the original matrix
            shift_matrix = ma.array(np.empty(shape=(input_matrix.shape[0], input_matrix.shape[1], 8)),
                                    mask=True, fill_value=1.e20, dtype=float)
            # up shift
            shift_matrix[: - 1, :, 0] = input_matrix[1:, :]
            # down shift
            shift_matrix[1:, :, 1] = input_matrix[0: - 1, :]
            # left shift
            shift_matrix[:, : - 1, 2] = input_matrix[:, 1:]
            # right shift
            shift_matrix[:, 1:, 3] = input_matrix[:, : - 1]
            # up-left shift
            shift_matrix[: - 1, : - 1, 4] = input_matrix[1:, 1:]
            # up-right shift
            shift_matrix[: - 1, 1:, 5] = input_matrix[1:, : - 1]
            # down-left shift
            shift_matrix[1:, : - 1, 6] = input_matrix[: - 1, 1:]
            # down-right shift
            shift_matrix[1:, 1:, 7] = input_matrix[: - 1, : - 1]
            # Mediate the shift matrix among the third dimension
            mean_matrix = ma.mean(shift_matrix, 2)
            # Replace input missing values with new ones belonging to the mean matrix
            input_matrix = ma.array(np.where(mean_matrix.mask + input_matrix.mask, mean_matrix, input_matrix),
                                    mask=mean_matrix.mask, fill_value=1.e20, dtype=float)
            input_matrix = ma.masked_where(mean_matrix.mask, input_matrix)
    return input_matrix
