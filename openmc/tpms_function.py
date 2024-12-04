import numpy as np

class FunctionForTPMS:

    def __init__(self):
        pass

    def evaluate(self, x, y, z):
        return 0.

class InterpolationForTPMS(FunctionForTPMS):

    def __init__(self, x_grid = np.array([0.]), y_grid = np.array([0.]), z_grid = np.array([0.]), matrix = np.array([[[0.]]])):
        assert np.ndim(x_grid) == 1
        assert np.ndim(y_grid) == 1
        assert np.ndim(z_grid) == 1
        assert np.ndim(matrix) == 3
        super().__init__()
        self.x_grid, self.y_grid, self.z_grid = x_grid, y_grid, z_grid
        self.matrix = np.array(matrix)

    def evaluate(self, x, y, z):
        def _interp_1d(coord, grid_coords):
            num_points = np.shape(grid_coords)[0]
            grid_indices = np.arange(num_points)
            interp_indices = np.interp(coord, grid_coords, grid_indices).astype(int)
            # Indices for interpolation
            lower_index = interp_indices
            upper_index = np.clip(interp_indices + 1, 0, num_points - 1)
            lower_coord, upper_coord = grid_coords[lower_index], grid_coords[upper_index]
            # Conditions for edge cases
            below_grid = (coord <= grid_coords[0])
            above_grid = (coord > grid_coords[-1])
            within_grid = ~below_grid & ~above_grid
            # Calculate interpolation ratios
            lower_weight = np.full_like(coord, 0.5)
            upper_weight = np.full_like(coord, 0.5)
            lower_weight[within_grid] = (upper_coord[within_grid] - coord[within_grid]) / (upper_coord[within_grid] - lower_coord[within_grid])
            upper_weight[within_grid] = (coord[within_grid] - lower_coord[within_grid]) / (upper_coord[within_grid] - lower_coord[within_grid])
            return lower_index, upper_index, lower_weight, upper_weight
        # Interpolate each dimension
        ix0, ix1, wx0, wx1 = _interp_1d(x, self.x_grid)
        iy0, iy1, wy0, wy1 = _interp_1d(y, self.y_grid)
        iz0, iz1, wz0, wz1 = _interp_1d(z, self.z_grid)
        # Calculate intermediate interpolation weights
        wxy00 = wx0 * wy0
        wxy01 = wx0 * wy1
        wxy10 = wx1 * wy0
        wxy11 = wx1 * wy1
        # Interpolate values
        value_000 = wxy00 * wz0 * self.matrix[ix0, iy0, iz0]
        value_001 = wxy00 * wz1 * self.matrix[ix0, iy0, iz1]
        value_010 = wxy01 * wz0 * self.matrix[ix0, iy1, iz0]
        value_011 = wxy01 * wz1 * self.matrix[ix0, iy1, iz1]
        value_100 = wxy10 * wz0 * self.matrix[ix1, iy0, iz0]
        value_101 = wxy10 * wz1 * self.matrix[ix1, iy0, iz1]
        value_110 = wxy11 * wz0 * self.matrix[ix1, iy1, iz0]
        value_111 = wxy11 * wz1 * self.matrix[ix1, iy1, iz1]
        # Sum contributions
        interpolated_value = value_000 + value_001 + value_010 + value_011 + value_100 + value_101 + value_110 + value_111    
        return interpolated_value