import numpy as np
from .wrapper import ffi, lib


def _float_ptr(x):
    arr = np.ascontiguousarray(x, dtype=np.float64)
    ptr = ffi.cast("double *", arr.ctypes.data)
    return arr, ptr


def _int_ptr(x):
    arr = np.ascontiguousarray(x, dtype=np.int32)
    ptr = ffi.cast("int *", arr.ctypes.data)
    return arr, ptr


class TensorSpline:
    def __init__(self, kind, dim, knots, data, exactdata):
        self.kind = int(kind)
        self.dim = int(dim)

        if len(knots) != self.dim:
            raise ValueError("len(knots) must equal dim")

        knot_sizes = [len(k) for k in knots]
        knots_flat = [v for row in knots for v in row]

        data = np.asarray(data, dtype=np.float64)
        exactdata = np.asarray(exactdata, dtype=np.float64)

        if data.ndim != 2 or data.shape[1] != self.dim + 1:
            raise ValueError("data must have shape (n_rows, dim+1)")

        if exactdata.ndim != 2 or exactdata.shape[1] != self.dim + 1:
            raise ValueError("exactdata must have shape (n_rows, dim+1)")

        self._knot_sizes_np, knot_sizes_ptr = _int_ptr(knot_sizes)
        self._knots_np, knots_ptr = _float_ptr(knots_flat)
        self._data_np, data_ptr = _float_ptr(data.ravel())
        self._exact_np, exact_ptr = _float_ptr(exactdata.ravel())

        self.id = lib.tensor_spline_create(
            self.kind,
            self.dim,
            knots_ptr,
            knot_sizes_ptr,
            data_ptr,
            data.shape[0],
            exact_ptr,
            exactdata.shape[0],
        )

        if self.id < 0:
            raise RuntimeError(f"tensor spline build failed with mode {-self.id}")

    def value(self, point):
        if len(point) != self.dim:
            raise ValueError("point must have length dim")

        self._point_np, point_ptr = _float_ptr(point)
        return lib.tensor_spline_value(self.id, point_ptr, self.dim)

    def value_der(self, point, var):
        if len(point) != self.dim:
            raise ValueError("point must have length dim")
        if not (0 <= var < self.dim):
            raise ValueError("var out of range")

        self._point_np, point_ptr = _float_ptr(point)
        return lib.tensor_spline_value_der(self.id, point_ptr, self.dim, int(var))

    def __del__(self):
        if hasattr(self, "id"):
            lib.tensor_spline_delete(self.id)