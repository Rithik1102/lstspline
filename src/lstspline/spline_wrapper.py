from .wrapper import convert_py_float_to_cffi, lib


class CubicSpline:
    def __init__(self, x, y, type=1, tau=0.0):
        self.x, x_ptr = convert_py_float_to_cffi(x)
        self.y, y_ptr = convert_py_float_to_cffi(y)

        if len(self.x) != len(self.y):
            raise ValueError("x and y must have the same length")

        self.n = len(self.x)
        self.id = lib.spline_create(x_ptr, y_ptr, self.n, type, tau)

    def value(self, t):
        return lib.spline_value(self.id, float(t))

    def __del__(self):
        if hasattr(self, "id"):
            lib.spline_delete(self.id)