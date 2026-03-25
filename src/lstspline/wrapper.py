import numpy as np
from pathlib import Path
import importlib.util
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BUILT_DIR = PROJECT_ROOT / "lstspline"

matches = list(BUILT_DIR.glob("_cffi*.pyd")) + list(BUILT_DIR.glob("_cffi*.so"))
if not matches:
    raise ModuleNotFoundError(
        "Could not find built CFFI module in the default location. Run 'python build.py' first."
    )

_cffi_path = matches[0]

spec = importlib.util.spec_from_file_location("lstspline_built._cffi", _cffi_path)
if spec is None or spec.loader is None:
    raise ImportError(f"Could not load CFFI module from {_cffi_path}")

_cffi = importlib.util.module_from_spec(spec)
sys.modules["lstspline_built._cffi"] = _cffi
spec.loader.exec_module(_cffi)

ffi = _cffi.ffi
lib = _cffi.lib


def convert_py_float_to_cffi(x):
    if x is not None:
        if isinstance(x, np.ndarray) and x.flags.c_contiguous and x.dtype == np.float64:
            px = x
        else:
            px = np.ascontiguousarray(x, dtype="float64")
        pxcffi = ffi.cast("double *", px.ctypes.data)
    else:
        px = np.array([0.0], dtype="float64")
        pxcffi = ffi.cast("double *", 0)
    return px, pxcffi


class MyArray:
    def __init__(self, arr):
        self.arr, c_ptr = convert_py_float_to_cffi(arr)
        self.n = len(self.arr)
        self.id = lib.create_array(c_ptr, self.n)

    def get(self, i):
        return lib.get_value(self.id, i)

    def sum(self):
        return lib.sum_array(self.id)

    def __del__(self):
        if hasattr(self, "id"):
            lib.delete_array(self.id)


class CubicSpline:
    def __init__(self, x, y, type=1, tau=0.0):
        self.x, x_ptr = convert_py_float_to_cffi(x)
        self.y, y_ptr = convert_py_float_to_cffi(y)

        self.n = len(self.x)

        self.id = lib.spline_create(x_ptr, y_ptr, self.n, type, tau)

    def value(self, t):
        return lib.spline_value(self.id, t)

    def __del__(self):
        if hasattr(self, "id"):
            lib.spline_delete(self.id)