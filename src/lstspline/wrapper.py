from cffi import FFI
import numpy as np
import os

ffi = FFI()

ffi.cdef("""
    void* create_array(double* arr, int n);
    double get_value(void* obj, int i);
    void delete_array(void* obj);
""")

# Detect OS and load correct library
if os.name == "nt":
    lib_name = "myarray.dll"
else:
    lib_name = "libmyarray.so"

lib_path = os.path.abspath(lib_name)
lib = ffi.dlopen(lib_path)


class MyArray:
    def __init__(self, arr):
        arr = np.array(arr, dtype='float64')
        self.n = len(arr)

        self.c_arr = ffi.cast("double*", arr.ctypes.data)
        self.obj = lib.create_array(self.c_arr, self.n)

    def get(self, i):
        return lib.get_value(self.obj, i)

    def __del__(self):
        if hasattr(self, "obj"):
            lib.delete_array(self.obj)