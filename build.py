from cffi import FFI
from pathlib import Path

ffibuilder = FFI()

ROOT = Path(__file__).resolve().parent

ffibuilder.cdef("""
    int create_array(double* arr, int n);
    double get_value(int id, int i);
    double sum_array(int id);
    void delete_array(int id);

    int spline_create(double* x, double* y, int n, int type, double tau);
    double spline_value(int id, double t);
    void spline_delete(int id);

    int tensor_spline_create(
        int kind,
        int dim,
        double* knots_flat,
        int* knot_sizes,
        double* data_flat,
        int data_rows,
        double* exact_flat,
        int exact_rows
    );
    double tensor_spline_value(int id, double* point, int dim);
    double tensor_spline_value_der(int id, double* point, int dim, int var);
    void tensor_spline_delete(int id);
""")

ffibuilder.set_source(
    "lstspline._cffi",
    """
    extern "C" {
        int create_array(double* arr, int n);
        double get_value(int id, int i);
        double sum_array(int id);
        void delete_array(int id);

        int spline_create(double* x, double* y, int n, int type, double tau);
        double spline_value(int id, double t);
        void spline_delete(int id);

        int tensor_spline_create(
            int kind,
            int dim,
            double* knots_flat,
            int* knot_sizes,
            double* data_flat,
            int data_rows,
            double* exact_flat,
            int exact_rows
        );
        double tensor_spline_value(int id, double* point, int dim);
        double tensor_spline_value_der(int id, double* point, int dim, int var);
        void tensor_spline_delete(int id);
    }
    """,
    sources=[
        str(ROOT / "c_src" / "wrapper.cpp"),
        str(ROOT / "c_src" / "spline_wrapper.cpp"),
        str(ROOT / "c_src" / "cubicsp.cpp"),
        str(ROOT / "c_src" / "tensor_spline_wrapper.cpp"),
        str(ROOT / "c_src" / "LSNSplineT.cpp"),
        str(ROOT / "c_src" / "NSpline.cpp"),
        str(ROOT / "c_src" / "Lsei.cpp"),
    ],
    source_extension=".cpp",
    extra_compile_args=["/std:c++17"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)