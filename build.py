from cffi import FFI
from pathlib import Path

ffibuilder = FFI()

ROOT = Path(__file__).resolve().parent
WRAPPER_CPP = ROOT / "c_src" / "wrapper.cpp"

ffibuilder.cdef("""
    int create_array(double* arr, int n);
    double get_value(int id, int i);
    double sum_array(int id);
    void delete_array(int id);
""")

ffibuilder.set_source(
    "lstspline._cffi",
    """
    extern "C" {
        int create_array(double* arr, int n);
        double get_value(int id, int i);
        double sum_array(int id);
        void delete_array(int id);
    }
    """,
    sources=[str(WRAPPER_CPP)],
    source_extension=".cpp",
    extra_compile_args=["/std:c++17"],
)

if __name__ == "__main__":
    ffibuilder.compile(
        verbose=True,
        tmpdir="build",
        target="src/lstspline/_cffi.*"
    )