# lstspline

Python package for numerical spline algorithms using C++ with CFFI.

This project demonstrates how to integrate Python with C++ for efficient
numerical computation, with a focus on cubic spline interpolation.

------------------------------------------------------------------------

## Architecture

Python → CFFI → C wrapper → C++ classes

-   Python provides user-facing API\
-   CFFI bridges Python and C\
-   C wrapper manages C++ objects via integer handlers\
-   C++ implements spline algorithms

------------------------------------------------------------------------

## Features

### Toy Example

-   `MyArray` class for demonstrating:
    -   Passing NumPy arrays to C++
    -   Accessing stored data
    -   Performing computations (`get`, `sum`)

### Cubic Spline (Core Implementation)

-   Dynamic memory using `std::vector<double>`
-   Removed fixed-size array limitations (`MAXKNOTS`)
-   Build spline from input data
-   Evaluate spline at arbitrary points

------------------------------------------------------------------------

## Project Structure

    lstspline/
    ├── c_src/
    │   ├── cubicsp.cpp
    │   ├── cubicsp.h
    │   ├── wrapper.cpp
    │   ├── spline_wrapper.cpp
    ├── src/lstspline/
    │   ├── __init__.py
    │   ├── wrapper.py
    │   ├── spline_wrapper.py
    ├── tests/
    │   └── test_basic.py
    ├── build.py

------------------------------------------------------------------------

## Build

    python build.py

------------------------------------------------------------------------
## Test
   
   python tests/test_basic.py  (For toy example)
   python tests/testcs.py      (For cubic spline implementation)
   python tests/test_tensor_spline.py (For Multivariate (Tensor) Spline implementation)
------------------------------------------------------------------------

## Usage

``` python
from lstspline import CubicSpline

x = [0,1,2,3,4,5]
y = [0,0.1,0.2,0.3,0.4,0.5]

s = CubicSpline(x, y)
print(s.value(2.5))   # ~0.25
```

------------------------------------------------------------------------

## Notes

-   Input arrays are converted to contiguous `float64`
-   C++ objects are managed using integer handlers
-   Supports interpolation between points

------------------------------------------------------------------------

## Multivariate (Tensor) Spline

The package also supports multivariate tensor-product splines.

### Example

```python
from lstspline import TensorSpline

knots = [
    [0.0, 0.5, 1.0],
    [0.0, 0.5, 1.0],
]

data = [
    [0.2, 0.3, 0.06],
    [0.7, 0.4, 0.28],
    [0.5, 0.5, 0.25],
    [0.9, 0.8, 0.72],
]

exactdata = [
    [0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0],
]

s = TensorSpline(kind=0, dim=2, knots=knots, data=data, exactdata=exactdata)

print(s.value([0.5, 0.5]))
print(s.value_der([0.5, 0.5], 0))
