# lstspline

A Python package for high-performance spline interpolation using C++ and CFFI, developed as part of a Minor Thesis research project investigating spline interpolation, monotonicity preservation, tension control, and multidimensional approximation.

---

# Overview

`lstspline` provides a collection of spline interpolation methods implemented in C++ and exposed to Python through CFFI. The project focuses on both implementation and empirical evaluation of spline behaviour across synthetic and real-world datasets.

Implemented methods include:

- Cubic spline interpolation
- Monotone spline interpolation
- Taut spline interpolation
- Tensor-product splines for multidimensional data
- Benchmarking against SciPy spline implementations

The project evaluates interpolation quality, shape preservation, monotonicity, numerical stability, and computational performance.

---

# Research Motivation

Spline interpolation is widely used in numerical analysis, scientific computing, machine learning, computer graphics, and data modelling.

Traditional cubic splines provide smooth approximations but may introduce:

- Overshooting
- Undershooting
- Oscillatory behaviour
- Violations of monotonicity

These problems become particularly noticeable when interpolating:

- Economic time series
- Population growth data
- Noisy measurements
- Piecewise monotone functions

Shape-preserving interpolation methods such as monotone splines, taut splines, and PCHIP were developed to address these limitations.

This project investigates how different spline methods behave under varying data characteristics and evaluates their suitability for practical applications.

---

# Architecture

```text
Python API
    ↓
CFFI Interface
    ↓
C Wrapper Layer
    ↓
C++ Spline Implementation
```

---

# Features

## Cubic Spline

- Piecewise cubic interpolation
- Continuous first and second derivatives
- Smooth approximation
- Suitable for smooth datasets

## Monotone Spline

- Preserves monotonicity
- Eliminates overshooting
- Shape-preserving interpolation
- Appropriate for growth trends and cumulative data

## Taut Spline

- Uses tension parameters to control curvature
- Reduces oscillatory behaviour
- Provides balance between flexibility and stability

## Tensor-Product Spline

- Multidimensional interpolation
- Surface approximation
- Scientific and engineering applications
- Supports higher-dimensional spline modelling

---

# Project Structure

```text
lstspline/
│
├── c_src/
├── src/lstspline/
├── tests/
├── plots/
├── generate_all_plots.py
├── benchmark_results.csv
├── build.py
└── README.md
```

---

# Installation

```bash
git clone https://github.com/Rithik1102/lstspline.git
cd lstspline
python build.py
```

---

# Running Tests

```bash
python tests/test_basic.py
python tests/testcs.py
python tests/test_tensor_spline.py
```

---

# Usage Example

```python
from lstspline import CubicSpline

x = [0, 1, 2, 3, 4, 5]
y = [0, 0.1, 0.2, 0.3, 0.4, 0.5]

spline = CubicSpline(x, y)

print(spline.value(2.5))
```

---

# Experimental Evaluation and Thesis Results

The project evaluates spline methods using both synthetic and real-world datasets.

## Synthetic Datasets

The following datasets were designed to investigate spline behaviour under controlled conditions:

- Slow monotone growth
- Plateau then sharp increase
- Irregular monotone increase
- Noisy sigmoid data
- Noisy step data
- Noisy piecewise growth
- Smooth monotone curves
- Tensor-product surface datasets

These datasets were selected to evaluate:

- Monotonicity preservation
- Overshooting behaviour
- Shape preservation
- Tension control
- Numerical stability

## Real-World Datasets

### Australian GDP Dataset

Used to analyse long-term economic growth trends and interpolation accuracy on smooth monotonic data.

### Melbourne Housing Dataset

Used to evaluate interpolation behaviour on nonlinear real-estate pricing data.

### Australian Inflation Dataset

Used to investigate interpolation on economic time-series data with local fluctuations.

### Australian Population Dataset

Used to assess shape-preserving interpolation for naturally increasing demographic trends.

## Tensor Spline Experiments

Multidimensional synthetic surfaces were generated to evaluate tensor-product spline interpolation.

The experiments examined:

- Surface reconstruction accuracy
- Smoothness
- Multivariate interpolation behaviour
- Numerical performance

---

# Evaluation Metrics

The spline methods were compared using:

- Mean Squared Error (MSE)
- Mean Absolute Error (MAE)
- Root Mean Squared Error (RMSE)
- Monotonicity Violations
- Runtime Performance

---

# Key Findings

## Cubic Splines

- Produce very smooth interpolations
- Can overshoot near sharp transitions
- May violate monotonicity

## Monotone Splines

- Preserve data ordering
- Eliminate overshooting
- Produce stable approximations

## Taut Splines

- Reduce excessive curvature
- Improve local control
- Provide compromise between smoothness and stability

## SciPy PCHIP

- Strong shape-preserving baseline
- Comparable performance to monotone spline methods

## Tensor Splines

- Successfully approximate multidimensional surfaces
- Extend spline interpolation beyond one-dimensional problems

## Real-World Dataset Results

Across the GDP, inflation, population, and housing datasets:

- Monotone splines consistently preserved underlying trends.
- Cubic splines occasionally introduced oscillations between sparse observations.
- PCHIP and monotone splines produced reliable shape-preserving approximations.
- Taut splines provided a useful compromise between smoothness and stability.
- No single method was optimal for every dataset.

---

# Thesis Outcomes

This research demonstrates that:

- No single spline method is optimal for all datasets.
- Cubic splines provide smooth interpolation but may introduce undesirable oscillations.
- Monotone splines consistently preserve underlying trends and eliminate overshoot.
- Taut splines improve shape control through tension adjustment.
- Tensor-product splines effectively extend interpolation techniques to multidimensional domains.
- The custom C++ implementation performs competitively against established SciPy implementations.

---

# Future Work

Potential future developments include:

- Adaptive tension parameter selection
- Higher-dimensional tensor spline extensions
- Automatic monotonicity detection
- GPU acceleration
- Integration with machine learning workflows
- Additional shape-preserving spline variants

---

# References

1. de Boor, C. *A Practical Guide to Splines*. Springer-Verlag, 1978.
2. Schumaker, L. L. *Spline Functions: Basic Theory*. Cambridge University Press, 2007.
3. Schumaker, L. L. "On Shape-Preserving Quadratic Spline Interpolation." SIAM Journal on Numerical Analysis, 1983.
4. Sapidis, N. S., and Kaklis, P. D. "An Algorithm for Constructing Convexity and Monotonicity Preserving Splines in Tension." Computer Aided Geometric Design, 1988.
5. Sapidis, N. S., Kaklis, P. D., and Loukakis, T. A. "A Method for Computing the Tension Parameters in Convexity-Preserving Spline-in-Tension Interpolation." Numerische Mathematik, 1988.
6. Fritsch, F. N., and Carlson, R. E. "Monotone Piecewise Cubic Interpolation." SIAM Journal on Numerical Analysis, 1980.

---

# Author

**Rithik Vishal Nair**

SIT792 Minor Thesis  
School of Information Technology  
Deakin University

**Supervisor:** Professor Gleb Beliakov