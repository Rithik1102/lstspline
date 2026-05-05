# lstspline

Python package for efficient spline interpolation using C++ (via CFFI), with a strong focus on **analysis, benchmarking, and visualisation of spline behaviour**.

---

## 🚀 Overview

This project implements and evaluates spline interpolation methods, including:

- Cubic spline interpolation
- Monotone-preserving spline interpolation
- Tensor-product splines (multivariate)
- Benchmarking against SciPy implementations

The goal is to study **accuracy, stability, and monotonicity preservation** under different types of data.

---

## 🧠 Motivation

Spline interpolation is widely used in numerical analysis and machine learning. However:

- Cubic splines may introduce **oscillations**
- Monotone splines improve **stability**
- Real-world datasets often contain **noise and discontinuities**

This project explores these behaviours through experiments and comparisons.

---

## 🏗 Architecture

Python → CFFI → C wrapper → C++ classes

- Python: User API and plotting
- CFFI: Interface layer
- C wrapper: Object management using handles
- C++: Core spline algorithms

---

## ✨ Features

### 🔹 Cubic Spline
- Smooth interpolation (C² continuity)
- Efficient C++ implementation
- Suitable for smooth data

### 🔹 Monotone Spline
- Preserves monotonicity
- Prevents overshooting
- Better for noisy / real-world data

### 🔹 Tensor (Multivariate) Spline
- Supports higher-dimensional approximation
- Used for 3D surface modelling

---

## 📊 Project Structure


lstspline/
├── c_src/
├── src/lstspline/
├── tests/
├── plots/
├── generate_all_plots.py
├── benchmark_results.csv
├── build.py


---

## ⚙️ Build

```bash
python build.py
▶️ Run Tests
python tests/test_basic.py
python tests/testcs.py
python tests/test_tensor_spline.py
📊 Generate Plots
python generate_all_plots.py
📈 Benchmarking & Evaluation

The project evaluates spline methods using:

Mean Squared Error (MSE)
Mean Absolute Error (MAE)
Monotonicity Violations
📊 Sample Results
Dataset	Method	MSE	MAE	Monotonicity Violations
Noisy sigmoid	Our monotone	0.000327	0.0143	181
Noisy sigmoid	SciPy PCHIP	0.000320	0.0141	181
Noisy step	Our cubic	0.1674	0.3229	230
Noisy step	SciPy PCHIP	0.1430	0.2979	229
Plateau	Monotone	—	—	0 violations

👉 Monotone splines consistently reduce instability and preserve shape.

📊 Visualisations
1D Spline Comparison

Observation:
Cubic splines introduce oscillations, while monotone splines preserve structure.

Step Data Behaviour

Observation:
Cubic splines overshoot near discontinuities.
Monotone splines provide stable approximation.

3D Tensor Spline

Observation:
Tensor splines approximate smooth surfaces effectively in 2D.

📦 Usage
from lstspline import CubicSpline

x = [0,1,2,3,4,5]
y = [0,0.1,0.2,0.3,0.4,0.5]

s = CubicSpline(x, y)
print(s.value(2.5))
📌 Key Insights
Cubic splines → smooth but unstable near discontinuities
Monotone splines → stable and shape-preserving
SciPy PCHIP → strong benchmark for monotone interpolation
Tensor splines → extend interpolation to higher dimensions
📚 References
R. L. Burden and J. D. Faires, Numerical Analysis
C. de Boor, A Practical Guide to Splines
Fritsch & Carlson, Monotone Piecewise Cubic Interpolation
SciPy Documentation
🎯 Conclusion

This project demonstrates that:

Choosing the correct spline method depends on data characteristics
Monotone splines outperform cubic splines on noisy or discontinuous datasets
Custom C++ implementation performs comparably to SciPy