# lstspline

Python interface for C++ numerical algorithms using CFFI.

This project demonstrates how to safely integrate Python with C++ by using a C wrapper and handler-based memory management.

---

## Architecture

Python → CFFI → C wrapper → C++ classes

- Python provides user-facing API
- CFFI bridges Python and C
- C wrapper exposes C++ functionality
- C++ manages objects and performs computations

---

## Features

- Pass NumPy arrays to C++ (float64, contiguous memory)
- Safe conversion using `np.ascontiguousarray`
- Handler-based memory management (no raw pointers)
- Multiple methods exposed from C++:
  - `get(i)` → access element
  - `sum()` → compute sum of array

---

## Build

Compile the C++ extension using CFFI:

```bash
python build.py

python tests/test_basic.py