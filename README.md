# lstspline

Toy Python package demonstrating C++ integration using CFFI.

## Architecture

Python → CFFI → C wrapper → C++

## Features

- Pass NumPy arrays to C++
- Store data in C++ class
- Access via Python methods

## Example

```python
from lstspline import MyArray

arr = MyArray([1,2,3])
print(arr.get(1))  # Output: 2