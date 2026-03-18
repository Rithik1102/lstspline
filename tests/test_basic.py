import sys
import os

# Add src to path
sys.path.append(os.path.abspath("src"))

from lstspline import MyArray

arr = MyArray([1, 2, 3, 4])

print("Array values:")
for i in range(4):
    print(f"Index {i}:", arr.get(i))