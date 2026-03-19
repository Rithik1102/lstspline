#include <vector>

class MyArray{
    public:
        std::vector<double> data;

        MyArray(double* arr, int n) {
            data = std::vector<double>(arr, arr + n);
        }

        double get(int i) {
            return data[i];
        }

        double sum() {
            double s = 0;
            for (double v : data) s += v;
            return s;
        }


};