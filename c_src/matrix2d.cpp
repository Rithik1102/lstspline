
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <span>
typedef vector<double>	dVec;
typedef vector<dVec*> VecArray;



template <typename T>
class Array2D {
private:
    std::vector<T> Data;
    size_t rows_;
    size_t cols_;

    // Proxy class for operator[]
    class RowProxy {
    private:
        T* row_data;
        size_t cols;
    public:
        RowProxy(T* data, size_t cols) : row_data(data), cols(cols) {}
        
        T& operator[](size_t col) {
            if (col >= cols) throw std::out_of_range("Column index out of bounds");
            return row_data[col];
        }
        
        const T& operator[](size_t col) const {
            if (col >= cols) throw std::out_of_range("Column index out of bounds");
            return row_data[col];
        }
    };

public:
    // Constructors
    Array2D() {};
    Array2D(size_t rows, size_t cols, const T& init_val = T())
        : rows_(rows), cols_(cols), Data(rows * cols, init_val) {}

    Array2D(std::initializer_list<std::initializer_list<T>> init) {
        rows_ = init.size();
        if (rows_ == 0) {
            cols_ = 0;
            return;
        }
        cols_ = init.begin()->size();
        Data.reserve(rows_ * cols_);
        
        for (const auto& row : init) {
            if (row.size() != cols_) {
                throw std::invalid_argument("All rows must have the same number of columns");
            }
            Data.insert(Data.end(), row.begin(), row.end());
        }
    }

    // Access methods with operator[]
    RowProxy operator[](size_t row) {
        if (row >= rows_) throw std::out_of_range("Row index out of bounds");
        return RowProxy(&Data[row * cols_], cols_);
    }

    const RowProxy operator[](size_t row) const {
        if (row >= rows_) throw std::out_of_range("Row index out of bounds");
        return RowProxy(const_cast<T*>(&Data[row * cols_]), cols_);
    }

    // Access methods with operator()
    T& operator()(size_t row, size_t col) {
        row--; col--;
        check_bounds(row, col);
        return Data[row * cols_ + col];
    }

    const T& operator()(size_t row, size_t col) const {
        row--; col--;
        check_bounds(row, col);
        return Data[row * cols_ + col];
    }

    // Size information
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }
    size_t size() const { return Data.size(); }

    size_t dim(int i) const { return  ((i == 1) ? rows_ : cols_); }

    void newsize(int i, int j) {
        Data.clear();
        Data.resize(i*j);
        rows_ = i;   cols_ = j;
    }

    T* data() { return Data.data(); }
    const T* operator()(size_t row)  {
        return &Data[row * cols_ ];
    }

    // Iterators
    typename std::vector<T>::iterator begin() { return Data.begin(); }
    typename std::vector<T>::iterator end() { return Data.end(); }
    typename std::vector<T>::const_iterator begin() const { return Data.begin(); }
    typename std::vector<T>::const_iterator end() const { return Data.end(); }

    // this one with a copy, how to clean up? will be destroyed automatically by vector?
    void MakeVecArray(VecArray& out)
    {
        dVec* r;
        for (int i = 0; i < rows_; i++) {
            r = new dVec;
            r->assign(&Data[i * cols_], &Data[(i ) * cols_ + cols_ -1] +1 );
            out.push_back(r);
        }
    }

    //GB for now working on that part
    void MakeVecArray1(VecArray& out)
    {
        for (int i = 0; i < rows_; i++) {
            // cout << i << endl;
            span<double> *r=new span<double>(&Data[i * cols_], &Data[(i)*cols_ + cols_ - 1] + 1);
            out.push_back( (dVec*)( r));
        }
    }
// we can use span without copying in C++20


private:
    void check_bounds(size_t row, size_t col) const {
        if (row >= rows_ || col >= cols_) {
            throw std::out_of_range("Array2D index out of bounds");
        }
    }
};

