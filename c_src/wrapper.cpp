#include "myclass.cpp"

extern "C" {

    void* create_array(double* arr, int n) {
        return new MyArray(arr, n);
    }

    double get_value(void* obj, int i){
        return ((MyArray*)obj)->get(i);
    }

    void delete_array(void* obj) {
        delete (MyArray*)obj;
    }
}