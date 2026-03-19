#include "myclass.cpp"
#include <map>


std::map<int, MyArray*> objects;
int current_id = 0;

extern "C" {

    int create_array(double* arr, int n) {
        MyArray* obj = new MyArray(arr, n);
        current_id++;
        objects[current_id] = obj;
        return current_id;
    }

    double get_value(int id, int i){
        return objects[id]->get(i);
    }

    double sum_array(int id){
        return objects[id]->sum();
    }

    void delete_array(int id) {
        delete objects[id];
        objects.erase(id);
    }
}