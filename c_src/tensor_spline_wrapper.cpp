#include "LSNSplineT.h"
#include <map>
#include <vector>

static std::map<int, LSNSplineT*> tensor_spline_objects;
static int tensor_spline_current_id = 0;

static VecArray make_knots_array_from_flat(double* flat, int* knot_sizes, int dim)
{
    VecArray out;
    out.reserve(dim);

    int offset = 0;
    for (int i = 0; i < dim; ++i) {
        dVec* row = new dVec;
        row->resize(knot_sizes[i]);
        for (int j = 0; j < knot_sizes[i]; ++j) {
            (*row)[j] = flat[offset + j];
        }
        offset += knot_sizes[i];
        out.push_back(row);
    }
    return out;
}

static VecArray make_data_array_from_flat(double* flat, int rows, int cols)
{
    VecArray out;
    out.reserve(rows);

    for (int i = 0; i < rows; ++i) {
        dVec* row = new dVec;
        row->resize(cols);
        for (int j = 0; j < cols; ++j) {
            (*row)[j] = flat[i * cols + j];
        }
        out.push_back(row);
    }
    return out;
}

extern "C" {

int tensor_spline_create(
    int kind,
    int dim,
    double* knots_flat,
    int* knot_sizes,
    double* data_flat,
    int data_rows,
    double* exact_flat,
    int exact_rows
) {
    LSNSplineT* s = new LSNSplineT();

    VecArray knots = make_knots_array_from_flat(knots_flat, knot_sizes, dim);
    VecArray data = make_data_array_from_flat(data_flat, data_rows, dim + 1);
    VecArray exactdata = make_data_array_from_flat(exact_flat, exact_rows, dim + 1);

    int mode = s->Build(kind, dim, knots, data, exactdata);

    deleteVecArray(knots);
    deleteVecArray(data);
    deleteVecArray(exactdata);

    if (mode != 0) {
        delete s;
        return -mode;
    }

    tensor_spline_current_id++;
    tensor_spline_objects[tensor_spline_current_id] = s;
    return tensor_spline_current_id;
}

double tensor_spline_value(int id, double* point, int dim)
{
    dVec t(dim);
    for (int i = 0; i < dim; ++i) {
        t[i] = point[i];
    }
    return tensor_spline_objects[id]->Value(t);
}

double tensor_spline_value_der(int id, double* point, int dim, int var)
{
    dVec t(dim);
    for (int i = 0; i < dim; ++i) {
        t[i] = point[i];
    }
    return tensor_spline_objects[id]->ValueDer(t, var);
}

void tensor_spline_delete(int id)
{
    auto it = tensor_spline_objects.find(id);
    if (it != tensor_spline_objects.end()) {
        delete it->second;
        tensor_spline_objects.erase(it);
    }
}

}