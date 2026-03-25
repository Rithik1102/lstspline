#include "cubicsp.h"
#include <map>


static std::map<int, cubic_spline*> spline_objects;
static int spline_current_id = 0;

extern "C" {

int spline_create(double* x, double* y, int n, int type, double tau) {
    cubic_spline* s = new cubic_spline();
    s->BuildSpline(x, y, n, type, tau);

    spline_current_id++;
    spline_objects[spline_current_id] = s;

    return spline_current_id;
}

double spline_value(int id, double t) {
    return spline_objects[id]->SplineValue(t);
}

void spline_delete(int id) {
    delete spline_objects[id];
    spline_objects.erase(id);
}

}