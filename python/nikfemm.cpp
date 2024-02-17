#include <iostream>
#include <string>
#include <functional>
#include <iomanip>

#include <nikfemm.hpp>

#include <Python.h>

typedef struct {
    PyObject_HEAD
    nikfemm::System<double>* system;
} nikfemm_System;

static void nikfemm_System_dealloc(nikfemm_System* self) {
    delete self->system;
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* nikfemm_System_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    nikfemm_System* self;

    self = (nikfemm_System*)type->tp_alloc(type, 0);

    return (PyObject*)self;
}

static PyTypeObject nikfemm_SystemType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "nikfemm.System",
    .tp_basicsize = sizeof(nikfemm_System),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)nikfemm_System_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_doc = PyDoc_STR("System of linear equations for current density"),
    .tp_new = nikfemm_System_new,
};

typedef struct {
    PyObject_HEAD
    nikfemm::MultiLayerCurrentDensitySimulation* simulation;
} nikfemm_MultiLayerCurrentDensitySimulation;

static void nikfemm_MultiLayerCurrentDensitySimulation_dealloc(nikfemm_MultiLayerCurrentDensitySimulation* self) {
    delete self->simulation;
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    nikfemm_MultiLayerCurrentDensitySimulation* self;

    // constructor must have 2 arguments

    self = (nikfemm_MultiLayerCurrentDensitySimulation*)type->tp_alloc(type, 0);
    uint32_t num_layers;
    PyObject* depths;

    if (self != NULL) {
        if (PyTuple_Size(args) == 2) {
            if (!PyArg_ParseTuple(args, "IO", &num_layers, &depths)) {
                Py_XDECREF(self);
                return NULL;
            }

            if (PyList_Size(depths) != num_layers) {
                PyErr_SetString(PyExc_ValueError, "depths.size() != num_layers");
                Py_XDECREF(self);
                return NULL;
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "constructor must have 2 arguments");
            Py_XDECREF(self);
            return NULL;
        }
    }

    self->simulation = new nikfemm::MultiLayerCurrentDensitySimulation(num_layers);

    for (uint32_t i = 0; i < num_layers; i++) {
        PyObject* depth_item = PyList_GetItem(depths, i);

        if (PyFloat_Check(depth_item)) {
            self->simulation->meshes[i].depth = PyFloat_AsDouble(depth_item);
        } else if (PyLong_Check(depth_item)) {
            // convert to double
            self->simulation->meshes[i].depth = PyLong_AsDouble(depth_item);
        } else {
            PyErr_SetString(PyExc_TypeError, "depth must be a double");
            Py_XDECREF(self);
            return NULL;
        }
    }

    for (uint32_t i = 0; i < num_layers; i++) {
        printf("depth[%d] = %f\n", i, self->simulation->meshes[i].depth);
    }

    self->simulation->interconnections = std::vector<nikfemm::CurrentDensityInterconnection>();

    return (PyObject*)self;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_generateSystem(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    // generateSystem may have from 0 to 2 arguments which are: bool, double
    PyObject* refine = Py_True;
    double max_triangle_area = 1;
    double min_angle = 33;

    if (PyTuple_Size(args) > 3) {
        PyErr_SetString(PyExc_TypeError, "generateSystem must have 0 to 3 arguments: bool, double, double");
        return NULL;
    }

    switch (PyTuple_Size(args)) {
        case 3:
            if (!PyFloat_Check(PyTuple_GetItem(args, 2)) && !PyLong_Check(PyTuple_GetItem(args, 2))) {
                PyErr_SetString(PyExc_TypeError, "min_angle must be a double");
                return NULL;
            }
            min_angle = PyFloat_AsDouble(PyTuple_GetItem(args, 2));
        case 2:
            if (!PyFloat_Check(PyTuple_GetItem(args, 1)) && !PyLong_Check(PyTuple_GetItem(args, 1))) {
                PyErr_SetString(PyExc_TypeError, "max_triangle_area must be a double");
                return NULL;
            }
            max_triangle_area = PyFloat_AsDouble(PyTuple_GetItem(args, 1));
        case 1:
            if (!PyBool_Check(PyTuple_GetItem(args, 0))) {
                PyErr_SetString(PyExc_TypeError, "refine must be a boolean");
                return NULL;
            }
            refine = PyTuple_GetItem(args, 0);
        default:
            break;
    }
    auto system = self->simulation->generateSystem(PyLong_AsLong(refine), max_triangle_area, min_angle);

    nikfemm_System* py_system = (nikfemm_System*)PyObject_CallObject((PyObject*)&nikfemm_SystemType, NULL);

    if (py_system == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create nikfemm_System object");
        return NULL;
    }

    try {
        py_system->system = new nikfemm::System(system);
    } catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }

    std::vector<uint64_t> layer_offsets;
    uint64_t layer_offset = 0;
    for (auto& mesh : self->simulation->meshes) {
        layer_offsets.push_back(layer_offset);
        layer_offset += mesh.data.numberofpoints;
    }

    // return triangles and vertices too
    PyObject* triangles = PyList_New(self->simulation->meshes.size());

    for (uint64_t i = 0; i < self->simulation->meshes.size(); i++) {
        PyObject* mytriangles = PyList_New(self->simulation->meshes[i].data.numberoftriangles);
        for (uint64_t j = 0; j < self->simulation->meshes[i].data.numberoftriangles; j++) {
            PyObject* triangle = PyTuple_Pack(3, PyLong_FromUnsignedLong(self->simulation->meshes[i].data.trianglelist[j][0] + layer_offsets[i]),
                                                 PyLong_FromUnsignedLong(self->simulation->meshes[i].data.trianglelist[j][1] + layer_offsets[i]),
                                                 PyLong_FromUnsignedLong(self->simulation->meshes[i].data.trianglelist[j][2] + layer_offsets[i]));
            PyList_SetItem(mytriangles, j, triangle);
        }
        PyList_SetItem(triangles, i, mytriangles);
    }

    PyObject* vertices = PyList_New(layer_offset);

    for (uint64_t i = 0; i < self->simulation->meshes.size(); i++) {
        for (uint64_t j = 0; j < self->simulation->meshes[i].data.numberofpoints; j++) {
            PyObject* vertex = PyTuple_Pack(2, PyFloat_FromDouble(self->simulation->meshes[i].data.pointlist[j].x),
                                              PyFloat_FromDouble(self->simulation->meshes[i].data.pointlist[j].y));
            PyList_SetItem(vertices, j + layer_offsets[i], vertex);
        }
    }

    PyObject* result = PyTuple_Pack(3, (PyObject*)py_system, triangles, vertices);

    return (PyObject*)result;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_solve(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    if (PyTuple_Size(args) != 1) {
        PyErr_SetString(PyExc_TypeError, "solve must have 1 argument of type System");
        return NULL;
    } else {
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 0), &nikfemm_SystemType)) {
            PyErr_SetString(PyExc_TypeError, "solve must have 1 argument of type System");
            return NULL;
        }
    }

    nikfemm::System<double>* system = ((nikfemm_System*)PyTuple_GetItem(args, 0))->system;
    std::vector<double> V = self->simulation->solve(*system);
    std::vector<double> P = self->simulation->computePowerDensity(V);

    // create list of tuples like this: [(int, double, double, double), ...]
    // where the first element is the layer id, and the next three are the x, y, and voltage

    PyObject* voltages = PyList_New(V.size());

    std::vector<uint64_t> layer_offsets;
    uint64_t layer_offset = 0;
    for (auto& mesh : self->simulation->meshes) {
        layer_offsets.push_back(layer_offset);
        layer_offset += mesh.data.numberofpoints;
    }

    for (uint64_t i = 0; i < self->simulation->meshes.size(); i++) {
        for (uint64_t j = 0; j < self->simulation->meshes[i].data.numberofpoints; j++) {
            PyObject* voltage = PyTuple_Pack(4, PyLong_FromUnsignedLong(i), 
                                                PyFloat_FromDouble(self->simulation->meshes[i].data.pointlist[j].x),
                                                PyFloat_FromDouble(self->simulation->meshes[i].data.pointlist[j].y), 
                                                PyFloat_FromDouble(V[j + layer_offsets[i]]));
            PyList_SetItem(voltages, j + layer_offsets[i], voltage);
        }
    }

    // create list of tuples like this: [(int, double), ...]
    // where the first element is the layer id, and the next is the power density

    PyObject* power_densities = PyList_New(P.size());

    for (uint64_t i = 0; i < P.size(); i++) {
        PyObject* power_density = PyTuple_Pack(2, PyLong_FromUnsignedLong(i), PyFloat_FromDouble(P[i]));
        PyList_SetItem(power_densities, i, power_density);
    }

    // create list of tuples like this with the triangles: [(int, int, int, int), ...]
    // where the first element is the layer id, and the next three are the indices of the vertices

    uint64_t total_number_of_triangles = 0;
    std::vector<uint64_t> triangle_offsets;
    for (auto& mesh : self->simulation->meshes) {
        triangle_offsets.push_back(total_number_of_triangles);
        total_number_of_triangles += mesh.data.numberoftriangles;
    }

    PyObject* triangles = PyList_New(total_number_of_triangles);

    for (uint64_t i = 0; i < self->simulation->meshes.size(); i++) {
        for (uint64_t j = 0; j < self->simulation->meshes[i].data.numberoftriangles; j++) {
            PyObject* triangle = PyTuple_Pack(4, PyLong_FromUnsignedLong(i), 
                                                 PyLong_FromUnsignedLong(self->simulation->meshes[i].data.trianglelist[j][0] + layer_offsets[i]),
                                                 PyLong_FromUnsignedLong(self->simulation->meshes[i].data.trianglelist[j][1] + layer_offsets[i]),
                                                 PyLong_FromUnsignedLong(self->simulation->meshes[i].data.trianglelist[j][2] + layer_offsets[i])
            );
            PyList_SetItem(triangles, j + triangle_offsets[i], triangle);
        }
    }

    PyObject* result = PyTuple_Pack(3, voltages, triangles, power_densities);

    return result;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_setVoltage(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    if (PyTuple_Size(args) != 4) {
        PyErr_SetString(PyExc_TypeError, "setVoltage must have 4 arguments: System, [double, double], double, uint64_t");
        return NULL;
    } else {
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 0), &nikfemm_SystemType)) {
            PyErr_SetString(PyExc_TypeError, "setVoltage must have the first argument of type System");
            return NULL;
        }
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 1), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "setVoltage must have the second argument of type [double, double]");
            return NULL;
        }
        if (!PyFloat_Check(PyTuple_GetItem(args, 2)) && !PyLong_Check(PyTuple_GetItem(args, 2))) {
            PyErr_SetString(PyExc_TypeError, "setVoltage must have the third argument of type double");
            return NULL;
        }
        if (!PyLong_Check(PyTuple_GetItem(args, 3))) {
            PyErr_SetString(PyExc_TypeError, "setVoltage must have the fourth argument of type uint64_t");
            return NULL;
        }
    }

    nikfemm::System<double>* system = ((nikfemm_System*)PyTuple_GetItem(args, 0))->system;
    double coords[2];
    for (int i = 0; i < 2; i++) {
        PyObject* item = PyList_GetItem(PyTuple_GetItem(args, 1), i);
        if (PyFloat_Check(item)) {
            coords[i] = PyFloat_AsDouble(item);
        } else if (PyLong_Check(item)) {
            coords[i] = PyLong_AsDouble(item);
        } else {
            PyErr_SetString(PyExc_TypeError, "setVoltage must have the second argument of type [double, double]");
            return NULL;
        }
    }
    nikfemm::Vector p = nikfemm::Vector(coords[0], coords[1]);

    double V;
    if (PyFloat_Check(PyTuple_GetItem(args, 2))) {
        V = PyFloat_AsDouble(PyTuple_GetItem(args, 2));
    } else if (PyLong_Check(PyTuple_GetItem(args, 2))) {
        V = PyLong_AsDouble(PyTuple_GetItem(args, 2));
    } else {
        PyErr_SetString(PyExc_TypeError, "setVoltage must have the third argument of type double");
        return NULL;
    }
    uint64_t layer_id = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 3));

    self->simulation->setVoltage(*system, p, V, layer_id);

    // printf("setVoltage: p = (%f, %f), V = %f, layer_id = %lu\n", p.x, p.y, V, layer_id);

    Py_RETURN_NONE;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_add_interconnection(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    // arguments should be [double, double], [double, double], uint64_t, uint64_t, double
    // the first two arguments are lists of doubles
    if (PyTuple_Size(args) != 5) {
        PyErr_SetString(PyExc_TypeError, "add_interconnection must have 5 arguments: [double, double], [double, double], uint64_t, uint64_t, double");
        return NULL;
    } else {
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 0), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "add_interconnection must have the first argument of type [double, double] or [int, int]");
            return NULL;
        }
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 1), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "add_interconnection must have the second argument of type [double, double] or [int, int]");
            return NULL;
        }
        if (!PyLong_Check(PyTuple_GetItem(args, 2))) {
            PyErr_SetString(PyExc_TypeError, "add_interconnection must have the third argument of type uint64_t");
            return NULL;
        }
        if (!PyLong_Check(PyTuple_GetItem(args, 3))) {
            PyErr_SetString(PyExc_TypeError, "add_interconnection must have the fourth argument of type uint64_t");
            return NULL;
        }
        if (!PyFloat_Check(PyTuple_GetItem(args, 4)) && !PyLong_Check(PyTuple_GetItem(args, 4))) {
            PyErr_SetString(PyExc_TypeError, "add_interconnection must have the fifth argument of type double");
            return NULL;
        }
    }

    // extract arguments
    // the first two arguments are lists of doubles, convert them to Vector
    double coords[2][2];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            PyObject* item = PyList_GetItem(PyTuple_GetItem(args, i), j);
            if (PyFloat_Check(item)) {
                coords[i][j] = PyFloat_AsDouble(item);
            } else if (PyLong_Check(item)) {
                coords[i][j] = PyLong_AsDouble(item);
            } else {
                PyErr_SetString(PyExc_TypeError, "add_interconnection must have the first two arguments of type [double, double]");
                return NULL;
            }
        }
    }
    // for each item if it is a long convert it to double, if it is a double use it
    nikfemm::Vector p1 = nikfemm::Vector(coords[0][0], coords[0][1]);
    nikfemm::Vector p2 = nikfemm::Vector(coords[1][0], coords[1][1]);
    uint64_t layer1_id = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 2));
    uint64_t layer2_id = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 3));
    double R;
    if (PyFloat_Check(PyTuple_GetItem(args, 4))) {
        R = PyFloat_AsDouble(PyTuple_GetItem(args, 4));
    } else if (PyLong_Check(PyTuple_GetItem(args, 4))) {
        R = PyLong_AsDouble(PyTuple_GetItem(args, 4));
    } else {
        PyErr_SetString(PyExc_TypeError, "add_interconnection must have the fifth argument of type double");
        return NULL;
    }

    // printf("add_interconnection: p1 = (%f, %f), p2 = (%f, %f), layer1_id = %lu, layer2_id = %lu, R = %f\n", p1.x, p1.y, p2.x, p2.y, layer1_id, layer2_id, R);

    // add interconnection
    self->simulation->interconnections.push_back({p1, p2, layer1_id, layer2_id, R});

    Py_RETURN_NONE;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_drawRectangle(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    // arguments should be (double, double), (double, double), uint64_t
    // the first two arguments are tuples of doubles
    if (PyTuple_Size(args) != 3) {
        PyErr_SetString(PyExc_TypeError, "drawRectangle must have 3 arguments: (double, double), (double, double), uint64_t");
        return NULL;
    } else {
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 0), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "drawRectangle must have the first argument of type (double, double)");
            return NULL;
        }
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 1), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "drawRectangle must have the second argument of type (double, double)");
            return NULL;
        }
        if (!PyLong_Check(PyTuple_GetItem(args, 2))) {
            PyErr_SetString(PyExc_TypeError, "drawRectangle must have the third argument of type uint64_t");
            return NULL;
        }
    }

    // extract arguments
    // the first two arguments are tuples of doubles, convert them to Vector
    double coords[2][2];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            PyObject* item = PyList_GetItem(PyTuple_GetItem(args, i), j);
            if (PyFloat_Check(item)) {
                coords[i][j] = PyFloat_AsDouble(item);
            } else if (PyLong_Check(item)) {
                coords[i][j] = PyLong_AsDouble(item);
            } else {
                PyErr_SetString(PyExc_TypeError, "drawRectangle must have the first two arguments of type (double, double)");
                return NULL;
            }
        }
    }
    nikfemm::Vector p1 = nikfemm::Vector(coords[0][0], coords[0][1]);
    nikfemm::Vector p2 = nikfemm::Vector(coords[1][0], coords[1][1]);
    uint64_t layer_id = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 2));

    // draw rectangle
    self->simulation->meshes[layer_id].drawing.drawRectangle(p1, p2);

    // printf("drawRectangle: p1 = (%f, %f), p2 = (%f, %f), layer_id = %lu\n", p1.x, p1.y, p2.x, p2.y, layer_id);

    Py_RETURN_NONE;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_drawRegion(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    // arguments should be (double, double), double, uint64_t
    if (PyTuple_Size(args) != 3) {
        PyErr_SetString(PyExc_TypeError, "drawRegion must have 3 arguments: (double, double), double, uint64_t");
        return NULL;
    } else {
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 0), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "drawRegion must have the first argument of type (double, double)");
            return NULL;
        }
        if (!PyFloat_Check(PyTuple_GetItem(args, 1))) {
            PyErr_SetString(PyExc_TypeError, "drawRegion must have the second argument of type double");
            return NULL;
        }
        if (!PyLong_Check(PyTuple_GetItem(args, 2))) {
            PyErr_SetString(PyExc_TypeError, "drawRegion must have the third argument of type uint64_t");
            return NULL;
        }
    }

    // extract arguments
    // the first argument is a tuple of doubles, convert it to Vector
    double coords[2];
    for (int i = 0; i < 2; i++) {
        PyObject* item = PyList_GetItem(PyTuple_GetItem(args, 0), i);
        if (PyFloat_Check(item)) {
            coords[i] = PyFloat_AsDouble(item);
        } else if (PyLong_Check(item)) {
            coords[i] = PyLong_AsDouble(item);
        } else {
            PyErr_SetString(PyExc_TypeError, "drawRegion must have the first argument of type (double, double)");
            return NULL;
        }
    }
    nikfemm::Vector p = nikfemm::Vector(coords[0], coords[1]);
    double material = PyFloat_AsDouble(PyTuple_GetItem(args, 1));
    uint64_t layer_id = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 2));

    // draw region
    self->simulation->meshes[layer_id].drawing.drawRegion(p, {material});

    // printf("drawRegion: p = (%f, %f), material = %f, layer_id = %lu\n", p.x, p.y, material, layer_id);

    Py_RETURN_NONE;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_drawPolygon(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    // arguments should be list of (double, double), uint64_t
    if (PyTuple_Size(args) != 2) {
        PyErr_SetString(PyExc_TypeError, "drawPolygon must have 2 arguments: [[double, double], ...], uint64_t");
        return NULL;
    } else {
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 0), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "drawPolygon must have the first argument of type [[double, double], ...]");
            return NULL;
        }
        if (!PyLong_Check(PyTuple_GetItem(args, 1))) {
            PyErr_SetString(PyExc_TypeError, "drawPolygon must have the second argument of type uint64_t");
            return NULL;
        }
    }

    // extract arguments
    // the first argument is a list of lists of two doubles, convert it to vector of Vector
    std::vector<nikfemm::Vector> points;
    for (int i = 0; i < PyList_Size(PyTuple_GetItem(args, 0)); i++) {
        PyObject* item = PyList_GetItem(PyTuple_GetItem(args, 0), i);
        if (!PyObject_TypeCheck(item, &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "drawPolygon must have the first argument of type [[double, double], ...]");
            return NULL;
        }
        if (PyList_Size(item) != 2) {
            PyErr_SetString(PyExc_TypeError, "drawPolygon must have the first argument of type [[double, double], ...]");
            return NULL;
        }
        double coords[2];
        for (int j = 0; j < 2; j++) {
            PyObject* subitem = PyList_GetItem(item, j);
            if (PyFloat_Check(subitem)) {
                coords[j] = PyFloat_AsDouble(subitem);
            } else if (PyLong_Check(subitem)) {
                coords[j] = PyLong_AsDouble(subitem);
            } else {
                PyErr_SetString(PyExc_TypeError, "drawPolygon must have the first argument of type [[double, double], ...]");
                return NULL;
            }
        }
        points.push_back(nikfemm::Vector(coords[0], coords[1]));
    }
    uint64_t layer_id = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 1));

    // draw polygon
    self->simulation->meshes[layer_id].drawing.drawPolygon(points);

    // printf("drawPolygon: layer_id = %lu\n", layer_id);

    Py_RETURN_NONE;
}

static PyMethodDef nikfemm_MultiLayerCurrentDensitySimulation_methods[] = {
    {"generate_system", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_generateSystem, METH_VARARGS, "Generate system of linear equations for current density"},
    {"solve", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_solve, METH_VARARGS, "Solve system of linear equations for current density"},
    {"set_voltage", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_setVoltage, METH_VARARGS, "Set voltage on a node"},
    {"add_interconnection", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_add_interconnection, METH_VARARGS, "Add interconnection between two points in two layers with a resistance"},
    {"draw_rectangle", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_drawRectangle, METH_VARARGS, "Draw a rectangle on a layer"},
    {"draw_region", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_drawRegion, METH_VARARGS, "Draw a region on a layer"},
    {"draw_polygon", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_drawPolygon, METH_VARARGS, "Draw a polygon on a layer"},
    {NULL} /* Sentinel */
};

static PyTypeObject nikfemm_MultiLayerCurrentDensitySimulationType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "nikfemm.MultiLayerCurrentDensitySimulation",
    .tp_basicsize = sizeof(nikfemm_MultiLayerCurrentDensitySimulation),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)nikfemm_MultiLayerCurrentDensitySimulation_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_doc = PyDoc_STR("Multi-layer current density simulation"),
    .tp_methods = nikfemm_MultiLayerCurrentDensitySimulation_methods,
    .tp_new = nikfemm_MultiLayerCurrentDensitySimulation_new,
};

static PyModuleDef nikfemm_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "nikfemm",
    .m_doc = "Python bindings for nikfemm",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_nikfemm(void) {
    PyObject* m;

    if (PyType_Ready(&nikfemm_SystemType) < 0) {
        return NULL;
    }

    if (PyType_Ready(&nikfemm_MultiLayerCurrentDensitySimulationType) < 0) {
        return NULL;
    }

    m = PyModule_Create(&nikfemm_module);
    if (m == NULL) {
        return NULL;
    }

    Py_INCREF(&nikfemm_SystemType);
    if (PyModule_AddObject(m, "System", (PyObject*)&nikfemm_SystemType) < 0) {
        Py_DECREF(&nikfemm_SystemType);
        Py_DECREF(m);
        return NULL;
    }

    Py_INCREF(&nikfemm_MultiLayerCurrentDensitySimulationType);
    if (PyModule_AddObject(m, "MultiLayerCurrentDensitySimulation", (PyObject*)&nikfemm_MultiLayerCurrentDensitySimulationType) < 0) {
        Py_DECREF(&nikfemm_MultiLayerCurrentDensitySimulationType);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}