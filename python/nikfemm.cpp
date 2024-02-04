#include <iostream>
#include <string>
#include <functional>
#include <iomanip>

#include <nikfemm.hpp>

#include <Python.h>

/*
namespace nikfemm {
    struct CurrentDensitySystem {
        BuildMatCOO<double> A;
        CV b;

        void addDirichletBoundaryCondition(uint32_t id, double value) {
            // https://community.freefem.org/t/implementation-of-dirichlet-boundary-condition-when-tgv-1/113
            // this function lets you set a Dirichlet boundary condition on a node

            // every element that has a row or column index equal to id is set to zero
            // for every element that has a column index equal to id, its value multiplied 
            // by the value of the boundary condition is subtracted from the corresponding element in the b vector
            // care must be taken because the matrix is symmetric but it is stored as upper triangular to save memory

            for (auto& elem : A.elems) {
                uint32_t m = elem.first >> 32;
                uint32_t n = elem.first & 0xFFFFFFFF;

                // skip diagonal elements
                if (m == n) continue;

                if (m == id) { // row index equal to id
                    b.val[n] -= elem.second * value;
                    elem.second = 0;
                }
                if (n == id) { // column index equal to id
                    b.val[m] -= elem.second * value;
                    elem.second = 0;
                }
            }

            // coo.elems[BuildMatCOO<int>::get_key(id, id)].setToConstant(1);
            A(id, id) = 1;
            b.val[id] = value;
        }
    };
}
*/

typedef struct {
    PyObject_HEAD
    nikfemm::CurrentDensitySystem* system;
} nikfemm_CurrentDensitySystem;

static void nikfemm_CurrentDensitySystem_dealloc(nikfemm_CurrentDensitySystem* self) {
    delete self->system;
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* nikfemm_CurrentDensitySystem_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    nikfemm_CurrentDensitySystem* self;

    self = (nikfemm_CurrentDensitySystem*)type->tp_alloc(type, 0);

    return (PyObject*)self;
}

static PyTypeObject nikfemm_CurrentDensitySystemType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "nikfemm.CurrentDensitySystem",
    .tp_basicsize = sizeof(nikfemm_CurrentDensitySystem),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)nikfemm_CurrentDensitySystem_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_doc = PyDoc_STR("System of linear equations for current density"),
    .tp_new = nikfemm_CurrentDensitySystem_new,
};

/*
namespace nikfemm {
    class MultiLayerCurrentDensitySimulation {
        public:
            std::vector<CurrentDensityMesh> meshes;
            CV V;
            std::vector<CurrentDensityInterconnection> interconnections;

            MultiLayerCurrentDensitySimulation(uint32_t num_layers, std::vector<double> depths, std::vector<double> max_triangle_areas);
            MultiLayerCurrentDensitySimulation(uint32_t num_layers);
            ~MultiLayerCurrentDensitySimulation();

            CurrentDensitySystem generateSystem(bool refine = true);
            void solve(CurrentDensitySystem& system);
            void setVoltage(CurrentDensitySystem& system, Vector p, double V, uint64_t layer_id);
        protected:
#ifdef NIKFEMM_USE_OPENCV
            void VplotRend(cv::Mat* image, double width, double height, uint64_t layer_id, double maxV, double minV);
        public:
            void Vplot(uint32_t width, uint32_t height);
#endif
    };
}

namespace nikfemm {
    MultiLayerCurrentDensitySimulation::MultiLayerCurrentDensitySimulation(uint32_t num_layers, std::vector<double> depths, std::vector<double> max_triangle_areas) {
        if (depths.size() != num_layers) {
            throw std::invalid_argument("depths.size() != num_layers");
        }
        if (max_triangle_areas.size() != num_layers) {
            throw std::invalid_argument("max_triangle_areas.size() != num_layers");
        }

        for (uint32_t i = 0; i < num_layers; i++) {
            meshes.push_back(CurrentDensityMesh());
            meshes[i].max_triangle_area = max_triangle_areas[i];
            meshes[i].depth = depths[i];
        }
    }

    MultiLayerCurrentDensitySimulation::MultiLayerCurrentDensitySimulation(uint32_t num_layers) {
        for (uint32_t i = 0; i < num_layers; i++) {
            meshes.push_back(CurrentDensityMesh());
            meshes[i].max_triangle_area = 1;
            meshes[i].depth = 1;
        }
    }

    MultiLayerCurrentDensitySimulation::~MultiLayerCurrentDensitySimulation() {

    }
}
*/

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

    // constructor must have 3 arguments

    self = (nikfemm_MultiLayerCurrentDensitySimulation*)type->tp_alloc(type, 0);
    uint32_t num_layers;
    PyObject* depths;
    PyObject* max_triangle_areas;
    if (self != NULL) {
        if (PyTuple_Size(args) == 3) {
            if (!PyArg_ParseTuple(args, "IOO", &num_layers, &depths, &max_triangle_areas)) {
                Py_XDECREF(self);
                return NULL;
            }

            if (!PyList_Check(depths) || !PyList_Check(max_triangle_areas)) {
                PyErr_SetString(PyExc_TypeError, "depths and max_triangle_areas must be lists");
                Py_XDECREF(self);
                return NULL;
            }

            if (PyList_Size(depths) != num_layers) {
                PyErr_SetString(PyExc_ValueError, "depths.size() != num_layers");
                Py_XDECREF(self);
                return NULL;
            }

            if (PyList_Size(max_triangle_areas) != num_layers) {
                PyErr_SetString(PyExc_ValueError, "max_triangle_areas.size() != num_layers");
                Py_XDECREF(self);
                return NULL;
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "constructor must have 3 arguments");
            Py_XDECREF(self);
            return NULL;
        }
    }

    self->simulation = new nikfemm::MultiLayerCurrentDensitySimulation(num_layers);

    for (uint32_t i = 0; i < num_layers; i++) {
        PyObject* max_triangle_area_item = PyList_GetItem(max_triangle_areas, i);
        PyObject* depth_item = PyList_GetItem(depths, i);

        // if floats use them
        if (PyFloat_Check(max_triangle_area_item)) {
            self->simulation->meshes[i].max_triangle_area = PyFloat_AsDouble(max_triangle_area_item);
        } else if (PyLong_Check(max_triangle_area_item)) {
            // convert to float
            self->simulation->meshes[i].max_triangle_area = PyLong_AsDouble(max_triangle_area_item);
        } else {
            PyErr_SetString(PyExc_TypeError, "max_triangle_area must be a float");
            Py_XDECREF(self);
            return NULL;
        }

        if (PyFloat_Check(depth_item)) {
            self->simulation->meshes[i].depth = PyFloat_AsDouble(depth_item);
        } else if (PyLong_Check(depth_item)) {
            // convert to float
            self->simulation->meshes[i].depth = PyLong_AsDouble(depth_item);
        } else {
            PyErr_SetString(PyExc_TypeError, "depth must be a float");
            Py_XDECREF(self);
            return NULL;
        }
    }

    for (uint32_t i = 0; i < num_layers; i++) {
        printf("max_triangle_area[%d] = %f\n", i, self->simulation->meshes[i].max_triangle_area);
        printf("depth[%d] = %f\n", i, self->simulation->meshes[i].depth);
    }

    self->simulation->interconnections = std::vector<nikfemm::CurrentDensityInterconnection>();

    return (PyObject*)self;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_generateSystem(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    // generateSystem may have 0 or 1 arguments
    PyObject* refine = Py_True;

    if (PyTuple_Size(args) == 0) {
        // do nothing, default value is already set
    } else if (PyTuple_Size(args) == 1) {
        // check type
        if (!PyBool_Check(PyTuple_GetItem(args, 0))) {
            PyErr_SetString(PyExc_TypeError, "refine must be a boolean");
            return NULL;
        }

        refine = PyTuple_GetItem(args, 0);
    } else if (PyTuple_Size(args) > 1) {
        PyErr_SetString(PyExc_TypeError, "generateSystem may have 0 or 1 arguments");
        return NULL;
    }

    auto system = self->simulation->generateSystem(PyLong_AsLong(refine));
    
    nikfemm_CurrentDensitySystem* py_system = (nikfemm_CurrentDensitySystem*)PyObject_CallObject((PyObject*)&nikfemm_CurrentDensitySystemType, NULL);

    if (py_system == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create nikfemm_CurrentDensitySystem object");
        return NULL;
    }

    try {
        py_system->system = new nikfemm::CurrentDensitySystem(system);
    } catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }

    return (PyObject*)py_system;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_solve(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    if (PyTuple_Size(args) != 1) {
        PyErr_SetString(PyExc_TypeError, "solve must have 1 argument of type CurrentDensitySystem");
        return NULL;
    } else {
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 0), &nikfemm_CurrentDensitySystemType)) {
            PyErr_SetString(PyExc_TypeError, "solve must have 1 argument of type CurrentDensitySystem");
            return NULL;
        }
    }

    nikfemm::CurrentDensitySystem* system = ((nikfemm_CurrentDensitySystem*)PyTuple_GetItem(args, 0))->system;
    self->simulation->solve(*system);

    Py_RETURN_NONE;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_setVoltage(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    if (PyTuple_Size(args) != 4) {
        PyErr_SetString(PyExc_TypeError, "setVoltage must have 4 arguments: CurrentDensitySystem, [float, float], float, uint64_t");
        return NULL;
    } else {
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 0), &nikfemm_CurrentDensitySystemType)) {
            PyErr_SetString(PyExc_TypeError, "setVoltage must have the first argument of type CurrentDensitySystem");
            return NULL;
        }
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 1), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "setVoltage must have the second argument of type [float, float]");
            return NULL;
        }
        if (!PyFloat_Check(PyTuple_GetItem(args, 2)) && !PyLong_Check(PyTuple_GetItem(args, 2))) {
            PyErr_SetString(PyExc_TypeError, "setVoltage must have the third argument of type float");
            return NULL;
        }
        if (!PyLong_Check(PyTuple_GetItem(args, 3))) {
            PyErr_SetString(PyExc_TypeError, "setVoltage must have the fourth argument of type uint64_t");
            return NULL;
        }
    }

    nikfemm::CurrentDensitySystem* system = ((nikfemm_CurrentDensitySystem*)PyTuple_GetItem(args, 0))->system;
    double coords[2];
    for (int i = 0; i < 2; i++) {
        PyObject* item = PyList_GetItem(PyTuple_GetItem(args, 1), i);
        if (PyFloat_Check(item)) {
            coords[i] = PyFloat_AsDouble(item);
        } else if (PyLong_Check(item)) {
            coords[i] = PyLong_AsDouble(item);
        } else {
            PyErr_SetString(PyExc_TypeError, "setVoltage must have the second argument of type [float, float]");
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
        PyErr_SetString(PyExc_TypeError, "setVoltage must have the third argument of type float");
        return NULL;
    }
    uint64_t layer_id = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 3));

    self->simulation->setVoltage(*system, p, V, layer_id);

    printf("setVoltage: p = (%f, %f), V = %f, layer_id = %lu\n", p.x, p.y, V, layer_id);

    Py_RETURN_NONE;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_add_interconnection(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    // arguments should be [float, float], [float, float], uint64_t, uint64_t, float
    // the first two arguments are lists of floats
    if (PyTuple_Size(args) != 5) {
        PyErr_SetString(PyExc_TypeError, "add_interconnection must have 5 arguments: [float, float], [float, float], uint64_t, uint64_t, float");
        return NULL;
    } else {
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 0), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "add_interconnection must have the first argument of type [float, float] or [int, int]");
            return NULL;
        }
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 1), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "add_interconnection must have the second argument of type [float, float] or [int, int]");
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
            PyErr_SetString(PyExc_TypeError, "add_interconnection must have the fifth argument of type float");
            return NULL;
        }
    }

    // extract arguments
    // the first two arguments are lists of floats, convert them to Vector
    double coords[2][2];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            PyObject* item = PyList_GetItem(PyTuple_GetItem(args, i), j);
            if (PyFloat_Check(item)) {
                coords[i][j] = PyFloat_AsDouble(item);
            } else if (PyLong_Check(item)) {
                coords[i][j] = PyLong_AsDouble(item);
            } else {
                PyErr_SetString(PyExc_TypeError, "add_interconnection must have the first two arguments of type [float, float]");
                return NULL;
            }
        }
    }
    // for each item if it is a long convert it to float, if it is a float use it
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
        PyErr_SetString(PyExc_TypeError, "add_interconnection must have the fifth argument of type float");
        return NULL;
    }

    printf("add_interconnection: p1 = (%f, %f), p2 = (%f, %f), layer1_id = %lu, layer2_id = %lu, R = %f\n", p1.x, p1.y, p2.x, p2.y, layer1_id, layer2_id, R);

    // add interconnection
    self->simulation->interconnections.push_back({p1, p2, layer1_id, layer2_id, R});

    Py_RETURN_NONE;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_drawRectangle(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    // arguments should be (float, float), (float, float), uint64_t
    // the first two arguments are tuples of floats
    if (PyTuple_Size(args) != 3) {
        PyErr_SetString(PyExc_TypeError, "drawRectangle must have 3 arguments: (float, float), (float, float), uint64_t");
        return NULL;
    } else {
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 0), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "drawRectangle must have the first argument of type (float, float)");
            return NULL;
        }
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 1), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "drawRectangle must have the second argument of type (float, float)");
            return NULL;
        }
        if (!PyLong_Check(PyTuple_GetItem(args, 2))) {
            PyErr_SetString(PyExc_TypeError, "drawRectangle must have the third argument of type uint64_t");
            return NULL;
        }
    }

    // extract arguments
    // the first two arguments are tuples of floats, convert them to Vector
    double coords[2][2];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            PyObject* item = PyList_GetItem(PyTuple_GetItem(args, i), j);
            if (PyFloat_Check(item)) {
                coords[i][j] = PyFloat_AsDouble(item);
            } else if (PyLong_Check(item)) {
                coords[i][j] = PyLong_AsDouble(item);
            } else {
                PyErr_SetString(PyExc_TypeError, "drawRectangle must have the first two arguments of type (float, float)");
                return NULL;
            }
        }
    }
    nikfemm::Vector p1 = nikfemm::Vector(coords[0][0], coords[0][1]);
    nikfemm::Vector p2 = nikfemm::Vector(coords[1][0], coords[1][1]);
    uint64_t layer_id = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 2));

    // draw rectangle
    self->simulation->meshes[layer_id].drawing.drawRectangle(p1, p2);

    printf("drawRectangle: p1 = (%f, %f), p2 = (%f, %f), layer_id = %lu\n", p1.x, p1.y, p2.x, p2.y, layer_id);

    Py_RETURN_NONE;
}

static PyObject* nikfemm_MultiLayerCurrentDensitySimulation_drawRegion(nikfemm_MultiLayerCurrentDensitySimulation* self, PyObject* args, PyObject* kwds) {
    // arguments should be (float, float), float, uint64_t
    if (PyTuple_Size(args) != 3) {
        PyErr_SetString(PyExc_TypeError, "drawRegion must have 3 arguments: (float, float), float, uint64_t");
        return NULL;
    } else {
        if (!PyObject_TypeCheck(PyTuple_GetItem(args, 0), &PyList_Type)) {
            PyErr_SetString(PyExc_TypeError, "drawRegion must have the first argument of type (float, float)");
            return NULL;
        }
        if (!PyFloat_Check(PyTuple_GetItem(args, 1))) {
            PyErr_SetString(PyExc_TypeError, "drawRegion must have the second argument of type float");
            return NULL;
        }
        if (!PyLong_Check(PyTuple_GetItem(args, 2))) {
            PyErr_SetString(PyExc_TypeError, "drawRegion must have the third argument of type uint64_t");
            return NULL;
        }
    }

    // extract arguments
    // the first argument is a tuple of floats, convert it to Vector
    double coords[2];
    for (int i = 0; i < 2; i++) {
        PyObject* item = PyList_GetItem(PyTuple_GetItem(args, 0), i);
        if (PyFloat_Check(item)) {
            coords[i] = PyFloat_AsDouble(item);
        } else if (PyLong_Check(item)) {
            coords[i] = PyLong_AsDouble(item);
        } else {
            PyErr_SetString(PyExc_TypeError, "drawRegion must have the first argument of type (float, float)");
            return NULL;
        }
    }
    nikfemm::Vector p = nikfemm::Vector(coords[0], coords[1]);
    double material = PyFloat_AsDouble(PyTuple_GetItem(args, 1));
    uint64_t layer_id = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 2));

    // draw region
    self->simulation->meshes[layer_id].drawing.drawRegion(p, {(float)material});

    printf("drawRegion: p = (%f, %f), material = %f, layer_id = %lu\n", p.x, p.y, material, layer_id);

    Py_RETURN_NONE;
}

static PyMethodDef nikfemm_MultiLayerCurrentDensitySimulation_methods[] = {
    {"generate_system", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_generateSystem, METH_VARARGS, "Generate system of linear equations for current density"},
    {"solve", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_solve, METH_VARARGS, "Solve system of linear equations for current density"},
    {"set_voltage", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_setVoltage, METH_VARARGS, "Set voltage on a node"},
    {"add_interconnection", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_add_interconnection, METH_VARARGS, "Add interconnection between two points in two layers with a resistance"},
    {"draw_rectangle", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_drawRectangle, METH_VARARGS, "Draw a rectangle on a layer"},
    {"draw_region", (PyCFunction)nikfemm_MultiLayerCurrentDensitySimulation_drawRegion, METH_VARARGS, "Draw a region on a layer"},
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

    if (PyType_Ready(&nikfemm_CurrentDensitySystemType) < 0) {
        return NULL;
    }

    if (PyType_Ready(&nikfemm_MultiLayerCurrentDensitySimulationType) < 0) {
        return NULL;
    }

    m = PyModule_Create(&nikfemm_module);
    if (m == NULL) {
        return NULL;
    }

    Py_INCREF(&nikfemm_CurrentDensitySystemType);
    if (PyModule_AddObject(m, "CurrentDensitySystem", (PyObject*)&nikfemm_CurrentDensitySystemType) < 0) {
        Py_DECREF(&nikfemm_CurrentDensitySystemType);
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