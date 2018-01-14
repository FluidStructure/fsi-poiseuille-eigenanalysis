//==============|
//  NAME        : velocitymodule.cpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 08.06.2010
//  DESCRIPTION : Py interface to the Fmm lambVortex, vortexelement and semi-infinite vortex sheet fmm functions. a fast multipole algorithm adapted from (Greengard, 1987) fast alg for coulombic particle systems
//  NOTES       : .
//  TODO        :
//==============|

#include <Python.h>
#include "treesolverlib.hpp"

PyObject * cfmm2py(double * contmem, int numTargets);
int cfmm2py(double * contmem, PyObject * outsequence, int numTargets);
int cfmm2pylist(double * contmem, PyObject * outsequence, int numTargets);
double * py2cfmm(PyObject * inobj);

static PyObject * lambVortex(PyObject*, PyObject *args)
{
    PyObject * py_x, * py_y, * py_str, * py_coresqr, * py_ex, * py_ey;
    PyObject * py_out_u, * py_out_v, * out_object;
    double precision;
    int maxThreads, targetsPerBox;
    // All inputs are in the tuple, need to get them out.
    bool no_errors = true;
    if (!PyArg_ParseTuple(args, "OOOOOOOOidi", &py_x, &py_y, &py_str, &py_coresqr, &py_ex, &py_ey, &py_out_u, &py_out_v, &maxThreads, &precision, &targetsPerBox)) { return NULL; };
    // as we know that input are lists or arrays need to get the data out;
    //is a tuple containing many tupples.
    // create temp c-array for input to function
    double * x, * y, * str, * coresqr, * ex, * ey, * out_u, * out_v;
    double * str_2PI, * core;
    int numTargets, numParticles;

    // pass pyobject and array reference in and perform conversion.
    // setting up output arrays first as this is easy.
    numParticles = PyObject_Length(py_x);
    numTargets = PyObject_Length(py_ex);
    if (numParticles < 1 || numTargets < 1) {no_errors = false;}

    // now targets
    // first checking they are all the same length
    if (PyObject_Length(py_ey) != numTargets) {no_errors = false;}
    if (PyObject_Length(py_y) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_str) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_coresqr) != numParticles) {no_errors = false;}

    if (no_errors)
    {
        // for particles first
        x = py2cfmm(py_x);
        y = py2cfmm(py_y);
        str = py2cfmm(py_str);
        core = py2cfmm(py_coresqr);
        // for targets now
        ex = py2cfmm(py_ex);
        ey = py2cfmm(py_ey);
        // for output fields
        out_u = new double[numTargets];
        out_v = new double[numTargets];
        // Do the actual computations using my c function
        // Cfmm requires inputs as a continuous c-array block of memory...thats
        //why all this conversion was necessary its format is:
        //static void run( double accuracy, size_t targetsPerBox, size_t
        //maxThreads, double * x, double * y,double * str, double * core,
        //double * evalx, double * evaly,  double * evalu, double * evalv,
        //size_t numparticles, size_t numtargets);
        double twoPi = 2.0 * M_PI;

        // Need to preprocess with coresquared and /2pi
        for (int index = 0; index < numParticles; index++)
        {
            core[index] *= core[index];
            str[index] /= twoPi;
        }
        coresqr = core;
        str_2PI = str;

        //timeval start,stop; gettimeofday(&start, NULL);
        bool withFMM = true;
        if (precision <= 0.0)
        {
            withFMM = false;
        }

        LambVortexFMM fmm;
        fmm.run( precision, targetsPerBox, maxThreads, x, y, str_2PI, coresqr, ex, ey, out_u, out_v, numParticles , numTargets, withFMM);
        //gettimeofday(&stop, NULL); double timeFull = timeDiff(start, stop);
        //cout << "full time taken " << timeFull << endl;

        // Retruning calculated values back to python now
        if ((cfmm2py(out_u, py_out_u, numTargets))) { no_errors = false; }
        if ((cfmm2py(out_v, py_out_v, numTargets))) { no_errors = false; }
         
        // clean up mess
        // void PyBuffer_Release(PyObject *obj, Py_buffer *view)¶
        delete out_u;
        delete out_v;
        delete x;
        delete y;
        delete str;
        delete coresqr;
        delete ex;
        delete ey;
        
        if (no_errors)
        {
            // under success will return here            
            //out_object = Py_BuildValue("OO", py_out_u, py_out_v);
            //	out_object = Py_BuildValue("i",0);
            Py_RETURN_NONE;
        }
        else
        {
            // py->c++ conversion failed must return null and trow an error?
            PyErr_SetString(PyExc_TypeError, "Output arrays were not of numTargets in length!");
            out_object = NULL;
        }
    }
    else
    {
        // py->c++ conversion failed must return null and trow an error?
        PyErr_SetString(PyExc_TypeError, "Input arrays were not of same or sufficient lengths!");
        out_object = NULL;
    }
    return out_object;
}



static PyObject * vortexElement(PyObject*, PyObject *args)
{
    PyObject * py_xl, * py_yl, * py_strl, * py_ex, * py_ey;
    PyObject * py_xr, * py_yr, * py_strr;
    PyObject * py_out_u, * py_out_v, * out_object;
    double precision, panelTolerance, assumePointLength;
    int maxThreads, targetsPerBox;
    bool isInfiniteSheet;
    // All inputs are in the tuple, need to get them out.
    bool no_errors = true;
    if (!PyArg_ParseTuple(args, "OOOOOOOOOOididdb", &py_xl, &py_yl, &py_strl, &py_xr, &py_yr, &py_strr, &py_ex, &py_ey, &py_out_u, &py_out_v, &maxThreads, &precision, &targetsPerBox, &panelTolerance, &assumePointLength, &isInfiniteSheet)) { return NULL; };
    // as we know that input are lists or arrays need to get the data out;
    //is a tuple containing many tupples.
    // create temp c-array for input to function
    double * xl, * yl, * strl, * xr, * yr, * strr, * ex, * ey, * out_u, * out_v;
    int numTargets, numParticles;

    // pass pyobject and array reference in and perform conversion.
    // setting up output arrays first as this is easy.
    numParticles = PyObject_Length(py_xl);
    numTargets = PyObject_Length(py_ex);
    if (numParticles < 1 || numTargets < 1) {no_errors = false;}

    // now targets
    // first checking they are all the same length
    if (PyObject_Length(py_ey) != numTargets) {no_errors = false;}
    if (PyObject_Length(py_yl) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_strl) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_xr) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_yr) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_strr) != numParticles) {no_errors = false;}

    if (no_errors)
    {
        // for particles first
        xl = py2cfmm(py_xl);
        yl = py2cfmm(py_yl);
        strl = py2cfmm(py_strl);
        xr = py2cfmm(py_xr);
        yr = py2cfmm(py_yr);
        strr = py2cfmm(py_strr);
        // for targets now
        ex = py2cfmm(py_ex);
        ey = py2cfmm(py_ey);
        // for output fields
        out_u = new double[numTargets];
        out_v = new double[numTargets];
        // Do the actual computations using my c function
        // Cfmm requires inputs as a continuous c-array block of memory...thats
        //why all this conversion was necessary its format is:
        //static void run( double accuracy, size_t targetsPerBox, size_t
        //maxThreads, double * x, double * y,double * str, double * core,
        //double * evalx, double * evaly,  double * evalu, double * evalv,
        //size_t numparticles, size_t numtargets);

        double twoPi = 2.0 * M_PI;

        // Need to preprocess with /2pi
        for (int index = 0; index < numParticles; index++)
        {
            strl[index] /= twoPi;
            strr[index] /= twoPi;
        }
        bool withFMM = true;
        if (precision <= 0.0)
        {
            withFMM = false;
        }

        VortexElementFMM fmm;
        fmm.run(precision, targetsPerBox, maxThreads, xl, yl, strl, xr, yr, strr, ex, ey, out_u, out_v, numParticles , numTargets, panelTolerance, assumePointLength, withFMM, isInfiniteSheet);

        // Retruning calculated values back to python now
        if ((cfmm2py(out_u, py_out_u, numTargets))) { no_errors = false; }
        if ((cfmm2py(out_v, py_out_v, numTargets))) { no_errors = false; }
         
        // clean up mess
        // void PyBuffer_Release(PyObject *obj, Py_buffer *view)¶
        delete out_u;
        delete out_v;
        delete xl;
        delete yl;
        delete strl;
        delete xr;
        delete yr;
        delete strr;
        delete ex;
        delete ey;

        if (no_errors)
        {
            // under success will return here            
            //out_object = Py_BuildValue("OO", py_out_u, py_out_v);
            //	out_object = Py_BuildValue("i",0);
            Py_RETURN_NONE;
        }
        else
        {
            // py->c++ conversion failed must return null and trow an error?
            PyErr_SetString(PyExc_TypeError, "Output arrays were not of numTargets in length!");
            out_object = NULL;
        }
    }
    else
    {
        // py->c++ conversion failed must return null and trow an error?
        PyErr_SetString(PyExc_TypeError, "Input arrays were not of same or sufficient lengths!");
        out_object = NULL;
    }
    return out_object;
}


static PyObject * lambVortexVort(PyObject*, PyObject *args)
{
    PyObject * py_x, * py_y, * py_str, * py_coresqr, * py_ex, * py_ey;
    PyObject * py_out_u, * py_out_v, * out_object;
    double precision;
    int maxThreads, targetsPerBox;
    // All inputs are in the tuple, need to get them out.
    bool no_errors = true;
    if (!PyArg_ParseTuple(args, "OOOOOOOOidi", &py_x, &py_y, &py_str, &py_coresqr,
                &py_ex, &py_ey, &py_out_u, &py_out_v, &maxThreads, &precision,
                &targetsPerBox)) { return NULL; };
    // as we know that input are lists or arrays need to get the data out;
    //is a tuple containing many tupples.
    // create temp c-array for input to function
    double * x, * y, * str, * coresqr, * ex, * ey, * out_u, * out_v;
    double * str_2PI, * core;
    int numTargets, numParticles;

    // pass pyobject and array reference in and perform conversion.
    // setting up output arrays first as this is easy.
    numParticles = PyObject_Length(py_x);
    numTargets = PyObject_Length(py_ex);
    if (numParticles < 1 || numTargets < 1) {no_errors = false;}

    // now targets
    // first checking they are all the same length
    if (PyObject_Length(py_ey) != numTargets) {no_errors = false;}
    if (PyObject_Length(py_y) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_str) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_coresqr) != numParticles) {no_errors = false;}

    if (no_errors)
    {
        // for particles first
        x = py2cfmm(py_x);
        y = py2cfmm(py_y);
        str = py2cfmm(py_str);
        core = py2cfmm(py_coresqr);
        // for targets now
        ex = py2cfmm(py_ex);
        ey = py2cfmm(py_ey);
        // for output fields
        out_u = new double[numTargets];
        out_v = new double[numTargets];
        // Do the actual computations using my c function
        // Cfmm requires inputs as a continuous c-array block of memory...thats
        //why all this conversion was necessary its format is:
        //static void run( double accuracy, size_t targetsPerBox, size_t
        //maxThreads, double * x, double * y,double * str, double * core,
        //double * evalx, double * evaly,  double * evalu, double * evalv,
        //size_t numparticles, size_t numtargets);
        double onePi = 1.0 * M_PI;

        // Need to preprocess with coresquared and /2pi
        for (int index = 0; index < numParticles; index++)
        {
            core[index] *= core[index];
            str[index] /= onePi;
        }
        coresqr = core;
        str_2PI = str;

        //timeval start,stop; gettimeofday(&start, NULL);
            LambVortexVortFMM fmm;
            fmm.run( precision, targetsPerBox, maxThreads, x, y, str_2PI, coresqr,
                    ex, ey, out_u, out_v, numParticles , numTargets);
            //gettimeofday(&stop, NULL); double timeFull = timeDiff(start, stop);
            //cout << "full time taken " << timeFull << endl;

            // Retruning calculated values back to python now
        if ((cfmm2py(out_u, py_out_u, numTargets))) { no_errors = false; }
        if ((cfmm2py(out_v, py_out_v, numTargets))) { no_errors = false; }

        // clean up mess
        // void PyBuffer_Release(PyObject *obj, Py_buffer *view)¶
        delete out_u;
        delete out_v;
        delete x;
        delete y;
        delete str;
        delete coresqr;
        delete ex;
        delete ey;

        if (no_errors)
        {
            // under success will return here
            //out_object = Py_BuildValue("OO", py_out_u, py_out_v);
            //	out_object = Py_BuildValue("i",0);
            Py_RETURN_NONE;
        }
        else
        {
            // py->c++ conversion failed must return null and trow an error?
            PyErr_SetString(PyExc_TypeError, "Output arrays were not of numTargets in length!");
            out_object = NULL;
        }
    }
    else
    {
        // py->c++ conversion failed must return null and trow an error?
        PyErr_SetString(PyExc_TypeError, "Input arrays were not of same or sufficient lengths!");
        out_object = NULL;
    }
    return out_object;
}

static PyObject * lambVortexMerge(PyObject*, PyObject *args)
{
    PyObject * py_x, * py_y, * py_str, * py_core, * py_prev_u, * py_prev_v;
    PyObject * py_out_x, * py_out_y, * py_out_str, * py_out_core, * py_out_prev_u, * py_out_prev_v, * out_object;
    int maxThreads, targetsPerBox;
    double mincoreratio, maxcoreratio, radiusratio, maxcoresizeratio;
    bool usingFMM;

    bool no_errors = true;
    if (!PyArg_ParseTuple(args, "ddddbiiOOOOOOOOOOOO", &mincoreratio, &maxcoreratio, &radiusratio, &maxcoresizeratio, &usingFMM, &targetsPerBox, &maxThreads, 
        &py_x, &py_y, &py_str, &py_core, &py_prev_u, &py_prev_v,
        &py_out_x, &py_out_y, &py_out_str, &py_out_core, &py_out_prev_u, &py_out_prev_v)) { return NULL; };
    // as we know that input are lists or arrays need to get the data out;
    //is a tuple containing many tupples.
    // create temp c-array for input to function
    double * x, * y, * str, * core, * prev_u, * prev_v,  * out_x, * out_y, * out_str, * out_core, * out_prev_u, * out_prev_v;
    int numParticles;
    size_t out_numParticles;

    // pass pyobject and array reference in and perform conversion.
    // setting up output arrays first as this is easy.
    numParticles = PyObject_Length(py_x);
    if (numParticles < 1) {no_errors = false;}

    // now targets
    // first checking they are all the same length
    if (PyObject_Length(py_y) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_str) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_core) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_prev_u) != numParticles) {no_errors = false;}
    if (PyObject_Length(py_prev_v) != numParticles) {no_errors = false;}

    if (no_errors)
    {
        // for particles first
        x = py2cfmm(py_x);
        y = py2cfmm(py_y);
        str = py2cfmm(py_str);
        core = py2cfmm(py_core);
        prev_u = py2cfmm(py_prev_u);
        prev_v = py2cfmm(py_prev_v);

        out_x = NULL;
        out_y = NULL;
        out_str = NULL;
        out_core = NULL;
        out_prev_u = NULL;
        out_prev_v = NULL;
        
        // Do the actual computations using my c function
        // Cfmm requires inputs as a continuous c-array block of memory...thats
        //why all this conversion was necessary its format is:

        //timeval start,stop; gettimeofday(&start, NULL);
        LambVortexMergeFMM fmm;
        fmm.run(mincoreratio, maxcoreratio, radiusratio, maxcoresizeratio, usingFMM, targetsPerBox, maxThreads, x, y,str, core,prev_u, prev_v,  numParticles, out_x, out_y, out_str, out_core, out_prev_u, out_prev_v, out_numParticles);
            //gettimeofday(&stop, NULL); double timeFull = timeDiff(start, stop);
            //cout << "full time taken " << timeFull << endl;

            // Retruning calculated values back to python now
        if ((cfmm2pylist(out_x, py_out_x, out_numParticles))) { no_errors = false; }
        if ((cfmm2pylist(out_y, py_out_y, out_numParticles))) { no_errors = false; }
        if ((cfmm2pylist(out_str, py_out_str, out_numParticles))) { no_errors = false; }
        if ((cfmm2pylist(out_core, py_out_core, out_numParticles))) { no_errors = false; }
        if ((cfmm2pylist(out_prev_u, py_out_prev_u, out_numParticles))) { no_errors = false; }
        if ((cfmm2pylist(out_prev_v, py_out_prev_v, out_numParticles))) { no_errors = false; }
        
        // clean up mess
        // void PyBuffer_Release(PyObject *obj, Py_buffer *view)¶
        delete out_x;
        delete out_y;
        delete out_str;
        delete out_core;
        delete out_prev_u;
        delete out_prev_v;
        delete x;
        delete y;
        delete str;
        delete core;
        delete prev_u;
        delete prev_v;

        if (no_errors)
        {
            // under success will return here
            //out_object = Py_BuildValue("OO", py_out_u, py_out_v);
            //	out_object = Py_BuildValue("i",0);
            Py_RETURN_NONE;
        }
        else
        {
            // py->c++ conversion failed must return null and trow an error?
            PyErr_SetString(PyExc_TypeError, "Output arrays were not of numTargets in length!");
            out_object = NULL;
        }
    }
    else
    {
        // py->c++ conversion failed must return null and trow an error?
        PyErr_SetString(PyExc_TypeError, "Input arrays were not of same or sufficient lengths!");
        out_object = NULL;
    }
    return out_object;
}

static PyMethodDef interfacedMethods[] = {
    {"lamb_vortex",  lambVortex, METH_VARARGS,
        "Use FMM to Calculate influence of lamb vortex on point inspace.\n\n\
            More details to follow..."},
    {"vortex_element",  vortexElement, METH_VARARGS,
        "Use FMM to Calculate influence of finite vortex sheets on point in space.\n\n\
            More details to follow..."},
    {"lamb_vortex_vort",  lambVortexVort, METH_VARARGS,
        "Use FMM to Calculate vorticity influence of lamb vortex on point inspace.\n\n\
            More details to follow..."},
    {"lamb_vortex_merge",  lambVortexMerge, METH_VARARGS,
        "Use FMM to merge vortex elements in a vormat.\n\n\
            Format for inputs are in order (\"ddddbiiOOOOOOOO\"): &mincoreratio, &maxcoreratio, &radiusratio, &maxcoresizeratio, &usingFMM, &targetsPerBox, &maxThreads,\
        &py_x, &py_y, &py_str, &py_core,\
        &py_out_x, &py_out_y, &py_out_str, &py_out_core"},
        {NULL, NULL, 0, NULL}        /* Sentinel */
};

/*// for python 3.1
static struct PyModuleDef interfacedModule = {
   PyModuleDef_HEAD_INIT,
   "velocitymodule",   // name of module
   NULL, // module documentation, may be NULL
   -1,       // size of per-interpreter state of the module,
             // or -1 if the module keeps state in global variables.
   interfacedMethods
};

PyMODINIT_FUNC PyInit_velocitymodule(void)
{
    return PyModule_Create(&interfacedModule);
}*/

// for python 2.6
PyMODINIT_FUNC initvelocitymodule(void)
{
    (void) Py_InitModule("velocitymodule", interfacedMethods);
    return;
}

double * py2cfmm(PyObject * inobj)
{
    double * convertedmem;
    int bfcheck = PyObject_CheckBuffer(inobj);
    if (bfcheck)
    {
        //careful this is python2.6 specific!
        cout << "Using buffer (2.6) interface" << endl;
        Py_buffer buffer_view;
        // get the buffer view of py_x in c-contiguous format
        int success = PyObject_GetBuffer(inobj, &buffer_view, PyBUF_C_CONTIGUOUS);
        cout << "buffer success" << success;
        convertedmem = (double *) buffer_view.buf;
    }
    else
    {
        int numTargets = PyObject_Length(inobj);
        convertedmem = new double[numTargets];
        PyObject * item;
        for (int index = 0; index < numTargets; index++)
        { /* get the element from the list/tuple */
            item = PySequence_GetItem(inobj, index);
            convertedmem[index] = PyFloat_AsDouble(item);
            Py_DECREF(item); // dont know if this is correct or going to cause a crash later
        }
    }
    return convertedmem;
}
int cfmm2py(double * contmem, PyObject * outsequence, int numTargets)
{
    // PyObject * py_converted = PyList_New(0);4
    // for (int index = 0; index < numTargets; index++)
    // {
    //     PyList_Append(py_converted, PyFloat_FromDouble(contmem [index]));
    // }
    int errors=0;
    int index = 0;
    PyObject * thedata;
    while (!errors && index < numTargets)
    {
        thedata = PyFloat_FromDouble(contmem [index]);
        errors = PySequence_SetItem(outsequence, index, thedata);
        Py_DECREF(thedata); // sequence setitem refuses to steal the reference, it creates a new one memleak used to be here!!!!!!
        index++;
    }
    return errors;
}
int cfmm2pylist(double * contmem, PyObject * py_outsequence, int numTargets)
{
    int errors=0;
    int index = 0;
    while (!errors && index < numTargets)
    {
        errors = PyList_Append(py_outsequence, PyFloat_FromDouble(contmem [index]));
        index++;
    }
    return errors;
}

/*
PyObject * cfmm2py(double * contmem, int numTargets)
{
    PyObject * py_converted = PyList_New(0);//4
    for (int index = 0; index < numTargets; index++)
    {
        PyList_Append(py_converted, PyFloat_FromDouble(contmem [index]));
    }
    return py_converted
}
*/
