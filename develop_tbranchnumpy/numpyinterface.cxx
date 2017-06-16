#include <deque>
#include <string>
#include <vector>

#include <Python.h>
#include <numpy/arrayobject.h>

#include <ROOT/TBulkBranchRead.hxx>
#include <TBranch.h>
#include <TBufferFile.h>
#include <TBuffer.h>
#include <TClass.h>
#include <TDataType.h>
#include <TFile.h>
#include <TLeafB.h>
#include <TLeafB.h>
#include <TLeafD.h>
#include <TLeafF.h>
#include <TLeaf.h>
#include <TLeafI.h>
#include <TLeafI.h>
#include <TLeafL.h>
#include <TLeafL.h>
#include <TLeafO.h>
#include <TLeafS.h>
#include <TLeafS.h>
#include <TObjArray.h>
#include <TTree.h>

/////////////////////////////////////////////////////// helper classes

// performance counters for diagnostics
static Long64_t baskets_loaded = 0;
static Long64_t bytes_loaded = 0;
static Long64_t baskets_copied = 0;
static Long64_t bytes_copied = 0;
static Long64_t items_scanned = 0;
static Long64_t items_copied = 0;

class BasketBuffer {
public:
  Long64_t entry_start;
  Long64_t entry_end;
  TBufferFile buffer;

  BasketBuffer() : entry_start(0), entry_end(0), buffer(TBuffer::kWrite, 32*1024) {}

  void read_basket(Long64_t entry, TBranch* branch) {
    entry_start = entry;
    entry_end = entry_start + branch->GetBulkRead().GetEntriesSerialized(entry, buffer);
    if (entry_end < entry_ start)
      entry_end = -1;
  }
};

class BranchData {
public:
  TBranch* branch;
  std::deque<BasketBuffer*> buffers;
  std::vector<char> extra_buffer;
  BranchData* counter;

  BranchData(TBranch* branch) : branch(branch) {
    buffers.push_back(new BasketBuffer);
  }

  ~BranchData() {
    while (!buffers.empty()) {
      delete buffers.front();
      buffers.pop_front();
    }
  }

  void* getdata(Long64_t &numbytes, Long64_t itemsize, Long64_t entry_start, Long64_t entry_end, Long64_t alignment) {
    if (counter == NULL) {
      Long64_t fill_end = -1;

      for (unsigned int i = 0;  i < buffers.size();  ++i) {
        BasketBuffer* buf = buffers[i];

        if (entry_start == buf->entry_start  &&  entry_end == buf->entry_end  &&  (alignment <= 0  ||  (size_t)buf->buffer.GetCurrent() % alignment == 0)) {
          // this whole buffer is exactly right, in terms of start/end and alignment; don't mess with extra_buffer, just send it (no copy)!
          numbytes = buf->buffer.BufferSize();
          return buf->buffer.GetCurrent();
        }
        else if (buf->entry_start <= entry_start  &&  entry_start < buf->entry_end) {
          if (entry_end <= buf->entry_end)
            fill_end = entry_end;
          else
            fill_end = buf->entry_end;

          // where *within this buffer* should we start and end the slice?
          Long64_t byte_start = (entry_start - buf->entry_start) * itemsize;
          Long64_t byte_end = (fill_end - buf->entry_start) * itemsize;

          // this is the first buffer in which we see the start, so we *replace* extra_buffer
          extra_buffer.resize(byte_end - byte_start);
          memcpy(extra_buffer.data(), &buf->buffer.GetCurrent()[byte_start], byte_end - byte_start);
        }

        else if (entry_start < buf->entry_start  &&  buf->entry_start < entry_end) {
          if (entry_end <= buf->entry_end)
            fill_end = entry_end;
          else
            fill_end = buf->entry_end;

          Long64_t byte_end = (fill_end - buf->entry_start) * itemsize;

          // this is not the first buffer with content that we want (may or may not be last), so we *append* extra_buffer
          size_type oldsize = extra_buffer.size();
          extra_buffer.resize(oldsize + byte_end);
          memcpy(&extra_buffer.data()[oldsize], buf->buffer.GetCurrent(), byte_end);
        }
      }

      if (fill_end != entry_end)
        return NULL;
      else
        return extra_buffer.data();
    }
    else {
      // die, just die
      return NULL;
    }
  }
};

/////////////////////////////////////////////////////// Python module

class ArrayInfo {
public:
  PyArray_Descr* dtype;
  int nd;
  std::vector<int> dims;
  bool varlen;
};

typedef struct {
  PyObject_HEAD
  Long64_t alignment;
  Long64_t num_entries;
  Long64_t entry_start;
  Long64_t entry_end;
  std::vector<BranchData> requested;
  std::vector<ArrayInfo> arrayinfo;
  std::vector<BranchData> extra_counters;
} BranchesIterator;

static PyObject* BranchesIterator_iter(PyObject* self);
static PyObject* BranchesIterator_next(PyObject* self);

static PyObject* iterate(PyObject* self, PyObject* args);

#if PY_MAJOR_VERSION >=3
#define Py_TPFLAGS_HAVE_ITER 0
#endif

static PyTypeObject BranchesIteratorType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "numpyinterface.BranchesIterator", /*tp_name*/
  sizeof(BranchesIterator), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  0,                         /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER, /* tp_flags */
  "Iterator over selected TTree branches, yielding a tuple of (entry_start, entry_end, *arrays) for each cluster.", /* tp_doc */
  0,                         /* tp_traverse */
  0,                         /* tp_clear */
  0,                         /* tp_richcompare */
  0,                         /* tp_weaklistoffset */
  BranchesIterator_iter, /* tp_iter: __iter__() method */
  BranchesIterator_next  /* tp_iternext: __next__() method */
};

static PyMethodDef module_methods[] = {
  {"iterate", (PyCFunctionWithKeywords)iterate, METH_VARARGS | METH_KEYWORDS, "Get an iterator over a selected set of TTree branches, yielding a tuple of (entry_start, entry_end, *arrays) for each cluster.\n\n    filePath (str): name of the TFile\n    treePath (str): name of the TTree\n    *branchNames (strs): name of requested branches\n\nAlternatively, TBranch objects from PyROOT may be supplied (FIXME).\n\n    alignment=0: if supplied and positive, guarantee that the data are aligned to this number of bytes, even if that means copying data."},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "numpyinterface",
  NULL,
  0,
  module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};

PyMODINIT_FUNC PyInit_numpyinterface(void) {
  BranchesIteratorType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&BranchesIteratorType) < 0)
    return NULL;

  PyObject* module = PyModule_Create(&moduledef);
  if (module == NULL)
    return NULL;

  import_array();

  Py_INCREF(&BranchesIteratorType);
  PyModule_AddObject(module, "BranchesIterator", (PyObject*)&BranchesIteratorType);

  return module;
}

#else // PY_MAJOR_VERSION <= 2

PyMODINIT_FUNC initnumpyinterface(void) {
  BranchesIteratorType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&BranchesIteratorType) < 0)
    return;

  PyObject* module = Py_InitModule3("numpyinterface", module_methods, "");
  if (module == NULL)
    return;

  Py_INCREF(&BranchesIteratorType);
  PyModule_AddObject(module, "BranchesIterator", (PyObject*)&BranchesIteratorType);

  if (module != NULL)
    import_array();
}

#endif

/////////////////////////////////////////////////////// utility functions

bool getfile(TFile* &file, char* filePath) {
  file = TFile::Open(filePath);
  if (file == NULL  ||  !file->IsOpen()) {
    PyErr_Format(PyExc_IOError, "could not open file \"%s\"", filePath);
    return false;
  }
  else
    return true;
}

bool gettree(TTree* &tree, TFile* file, char* filePath, char* treePath) {
  file->GetObject(treePath, tree);
  if (tree == NULL) {
    PyErr_Format(PyExc_IOError, "could not read tree \"%s\" from file \"%s\"", treePath, filePath);
    return false;
  }
  else
    return true;
}

bool getbranch(TBranch* &branch, TTree* tree, char* filePath, char* treePath, char* branchName) {
  branch = tree->GetBranch(branchName);
  if (branch == NULL) {
    PyErr_Format(PyExc_IOError, "could not read branch \"%s\" from tree \"%s\" from file \"%s\"", branchName, treePath, filePath);
    return false;
  }
  else
    return true;
}

const char* leaftype(TLeaf* leaf) {
  if (leaf->IsA() == TLeafO::Class()) {
    return "bool";
  }
  else if (leaf->IsA() == TLeafB::Class()  &&  leaf->IsUnsigned()) {
    return "u1";
  }
  else if (leaf->IsA() == TLeafB::Class()) {
    return "i1";
  }
  else if (leaf->IsA() == TLeafS::Class()  &&  leaf->IsUnsigned()) {
    return ">u2";
  }
  else if (leaf->IsA() == TLeafS::Class()) {
    return ">i2";
  }
  else if (leaf->IsA() == TLeafI::Class()  &&  leaf->IsUnsigned()) {
    return ">u4";
  }
  else if (leaf->IsA() == TLeafI::Class()) {
    return ">i4";
  }
  else if (leaf->IsA() == TLeafL::Class()  &&  leaf->IsUnsigned()) {
    return ">u8";
  }
  else if (leaf->IsA() == TLeafL::Class()) {
    return ">i8";
  }
  else if (leaf->IsA() == TLeafF::Class()) {
    return ">f4";
  }
  else if (leaf->IsA() == TLeafD::Class()) {
    return ">f8";
  }
  else {
    TClass* expectedClass;
    EDataType expectedType;
    leaf->GetBranch()->GetExpectedType(expectedClass, expectedType);
    switch (expectedType) {
      case kBool_t:     return "bool";
      case kUChar_t:    return "u1";
      case kchar:       return "i1";
      case kChar_t:     return "i1";
      case kUShort_t:   return ">u2";
      case kShort_t:    return ">i2";
      case kUInt_t:     return ">u4";
      case kInt_t:      return ">i4";
      case kULong_t:    return ">u8";
      case kLong_t:     return ">i8";
      case kULong64_t:  return ">u8";
      case kLong64_t:   return ">i8";
      case kFloat_t:    return ">f4";
      case kDouble32_t: return ">f4";
      case kDouble_t:   return ">f8";
    }
  }
  return NULL;
}

void getdim(TLeaf* leaf, std::vector<int>& dims, std::vector<std::string>& counters) {
  const char* title = leaf->GetTitle();
  bool iscounter = false;

  for (const char* c = title;  *c != 0;  c++) {
    if (*c == '[') {
      dims.push_back(0);
      counters.push_back(std::string());
      iscounter = false;
    }

    else if (*c == ']') {
      if (iscounter)           // a dimension either fills int-valued dims or string-valued counters
        dims.pop_back();       // because we don't handle both
      else
        counters.pop_back();
    }

    else if (!dims.empty()) {
      if ('0' <= *c  &&  *c <= '9')
        dims.back() = dims.back() * 10 + (*c - '0');
      else
        iscounter = true;
      counters.back() = counters.back() + *c;   // accumulate any char that isn't '[' or ']'
    }

    // else this is part of the TLeaf name (before the first '[')
  }
}

/////////////////////////////////////////////////////// Python iterator functions

static PyObject* BranchesIterator_iter(PyObject* self) {
  Py_INCREF(self);
  return self;
}

bool update_BranchesIterator(BranchesIterator* thyself, const char* &error_string) {





}


static PyObject* BranchesIterator_next(PyObject* self) {
  BranchesIterator* thyself = reinterpret_cast<BranchesIterator*>(self);

  const char* error_string = NULL;
  bool done = update_BranchesIterator(thyself, error_string);

  if (error_string != NULL) {
    PyErr_SetString(PyExc_IOError, error_string);
    return NULL;
  }
  else if (done) {
    PyErr_SetNone(PyExc_StopIteration);
    return NULL;
  }
  else {
    PyObject* out = PyTuple_New(2 + thyself->requested.size());
    PyTuple_SET_ITEM(out, 0, PyLong_FromLong(thyself->entry_start));
    PyTuple_SET_ITEM(out, 1, PyLong_FromLong(thyself->entry_end));

    for (unsigned int i = 0;  i < thyself->requested.size();  i++) {
      Long64_t numbytes;
      void* data = thyself->requested[i].getdata(numbytes, thyself->dtypes[i].elsize, thyself->entry_start, thyself->entry_end, thyself->alignment);

      npy_intp dims[1];
      dims[0] = numbytes / thyself->dtypes[i].elsize;

      int flags = NPY_ARRAY_C_CONTIGUOUS;
      if (thyself->alignment > 0)
        flags |= NPY_ARRAY_ALIGNED;

      PyObject* array = PyArray_NewFromDescr(&PyArray_Type, thyself->dtypes[i], 1, dims, NULL, data, flags, NULL);

      PyTuple_SET_ITEM(out, i + 2, array);
    }

    return out;
  }
}

/////////////////////////////////////////////////////// iterator-making function and its helpers

bool dtypedim(PyArray_Descr* &dtype, TLeaf* leaf) {
  dtype = PyArray_DescrFromType(0);

  const char* asstring = leaftype(leaf);
  if (asstring == NULL) {
    PyErr_Format(PyExc_ValueError, "cannot convert type of TLeaf \"%s\" to Numpy", leaf->GetName());
    return false;
  }

  if (!PyArray_DescrConverter(asstring, &dtype)) {
    PyErr_SetString(PyExc_ValueError, "cannot create a dtype");
    return NULL;
  }

  return true;
}

bool dtypedim_unileaf(ArrayInfo &arrayinfo, TLeaf* leaf) {
  std::vector<std::string> counters;
  getdim(leaf, arrayinfo->dims, counters);
  arrayinfo->nd = 1 + arrayinfo->dims.size();    // first dimension is for the set of entries itself
  arrayinfo->varlen = !counters.empty();

  if (arrayinfo->nd > 1  &&  arrayinfo->varlen) {
    PyErr_Format("TLeaf \"%s\" has both fixed-length dimensions and variable-length dimensions", leaf->GetTitle());
    return false;
  }

  return dtypedim(arrayinfo->dtype, leaf);
}

bool dtypedim_multileaf(ArrayInfo &arrayinfo, TObjArray* leaves) {
  // TODO: Numpy recarray dtype
  PyErr_SetString(PyExc_NotImplementedError, "multileaf");
  return false;
}

bool dtypedim_branch(ArrayInfo &arrayinfo, TBranch* branch) {
  TObjArray* subbranches = branch->GetListOfBranches();
  if (subbranches->GetEntries() != 0) {
    PyErr_Format(PyErr_ValueError, "TBranch \"%s\" has subbranches; only branches of TLeaves are allowed", branch->GetName());
    return false;
  }

  TObjArray* leaves = branch->GetListOfLeaves();
  if (leaves->GetEntries() == 1)
    return dtypedim_unileaf(arrayinfo, dynamic_cast<TLeaf*>(leaves->First()));
  else
    return dtypedim_multileaf(arrayinfo, leaves);
}

const char* gettuplestring(PyObject* p, Py_ssize_t pos) {
  PyObject* obj = PyTuple_GET_ITEM(p, pos);
  if (PyString_Check(obj))
    return PyString_AsString(obj);
  else {
    PyErr_Format(PyExc_TypeError, "expected a string in argument %d", pos);
    return NULL;
  }
}

static PyObject* iterate(PyObject* self, PyObject* args, PyObject* kwds) {
  std::vector<TBranch*> branches;
  Long64_t alignment = 0;

  if (PyTuple_GET_SIZE(args) < 1) {
    PyErr_SetString(PyExc_TypeError, "at least one argument is required");
    return NULL;
  }

  if (PyString_Check(PyTuple_GET_ITEM(args, 0))) {
    // first argument is a string: filePath, treePath, branchNames... signature

    if (PyTuple_GET_SIZE(args) < 3) {
      PyErr_SetString(PyExc_TypeError, "in the string-based signture, at least three arguments are required");
      return NULL;
    }

    const char* filePath = gettuplestring(args, 0);
    const char* treePath = gettuplestring(args, 1);
    if (filePath == NULL  ||  treePath == NULL)
      return NULL;

    TFile* file;
    if (!getfile(file, filePath)) return NULL;

    TTree* tree;
    if (!gettree(tree, file, filePath, treePath)) return NULL;

    for (int i = 2;  i < PyTuple_GET_SIZE(args);  i++) {
      const char* branchName = gettuplestring(args, i);
      TBranch* branch;
      if (!getbranch(branch, tree, filePath, treePath, branchName)) return NULL;
      branches.push_back(branch);
    }

    if (kwds != NULL) {
      PyObject* py_alignment = PyDict_GetItemString(kwds, "alignment");
      if (py_alignment != NULL) {
        if (!PyInt_Check(py_alignment)) {
          PyErr_SetString(PyExc_TypeError, "alignment must be an integer");
          return NULL;
        }
        alignment = PyInt_AsLong(py_alignment);

        if (PyDict_Size(kwds) != 1) {
          PyErr_SetString(PyExc_TypeError, "only one keyword expected");
          return NULL;
        }
      }
      else if (PyDict_Size(kwds) != 0) {
        PyErr_SetString(PyExc_TypeError, "only one keyword expected");
        return NULL;
      }
    }
  }
  else {
    // first argument is an object: TBranch, TBranch, TBranch... signature
    // TODO: insist that all branches come from the same TTree
    PyErr_SetString(PyExc_NotImplementedError, "FIXME: accept PyROOT TBranches");
    return NULL;
  }

  BranchesIterator* out = PyObject_New(BranchesIterator, &BranchesIteratorType);

  if (!PyObject_Init((PyObject*)out, &BranchesIteratorType)) {
    Py_DECREF(out);
    return NULL;
  }

  out->alignment = alignment;
  out->num_entries = branches.front()->GetTree()->GetEntries();
  out->entry_start = 0;
  out->entry_end = 0;

  for (unsigned int i = 0;  i < branches.size();  i++) {
    out->requested.push_back(BranchData(branches[i]));
    out->arrayinfo.push_back(ArrayInfo());

    if (!dtypedim_branch(out->arrayinfo.back(), branches[i]))
      return NULL;
  }
  
  // FIXME: find all the counters, link them up, and put any missing from "requested" into "extra_counters"

  return (PyObject*)out;
}
