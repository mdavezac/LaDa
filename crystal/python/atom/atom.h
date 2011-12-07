#include <Python.h>
#include "../../atom.h"

//! Structure holding shared pointer to an atom.
extern "C" struct Sequence
{ 
  PyObject_HEAD
  boost::shared_ptr< LaDa::crystal::AtomData< std::vector<std::string> > > ptr_atom;
  std::vector< std::string > *ptr_seq;
};

//! Structure holding shared pointer to an atom.
extern "C" struct AtomStr
{
  PyObject_HEAD
  PyObject* weakreflist;
  boost::shared_ptr< LaDa::crystal::AtomData< std::string > > atom;
};
//! Structure holding shared pointer to an atom.
extern "C" struct AtomSequence
{
  PyObject_HEAD
  Sequence *sequence;
  PyObject* weakreflist;
  boost::shared_ptr< LaDa::crystal::AtomData< std::vector<std::string> > > atom;
};

//! Type instance for string atoms.
extern "C" PyTypeObject atomstr_type;
//! Type instance for sequence atoms.
extern "C" PyTypeObject atomsequence_type;
//! Type instance for sequences.
extern "C" PyTypeObject sequence_type;
//! Returns true if an object is a string atom or subtype.
#define PyAtomStr_Check(object) PyObject_TypeCheck(object, (PyTypeObject*)&atomstr_type)
//! Returns true if an object is a string atom.
#define PyAtomStr_CheckExact(object) object->ob_type == (PyTypeObject*)&atomstr_type
//! Returns true if an object is a sequence atom or subtype.
#define PyAtomSequence_Check(object) PyObject_TypeCheck(object, (PyTypeObject*)&atomsequence_type)
//! Returns true if an object is a sequence atom.
#define PyAtomSequence_CheckExact(object) object->ob_type == (PyTypeObject*)&atomsequence_type

//! \brief Extracts a string atom from a python object.
//! \details If not a valid atom object, an exception is issued and NULL is returned.
//!          The life of the pointer should not exceed that of _atom. 
extern "C" inline LaDa::crystal::Atom<std::string> PyAtomStr_AsAtom(PyObject *_atom)
{
  if(PyAtomStr_Check(_atom)) return LaDa::crystal::Atom<std::string>(((AtomStr*)_atom)->atom);
  LADA_PYERROR(TypeError, "Object is not an atom.");
  return LaDa::crystal::Atom<std::string>();
}
//! \brief Extracts a sequence atom from a python object.
//! \details If not a valid atom object, an exception is issued and NULL is returned.
//!          The life of the pointer should not exceed that of _atom. 
extern "C" inline LaDa::crystal::Atom< std::vector<std::string> >
  PyAtomSequence_AsAtom(PyObject *_atom)
  {
    if(PyAtomSequence_Check(_atom))
      return LaDa::crystal::Atom< std::vector<std::string> >(((AtomSequence*)_atom)->atom);
    LADA_PYERROR(TypeError, "Object is not an atom.");
    return LaDa::crystal::Atom< std::vector<std::string> >();
  }
//! Extracts a string atom from a python object. No type check.
#define PyAtomStr_AS_ATOM(object) LaDa::crystal::Atom<std::string>(((AtomStr*)_atom)->atom);
//! Extracts a sequence atom from a python object.
#define PyAtomSequence_AS_ATOM(object) &((AtomSequence*)_atom)->atom;

//! \brief Creates an atom wrapper instance.
//! \brief The attribute dictionary is created on call only. It may not yet
//!         exist. g++ warns against non-template friend functions of template
//!         classes. Hence using templated form here although it doesn't make
//!         much sense.
template<class T> PyObject* PyAtom_FromAtom(LaDa::crystal::Atom<T> const &_atom);

//! Creates an atom wrapper instance around a pre-existing atom.
template<> PyObject* PyAtom_FromAtom<std::string>(LaDa::crystal::Atom<std::string> const &_atom)
{
  if(not _atom.atom_) 
  {
    LADA_PYERROR(internal, "Invalid atom.\n");
    return NULL;
  }
  AtomStr* result = (AtomStr*)atomsequence_type.tp_alloc(&atomstr_type, 0);
  if(result == NULL) return NULL;
  
  // set everything to null, just in case we exit to fast.
  result->weakreflist = NULL;
  // Now starts setting things up.
  typedef LaDa::crystal::AtomData<std::string> t_Atom;
  // According to boost, should never throw.
  new(&result->atom) boost::shared_ptr<t_Atom>(_atom.atom_); 
  
  return (PyObject*) result;
}
//! Creates an atom wrapper instance around a pre-existing atom.
template<> PyObject* PyAtom_FromAtom< std::vector<std::string> >
  (LaDa::crystal::Atom< std::vector<std::string> > const &_atom)
{
  if(not _atom.atom_) 
  {
    LADA_PYERROR(internal, "Invalid atom.\n");
    return NULL;
  }
  AtomSequence* result = (AtomSequence*)atomsequence_type.tp_alloc(&atomsequence_type, 0);
  if(result == NULL) return NULL;
  
  // set everything to null, just in case we exit to fast.
  result->weakreflist = NULL;
  result->sequence = NULL;
  // Now starts setting things up.
  typedef LaDa::crystal::AtomData< std::vector<std::string> > t_Atom;
  // According to boost, should never throw.
  new(&result->atom) boost::shared_ptr<t_Atom>(_atom.atom_); 
  
  result->sequence = (Sequence*)PyObject_New(Sequence, &sequence_type);
  if(result->sequence == NULL)
  {
    Py_DECREF(result);
    return NULL;
  }
  result->sequence->ptr_seq = &result->atom->type;
  typedef LaDa::crystal::AtomData< std::vector<std::string> > t_Atom;
  // something weird here. Seems that the constructor for the boost shared pointer was never called.
  // Perfome in-place construction using new.
  // According to boost, should never throw.
  new(&result->sequence->ptr_atom) boost::shared_ptr<t_Atom>(result->atom); 

  return (PyObject*) result;
}

//! Creates a new atom and its wrapper.
extern "C" PyObject* PyAtomStr_New();
//! Creates a new atom and its wrapper.
extern "C" PyObject* PyAtomSequence_New();
