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
//! \details If not a valid atom object, an exception is issued and a default
//!          atom is returned.
extern "C" inline LaDa::crystal::Atom<std::string> PyAtomStr_AsAtom(PyObject *_atom)
{
  if(PyAtomStr_Check(_atom)) return LaDa::crystal::Atom<std::string>(((AtomStr*)_atom)->atom);
  LADA_PYERROR(TypeError, "Object is not an atom.");
  return LaDa::crystal::Atom<std::string>();
}
//! \brief Extracts a sequence atom from a python object.
//! \details If not a valid atom object, an exception is issued and a default
//!          atom is returned.
extern "C" inline LaDa::crystal::Atom< std::vector<std::string> >
  PyAtomSequence_AsAtom(PyObject *_atom)
  {
    if(PyAtomSequence_Check(_atom))
      return LaDa::crystal::Atom< std::vector<std::string> >(((AtomSequence*)_atom)->atom);
    LADA_PYERROR(TypeError, "Object is not an atom.");
    return LaDa::crystal::Atom< std::vector<std::string> >();
  }
//! Extracts a string atom from a python object. No type check.
#define PyAtomStr_AS_ATOM(object) LaDa::crystal::Atom<std::string>(((AtomStr*)object)->atom);
//! Extracts a sequence atom from a python object.
#define PyAtomSequence_AS_ATOM(object) &((AtomSequence*)object)->atom;

//! \brief Creates an atom wrapper instance.
//! \brief The attribute dictionary is created on call only. It may not yet
//!         exist. g++ warns against non-template friend functions of template
//!         classes. Hence using templated form here although it doesn't make
//!         much sense.
template<class T> PyObject* PyAtom_FromAtom(LaDa::crystal::Atom<T> const &_atom);
//! \brief Creates an atom wrapper instance.
//! \brief The attribute dictionary is created on call only. It may not yet
//!         exist. g++ warns against non-template friend functions of template
//!         classes. Hence using templated form here although it doesn't make
//!         much sense.
template<class T> PyObject* PyAtom_FromAtom(boost::shared_ptr< LaDa::crystal::AtomData<T> > const &_atom);

// The following is very ad-hoc, but then I kept screwing up..
# ifdef LADA_DECLARE
#   error LADA_DECLARE already defined.
# endif
# ifdef LADA_SEQUENCE_NULL
#   error LADA_SEQUENCE_NULL already defined.
# endif
# ifdef LADA_SEQUENCE_SET
#   error LADA_SEQUENCE_SET already defined.
# endif
# define LADA_DECLARE(TYPE, name, Name, alltype, dotatom)                                  \
  template<> PyObject* PyAtom_FromAtom< TYPE >(alltype const &_atom)                       \
  {                                                                                        \
    if(not _atom dotatom)                                                                  \
    {                                                                                      \
      LADA_PYERROR(internal, "Invalid atom.\n");                                           \
      return NULL;                                                                         \
    }                                                                                      \
    /* check for existence of a wrapper first. */                                          \
    if(_atom->pyself) { Py_INCREF(_atom->pyself); return _atom->pyself; }                  \
                                                                                           \
    Atom ## Name * result = (Atom ## Name*)atom ## name ## _type                           \
                                  .tp_alloc(&atom ## name ## _type, 0);                    \
    if(result == NULL) return NULL;                                                        \
                                                                                           \
    /* set everything to null, just in case we exit to fast. */                            \
    result->weakreflist = NULL;                                                            \
    /* Now starts setting things up. */                                                    \
    typedef LaDa::crystal::AtomData< TYPE > t_Atom;                                        \
    LADA_SEQUENCE_NULL;                                                                    \
    /* According to boost, should never throw. */                                          \
    new(&result->atom) boost::shared_ptr<t_Atom>(_atom dotatom);                           \
                                                                                           \
    LADA_SEQUENCE_SET;                                                                     \
                                                                                           \
    /* Now sets internal reference, but no incref'ing. */                                  \
    _atom dotatom->pyself = (PyObject*)result;                                             \
    return (PyObject*) result;                                                             \
  }

# define LADA_SEQUENCE_NULL
# define LADA_SEQUENCE_SET
//! Creates an atom wrapper instance around a pre-existing atom.
LADA_DECLARE(std::string, str, Str, LaDa::crystal::Atom<std::string>, .atom_);
//! Creates an atom wrapper instance around a pre-existing atom.
LADA_DECLARE(std::string, str, Str, boost::shared_ptr< LaDa::crystal::AtomData<std::string> >, );
# undef LADA_SEQUENCE_NULL
# undef LADA_SEQUENCE_SET

# define LADA_SEQUENCE_NULL result->sequence = NULL;
# define LADA_SEQUENCE_SET                                                    \
    result->sequence = (Sequence*)PyObject_New(Sequence, &sequence_type);     \
    if(result->sequence == NULL)                                              \
    {                                                                         \
      Py_DECREF(result);                                                      \
      return NULL;                                                            \
    }                                                                         \
    result->sequence->ptr_seq = &result->atom->type;                          \
    /* something weird here. Seems that the constructor for the boost shared  \
       pointer was never called.  Perfome in-place construction using new.    \
       According to boost, should never throw. */                             \
    new(&result->sequence->ptr_atom) boost::shared_ptr<t_Atom>(result->atom); 
//! Creates an atom wrapper instance around a pre-existing atom.
LADA_DECLARE(std::vector<std::string>, sequence, Sequence, 
             LaDa::crystal::Atom< std::vector<std::string> >, .atom_);
//! Creates an atom wrapper instance around a pre-existing atom.
LADA_DECLARE(std::vector<std::string>, sequence, Sequence, 
             boost::shared_ptr< LaDa::crystal::AtomData< std::vector<std::string> > >, )
# undef LADA_DECLARE
# undef LADA_SEQUENCE_NULL
# undef LADA_SEQUENCE_SET



//! Creates a new atom and its wrapper.
extern "C" PyObject* PyAtomStr_New();
//! Creates a new atom and its wrapper.
extern "C" PyObject* PyAtomSequence_New();
