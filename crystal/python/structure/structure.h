#ifndef LADA_CRYSTAL_PYTHON_STRUCTURE_STRUCTURE_H
#define LADA_CRYSTAL_PYTHON_STRUCTURE_STRUCTURE_H
#include <Python.h>
#include "../../structure.h"

//! Structure holding shared pointer to a structure.
extern "C" struct StructureStr
{
  PyObject_HEAD
  PyObject* weakreflist;
  boost::shared_ptr< LaDa::crystal::StructureData< std::string > > structure;
};
//! Structure holding shared pointer to a structure.
extern "C" struct StructureSequence
{
  PyObject_HEAD
  PyObject* weakreflist;
  boost::shared_ptr< LaDa::crystal::StructureData< std::vector<std::string> > > structure;
};

//! Type instance for string structures.
extern "C" PyTypeObject structurestr_type;
//! Type instance for sequence structures.
extern "C" PyTypeObject structuresequence_type;
//! Returns true if an object is a string structure or subtype.
#define PyStructureStr_Check(object) PyObject_TypeCheck(object, (PyTypeObject*)&structurestr_type)
//! Returns true if an object is a string structure.
#define PyStructureStr_CheckExact(object) object->ob_type == (PyTypeObject*)&structurestr_type
//! Returns true if an object is a sequence structure or subtype.
#define PyStructureSequence_Check(object) PyObject_TypeCheck(object, (PyTypeObject*)&structuresequence_type)
//! Returns true if an object is a sequence structure.
#define PyStructureSequence_CheckExact(object) object->ob_type == (PyTypeObject*)&structuresequence_type

//! \brief Extracts a string structure from a python object.
//! \details If not a valid structure object, an exception is issued and and
//!          empty Structure is returned.
extern "C" inline LaDa::crystal::Structure<std::string> PyStructureStr_AsStructure(PyObject *_structure)
{
  if(PyStructureStr_Check(_structure))
    return LaDa::crystal::Structure<std::string>(((StructureStr*)_structure)->structure);
  LADA_PYERROR(TypeError, "Object is not a structure.");
  return LaDa::crystal::Structure<std::string>();
}
//! \brief Extracts a sequence structure from a python object.
//! \details If not a valid structure object, an exception is issued and and
//!          empty Structure is returned.
extern "C" inline LaDa::crystal::Structure< std::vector<std::string> >
  PyStructureSequence_AsStructure(PyObject *_structure)
  {
    typedef LaDa::crystal::Structure< std::vector<std::string> > t_Structure;
    if(PyStructureSequence_Check(_structure))
      return t_Structure(((StructureSequence*)_structure)->structure);
    LADA_PYERROR(TypeError, "Object is not an structure.");
    return t_Structure();
  }
//! Extracts a string structure from a python object. No type check.
#define PyStructureStr_AS_STRUCTURE(object)\
  LaDa::crystal::Structure<std::string>(((StructureStr*)object)->structure);
//! Extracts a sequence structure from a python object.
#define PyStructureSequence_AS_STRUCTURE(object) &((StructureSequence*)object)->structure;

//! \brief Creates an structure wrapper instance.
//! \brief The attribute dictionary is created on call only. It may not yet
//!         exist. g++ warns against non-template friend functions of template
//!         classes. Hence using templated form here although it doesn't make
//!         much sense.
template<class T> PyObject* PyStructure_FromStructure(LaDa::crystal::Structure<T> const &_structure);
//! \brief Creates an structure wrapper instance.
//! \brief The attribute dictionary is created on call only. It may not yet
//!         exist. g++ warns against non-template friend functions of template
//!         classes. Hence using templated form here although it doesn't make
//!         much sense.
template<class T> PyObject* PyStructure_FromStructure
   (boost::shared_ptr< LaDa::crystal::StructureData<T> > const &_structure);

// The following is very ad-hoc, but then I kept screwing up Sequence and Str
// in four functions below.
# ifdef LADA_DECLARE
#   error LADA_DECLARE already defined.
# endif
# define LADA_DECLARE(TYPE, name, Name, dotstructure)                                      \
    /* First check the shared pointer is not null. Should never actually happen, but... */ \
    if(not _structure dotstructure)                                                        \
    {                                                                                      \
      LADA_PYERROR(internal, "Invalid structure.\n");                                      \
      return NULL;                                                                         \
    }                                                                                      \
    /* check for existence of a wrapper first. */                                          \
    if(_structure->pyself) { Py_INCREF(_structure->pyself); return _structure->pyself; }   \
                                                                                           \
    Structure ## Name * result = (Structure ## Name*)structure ## name ## _type            \
                                  .tp_alloc(&structure ## name ## _type, 0);               \
    if(result == NULL) return NULL;                                                        \
                                                                                           \
    /* set everything to null, just in case we exit too fast. */                           \
    result->weakreflist = NULL;                                                            \
    /* Now starts setting things up. */                                                    \
    typedef LaDa::crystal::StructureData< TYPE > t_Structure;                              \
    /* According to boost, should never throw. */                                          \
    new(&result->structure) boost::shared_ptr<t_Structure>(_structure dotstructure);       \
                                                                                           \
    /* Now sets internal reference, but no incref'ing. */                                  \
    /* Note that using dotstructure (.impl_) allows access to an otherwise                 \
     * constant variable. The beauty of using pointers. */                                 \
    _structure dotstructure->pyself = (PyObject*)result;                                   \
    return (PyObject*) result;

//! Creates a structure wrapper instance around a pre-existing structure.
template<> PyObject* PyStructure_FromStructure<std::string>
  (LaDa::crystal::Structure<std::string> const &_structure)
{
  LADA_DECLARE(std::string, str, Str, .impl_);
}
//! Creates a structure wrapper instance around a pre-existing structure.
template<> PyObject* PyStructure_FromStructure<std::string>
  (boost::shared_ptr< LaDa::crystal::StructureData<std::string> > const &_structure)
{
  LADA_DECLARE(std::string, str, Str, );
}
//! Creates a structure wrapper instance around a pre-existing structure.
template<> PyObject* PyStructure_FromStructure< std::vector<std::string> >
  (LaDa::crystal::Structure< std::vector<std::string> > const &_structure)
{
  LADA_DECLARE( std::vector<std::string>, sequence, Sequence, .impl_);
}
//! Creates a structure wrapper instance around a pre-existing structure.
template<> PyObject* PyStructure_FromStructure< std::vector<std::string> >
  (boost::shared_ptr< LaDa::crystal::StructureData< std::vector<std::string> > > const &_structure)
{
  LADA_DECLARE( std::vector<std::string>, sequence, Sequence,);
}
# undef LADA_DECLARE



//! Creates a new structure and its wrapper.
extern "C" PyObject* PyStructureStr_New();
//! Creates a new structure and its wrapper.
extern "C" PyObject* PyStructureSequence_New();
#endif
