#include "LaDaConfig.h"

#ifdef LADA_DO_PYTHON
# include "Python.h"
#endif

#include <boost/regex.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "cast.h"


namespace LaDa
{
  namespace crystal
  {
    // \brief Casts a sequence atom to a string atom.
    // \param[in] _in: Atom to cast.
    // \param[inout] _in: resulting cast atom.
    // \param[_in] _sep: Separator between atomic type in resulting atom.
    void cast( AtomData< std::vector<std::string> > const &_in, 
               AtomData< std::string> &_out, std::string const &_sep)
    {
#     ifdef LADA_DO_PYTHON
        if(_out.pydict != NULL) 
        {
          Py_DECREF(_out.pydict);
          _out.pydict = NULL;
        }
#     endif
      _out.pos =    _in.pos;
      _out.site =   _in.site;
      _out.freeze = _in.freeze;
      if(_in.type.size() > 0)
      {
        std::vector<std::string>::const_iterator i_first = _in.type.begin();
        std::vector<std::string>::const_iterator const i_end = _in.type.end();
        std::ostringstream sstr;
        sstr << *i_first;
        for(++i_first; i_first != i_end; ++i_first)
          sstr << _sep << *i_first;
        _out.type = sstr.str();
      }
#     ifdef LADA_DO_PYTHON
        if(_in.pydict != NULL)
        {
          if(PyDict_Size(_in.pydict) > 0)
          {  
            PyObject* copymod = PyImport_ImportModule("copy");
            if(copymod == NULL) return;
            PyObject *deepcopystr = PyString_FromString("deepcopy");
            if(not deepcopystr) { Py_DECREF(copymod); return; }
            _out.pydict = PyObject_CallMethodObjArgs(copymod, deepcopystr, _in.pydict, NULL);
            Py_DECREF(copymod);
            Py_DECREF(deepcopystr);
          }
        }
#     endif
    }
    // \brief Casts a sequence atom to a string atom.
    // \param[in] _in: Atom to cast.
    // \param[_in] _sep: Separator between atomic type in resulting atom.
    Atom<std::string> cast(Atom< std::vector<std::string> > const &_in, std::string const &_sep)
    {
      boost::shared_ptr< AtomData<std::string> > result(new AtomData<std::string>);
      cast(*_in.get(), *result, _sep);
#     ifdef LADA_DO_PYTHON
        if(PyErr_Occurred() != NULL) return Atom<std::string>();
#     endif
      return Atom<std::string>(result);
    }
    // \brief Casts a string atom to a sequence atom.
    // \param[in] _in: Atom to cast.
    // \param[inout] _in: resulting cast atom.
    //! \param[_in] _sep: Separator between atomic type in resulting atom.
    void cast( AtomData<std::string> const &_in, 
               AtomData< std::vector<std::string> > &_out, std::string const &_sep)
    {
#     ifdef LADA_DO_PYTHON
        if(_out.pydict != NULL) 
        {
          Py_DECREF(_out.pydict);
          _out.pydict = NULL;
        }
#     endif
      _out.pos =    _in.pos;
      _out.site =   _in.site;
      _out.freeze = _in.freeze;
      if(_in.type.size() > 0)
      {
        _out.type.clear();
        std::string sep(_sep);
        boost::trim(sep);
        boost::regex re(sep);
        boost::sregex_token_iterator i(_in.type.begin(), _in.type.end(), re, -1);
        boost::sregex_token_iterator j;
        while(i != j)
        {
          _out.type.push_back(*i++);
          boost::trim(_out.type.back());
        }
      }
#     ifdef LADA_DO_PYTHON
        if(_in.pydict != NULL)
        {
          if(PyDict_Size(_in.pydict) > 0)
          {  
            PyObject* copymod = PyImport_ImportModule("copy");
            if(copymod == NULL) return;
            PyObject *deepcopystr = PyString_FromString("deepcopy");
            if(not deepcopystr) { Py_DECREF(copymod); return; }
            _out.pydict = PyObject_CallMethodObjArgs(copymod, deepcopystr, _in.pydict, NULL);
            Py_DECREF(copymod);
            Py_DECREF(deepcopystr);
          }
        }
#     endif
    }
    // \brief Casts a string atom to a sequence.
    // \param[in] _in: Atom to cast.
    // \param[_in] _sep: Different species will be split using this separator.
    Atom< std::vector<std::string> > cast(Atom<std::string> const &_in, std::string const &_sep)
    {
      boost::shared_ptr< AtomData< std::vector<std::string> > >
        result(new AtomData< std::vector<std::string> >);
      cast(*_in.get(), *result, _sep);
#     ifdef LADA_DO_PYTHON
        if(PyErr_Occurred() != NULL) return Atom< std::vector<std::string> >();
#     endif
      return Atom< std::vector<std::string> >(result);
    }
    // \brief Casts a string structure to a sequence structure.
    // \param[in] _in: Structure to cast.
    // \param[inout] _in: resulting cast structure.
    // \param[_in] _sep: Separator between atomic type in resulting structure.
    void cast( StructureData<std::string> const &_in, 
               StructureData< std::vector<std::string> > &_out, std::string const &_sep)
    {
#     ifdef LADA_DO_PYTHON
        if(_out.pydict != NULL) 
        {
          Py_DECREF(_out.pydict);
          _out.pydict = NULL;
        }
#     endif

      _out.cell   = _in.cell;
      _out.freeze = _in.freeze;
      _out.energy = _in.energy;
      _out.weight = _in.weight;
      _out.scale  = _in.scale;
      _out.name   = _in.name;
      _out.atoms.resize(_in.atoms.size());
      StructureData<std::string>::t_Atoms::const_iterator i_first = _in.atoms.begin();
      StructureData<std::string>::t_Atoms::const_iterator const i_end = _in.atoms.end();
      StructureData< std::vector<std::string> >::t_Atoms::iterator i_out = _out.atoms.begin();
      for(; i_first != i_end; ++i_first, ++i_out)
      {
        *i_out = cast(*i_first, _sep); 
#       ifdef LADA_DO_PYTHON
          if(PyErr_Occurred() != NULL) return;
#       endif
      }
#     ifdef LADA_DO_PYTHON
        if(_in.pydict != NULL)
        {
          if(PyDict_Size(_in.pydict) > 0)
          {  
            PyObject* copymod = PyImport_ImportModule("copy");
            if(copymod == NULL) return;
            PyObject *deepcopystr = PyString_FromString("deepcopy");
            if(not deepcopystr) { Py_DECREF(copymod); return; }
            _out.pydict = PyObject_CallMethodObjArgs(copymod, deepcopystr, _in.pydict, NULL);
            Py_DECREF(copymod);
            Py_DECREF(deepcopystr);
          }
        }
#     endif
    }
    // \brief Casts a sequence structure to a string structure.
    // \param[in] _in: Structure to cast.
    // \param[_in] _sep: Separator between atomic type in resulting structure.
    Structure<std::string> cast( Structure< std::vector<std::string> > const &_in,
                                 std::string const &_sep)
    {
      boost::shared_ptr< StructureData<std::string> > result(new StructureData<std::string>);
      cast(*_in.get(), *result, _sep);
#     ifdef LADA_DO_PYTHON
        if(PyErr_Occurred() != NULL) return Structure<std::string>();
#     endif
      return Structure<std::string>(result);
    }
    // \brief Casts a sequence structure to a string structure.
    // \param[in] _in: Structure to cast.
    // \param[inout] _in: resulting cast structure.
    // \param[_in] _sep: Different species will be split using this separator.
    void cast( StructureData< std::vector<std::string> > const &_in, 
               StructureData<std::string> &_out, std::string const &_sep)
    {
#     ifdef LADA_DO_PYTHON
        if(_out.pydict != NULL) 
        {
          Py_DECREF(_out.pydict);
          _out.pydict = NULL;
        }
#     endif

      _out.cell   = _in.cell;
      _out.freeze = _in.freeze;
      _out.energy = _in.energy;
      _out.weight = _in.weight;
      _out.scale  = _in.scale;
      _out.name   = _in.name;
      _out.atoms.resize(_in.atoms.size());
      StructureData< std::vector<std::string> >::t_Atoms::const_iterator i_first = _in.atoms.begin();
      StructureData< std::vector<std::string> >::t_Atoms::const_iterator const i_end = _in.atoms.end();
      StructureData<std::string>::t_Atoms::iterator i_out = _out.atoms.begin();
      for(; i_first != i_end; ++i_first, ++i_out)
      {
        *i_out = cast(*i_first, _sep); 
#       ifdef LADA_DO_PYTHON
          if(PyErr_Occurred() != NULL) return;
#       endif
      }
#     ifdef LADA_DO_PYTHON
        if(_in.pydict != NULL)
        {
          if(PyDict_Size(_in.pydict) > 0)
          {  
            PyObject* copymod = PyImport_ImportModule("copy");
            if(copymod == NULL) return;
            PyObject *deepcopystr = PyString_FromString("deepcopy");
            if(not deepcopystr) { Py_DECREF(copymod); return; }
            _out.pydict = PyObject_CallMethodObjArgs(copymod, deepcopystr, _in.pydict, NULL);
            Py_DECREF(copymod);
            Py_DECREF(deepcopystr);
          }
        }
#     endif
    }
    // \brief Casts a string structure to a sequence structure.
    // \param[in] _in: Structure to cast.
    // \param[_in] _sep: Different species will be split using this separator.
    Structure< std::vector<std::string> > cast(Structure<std::string> const &_in, std::string const &_sep)
    {
      boost::shared_ptr< StructureData< std::vector<std::string> > >
        result(new StructureData< std::vector<std::string> >);
      cast(*_in.get(), *result, _sep);
#     ifdef LADA_DO_PYTHON
        if(PyErr_Occurred() != NULL) return Structure< std::vector<std::string> >();
#     endif
      return Structure< std::vector<std::string> >(result);
    }
  }
}
