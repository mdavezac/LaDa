namespace LaDa
{
  namespace crystal
  {
    extern "C" 
    {
      //! Returns a representation of the object.
      static PyObject* structure_repr(StructureData* _self);
      //! Returns a deepcopy of the atom.
      static PyObject* structure_copy(StructureData* _self)
        { return (PyObject*) PyStructure_Copy(_self, NULL); }
      //! Implements deepcopy.
      static PyObject* structure_deepcopy(StructureData* _self, PyObject* _memo)
        { return (PyObject*) PyStructure_Copy(_self, _memo); }
      //! Implements shallow copy.
      static PyObject* structure_shallowcopy(StructureData* _self)
        { Py_INCREF(_self); return (PyObject*)_self; }
      //! Returns a dictionary with same info as atom.
      static PyObject* structure_to_dict(StructureData* _self);
      //! Implements getstate for pickling.
      static PyObject* structure_getstate(StructureData* _self);
      //! Implements setstate for pickling.
      static PyObject* structure_setstate(StructureData* _self, PyObject *_dict);
      //! Implements reduce for pickling.
      static PyObject* structure_reduce(StructureData* _self);
      //! Implements add atom.
      static PyObject* structure_add_atom(StructureData* _self, PyObject* _args, PyObject* _kwargs);
    }

    //! Returns a representation of the object.
    PyObject* structure_repr(StructureData* _self)
    {
      std::string name(_self->ob_type->tp_name);
      name = name.substr(name.rfind('.')+1);
      std::ostringstream result;
      // Create structure variable first.
      result << "structure = " << name << "( "
             <<         _self->cell(0, 0)
             << ", " << _self->cell(0, 1) 
             << ", " << _self->cell(0, 2) 
             << ",\\\n" << std::string(name.size()+14, ' ')
             <<         _self->cell(1, 0)
             << ", " << _self->cell(1, 1) 
             << ", " << _self->cell(1, 2) 
             << ",\\\n" << std::string(name.size()+14, ' ')
             <<         _self->cell(2, 0)
             << ", " << _self->cell(2, 1) 
             << ", " << _self->cell(2, 2)
             << ",\\\n" << std::string(name.size()+14, ' ')
             << "scale=" << _self->scale;
          
      // Including python dynamic attributes.
      if(_self->pydict != NULL)
      {
        if(PyDict_Size(_self->pydict) > 0)
        {
          PyObject *key, *value;
          Py_ssize_t pos = 0;
          while (PyDict_Next(_self->pydict, &pos, &key, &value)) 
          {
            PyObject* repr = PyObject_Repr(value);
            if(repr == NULL) return NULL;
            result << ", " << PyString_AsString(key);
            if(PyErr_Occurred() != NULL) {Py_DECREF(repr); return NULL;}
            result << "=" << PyString_AsString(repr);
            Py_DECREF(repr);
            if(PyErr_Occurred() != NULL) return NULL;
          }
        }
      }
      result << " )";
      // Then add atoms.
      if(_self->atoms.size() > 0)
      {
        std::vector<Atom>::const_iterator i_first = _self->atoms.begin();
        std::vector<Atom>::const_iterator const i_end = _self->atoms.end();
        result << "\n" << "structure.add_atom";
        {
          python::Object const atomstr = PyObject_Repr(*i_first);
          if(not atomstr) return NULL;
          std::string const atom = PyString_AsString(atomstr);
          if(atom.empty()) return NULL;
          else if(atom.substr(0, 4) == "Atom") result << atom.substr(4);
          else result << "(" << atom << ")";
        }
        std::string const empty("\\\n" + std::string(18, ' '));
        for(++i_first; i_first != i_end; ++i_first)
        {
          python::Object atomstr = PyObject_Repr(*i_first);
          if(not atomstr) return NULL;
          std::string atom = PyString_AsString(atomstr);
          if(atom.empty()) return NULL;
          else if(atom.substr(0, 4) == "Atom") result << empty << atom.substr(4);
          else result << empty << "(" << atom << ")";
        }
      }
      return PyString_FromString(result.str().c_str());
    }

    // Creates dictionary from atom with shallow copies.
    PyObject *structure_to_dict(StructureData* _self)
    {
      python::Object result = PyDict_New();
      if(not result) return NULL;
  
      python::Object const cell = structure_getcell(_self, NULL);
      if(not cell) return NULL;
      if(PyDict_SetItemString(result, "cell", cell) < 0) return NULL;
      python::Object const scale = structure_getscale(_self, NULL);
      if(not scale) return NULL;
      if(PyDict_SetItemString(result, "scale", scale) < 0) return NULL;
  
      std::vector<Atom>::const_iterator i_atom = _self->atoms.begin();
      std::vector<Atom>::const_iterator const i_end = _self->atoms.end();
      char mname[] = "to_dict";
      for(long i(0); i_atom != i_end; ++i_atom, ++i)
      {
        // Gets dictionary description.
        python::Object const item = PyObject_CallMethod((PyObject*)&(*i_atom), mname, NULL);
        if(not item) return NULL;
        // Then create pyobject index.
        python::Object const index = PyInt_FromLong(i);
        if(not index) return NULL;
        // finally, adds to dictionary.
        if(PyDict_SetItem(result, index, item) < 0) return NULL;
      }
      // Merge attribute dictionary if it exists.
      if(_self->pydict != NULL and PyDict_Merge(result, _self->pydict, 1) < 0) return NULL;
  
      return result.release();
    }

    // Implements __reduce__ for pickling.
    PyObject* structure_reduce(StructureData* _self)
    {
      // Creates return tuple of three elements.
      python::Object type = PyObject_Type((PyObject*)_self);
      if(not type) return NULL;
      // Second element is a null tuple, argument to the callable type above.
      python::Object tuple = PyTuple_New(0);
      if(not tuple) return NULL;
      // Third element is the state of this object.
      char getstate[] = "__getstate__";
      python::Object state = PyObject_CallMethod((PyObject*)_self, getstate, NULL, NULL);
      if(not state) return NULL;

      return PyTuple_Pack(3, type.borrowed(), tuple.borrowed(), state.borrowed());
    }

    // Implements getstate for pickling.
    PyObject* structure_getstate(StructureData* _self)
    {
      // get cell attribute.
      python::Object cell = structure_getcell(_self, NULL);
      if(not cell) return NULL;
      // get scale attribute.
      python::Object scale = structure_getscale(_self, NULL);
      if(not scale) return NULL;
      // get python dynamic attributes.
      python::Object dict = _self->pydict == NULL ? Py_None: PyDict_New();
      if(not dict) return NULL;
      if(_self->pydict != NULL and PyDict_Merge(dict, _self->pydict, 0) < 0) return NULL;
      // get all atoms.
      python::Object atoms = PyTuple_New(_self->atoms.size());
      if(not atoms) return NULL;
      std::vector<Atom>::const_iterator i_first = _self->atoms.begin();
      std::vector<Atom>::const_iterator const i_end = _self->atoms.end();
      char mname[] = "__getstate__";
      for(Py_ssize_t i(0); i_first != i_end; ++i_first, ++i) 
      {
        python::Object item = PyObject_CallMethod((PyObject*)&(*i_first), mname, NULL);
        if(not item) return NULL;
        if(PyTuple_SET_ITEM((PyObject*)atoms, i, item.release()) < 0) return NULL;
      }

      return PyTuple_Pack(4, cell.borrowed(), scale.borrowed(), dict.borrowed(), atoms.borrowed());
    }

    // Implements setstate for pickling.
    PyObject* structure_setstate(StructureData* _self, PyObject *_tuple)
    {
      if(not PyTuple_Check(_tuple))
      {
        LADA_PYERROR(TypeError, "Expected state to be a tuple.");
        return NULL;
      }
      if(PyTuple_Size(_tuple) != 4)
      {
        LADA_PYERROR(TypeError, "Expected state to be a 4-tuple.");
        return NULL;
      }
      // first cell and scale.
      if(structure_setcell(_self, PyTuple_GET_ITEM(_tuple, 0), NULL) < 0) return NULL;
      if(structure_setscale(_self, PyTuple_GET_ITEM(_tuple, 1), NULL) < 0) return NULL;

      PyObject *atoms = PyTuple_GetItem(_tuple, 3);
      if(not PyTuple_Check(atoms)) 
      {
        LADA_PYERROR(TypeError, "Expected atom container to be a tuple.");
        return NULL;
      }
      // then atoms.
      Py_ssize_t const N(PyTuple_Size(atoms));
      _self->atoms.reserve(N);
      for(Py_ssize_t i(0); i < N; ++i)
      {
        AtomData* atom = (AtomData*)PyTuple_GET_ITEM(atoms, i);
        if(not PyAtom_Check(atom)) 
        {
          LADA_PYERROR(TypeError, "Expected an atom when unpickling.");
          return NULL;
        }
        Py_INCREF(atom);
        _self->atoms.push_back(Atom(atom));
      }
      // finally, dictionary, so we can return without issue on error.
      PyObject *dict = PyTuple_GET_ITEM(_tuple, 2);
      if(dict == Py_None) { Py_RETURN_NONE; }
      if(_self->pydict == NULL)
      {
        _self->pydict = PyDict_New();
        if(_self->pydict == NULL) return NULL;
      }
      if(PyDict_Merge(_self->pydict, dict, 0) < 0) return NULL;
      Py_RETURN_NONE;
    }

    // Implements add atom.
    static PyObject* structure_add_atom(StructureData* _self, PyObject* _args, PyObject* _kwargs)
    {
      // Check first that _args is not a tuple containing an atom.
      if(PyTuple_Size(_args) == 1)
      {
        AtomData* wrapper = (AtomData*)PyTuple_GET_ITEM(_args, 0);
        if(PyAtom_Check(wrapper)) 
        {
          if(_kwargs != NULL)
          {
            LADA_PYERROR(TypeError, "Cannot insert an atom and motify in-place.");
            return NULL;
          }
          if(not wrapper)
          {
            LADA_PYERROR(internal, "Should never find an empty atom. Internal bug.");
            return NULL;
          }
          Py_INCREF(wrapper);
          _self->atoms.push_back(Atom(wrapper));
          return PyObject_GetAttrString((PyObject*)_self, "add_atom");
        }
      }
      // create new atom and its wrapper.
      Atom atom(_args, _kwargs);
      if(not atom) return NULL;
      // Add it to the container.
      _self->atoms.push_back(atom);
      // Finally, returns this very function for chaining.
      return PyObject_GetAttrString((PyObject*)_self, "add_atom");
    }
  }
}
