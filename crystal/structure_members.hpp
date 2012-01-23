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
      //! Implements structure affine transformation.
      static PyObject* structure_transform(StructureData* _self, PyObject* _args);
    }

    //! Returns a representation of the object.
    PyObject* structure_repr(StructureData* _self)
    {
      std::string name(_self->ob_type->tp_name);
      name = name.substr(name.rfind('.')+1);
      std::ostringstream result;
      // Create structure variable first.
      result << name << "( "
             <<         _self->cell(0, 0)
             << ", " << _self->cell(0, 1) 
             << ", " << _self->cell(0, 2) 
             << ",\\\n" << std::string(name.size()+2, ' ')
             <<         _self->cell(1, 0)
             << ", " << _self->cell(1, 1) 
             << ", " << _self->cell(1, 2) 
             << ",\\\n" << std::string(name.size()+2, ' ')
             <<         _self->cell(2, 0)
             << ", " << _self->cell(2, 1) 
             << ", " << _self->cell(2, 2)
             << ",\\\n" << std::string(name.size()+2, ' ')
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
        result << "\\\n  .add_atom";
        {
          python::Object const atomstr = PyObject_Repr(i_first->borrowed());
          if(not atomstr) return NULL;
          std::string const atom = PyString_AsString(atomstr.borrowed());
          if(atom.empty()) return NULL;
          else if(atom.substr(0, 4) == "Atom") result << atom.substr(4);
          else result << "(" << atom << ")";
        }
        std::string const prefix("\\\n  .add_atom");
        for(++i_first; i_first != i_end; ++i_first)
        {
          python::Object atomstr = PyObject_Repr(i_first->borrowed());
          if(not atomstr) return NULL;
          std::string atom = PyString_AsString(atomstr.borrowed());
          if(atom.empty()) return NULL;
          else if(atom.substr(0, 4) == "Atom") result << prefix << atom.substr(4);
          else result << prefix << "(" << atom << ")";
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
      if(PyDict_SetItemString(result.borrowed(), "cell", cell.borrowed()) < 0) return NULL;
      python::Object const scale = structure_getscale(_self, NULL);
      if(not scale) return NULL;
      if(PyDict_SetItemString(result.borrowed(), "scale", scale.borrowed()) < 0) return NULL;
  
      std::vector<Atom>::const_iterator i_atom = _self->atoms.begin();
      std::vector<Atom>::const_iterator const i_end = _self->atoms.end();
      char mname[] = "to_dict";
      for(long i(0); i_atom != i_end; ++i_atom, ++i)
      {
        // Gets dictionary description.
        python::Object const item = PyObject_CallMethod(i_atom->borrowed(), mname, NULL);
        if(not item) return NULL;
        // Then create pyobject index.
        python::Object const index = PyInt_FromLong(i);
        if(not index) return NULL;
        // finally, adds to dictionary.
        if(PyDict_SetItem(result.borrowed(), index.borrowed(), item.borrowed()) < 0) return NULL;
      }
      // Merge attribute dictionary if it exists.
      if(_self->pydict != NULL and PyDict_Merge(result.borrowed(), _self->pydict, 1) < 0) return NULL;
  
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
      python::Object iterator = structureiterator_create(_self);
      if(not iterator) return NULL;

      return PyTuple_Pack(4, type.borrowed(), tuple.borrowed(), state.borrowed(), iterator.borrowed());
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
      python::Object dict = _self->pydict == NULL ? python::Object::acquire(Py_None): PyDict_New();
      if(not dict) return NULL;
      if(_self->pydict != NULL and PyDict_Merge(dict.borrowed(), _self->pydict, 0) < 0) return NULL;

      return PyTuple_Pack(3, cell.borrowed(), scale.borrowed(), dict.borrowed());
    }

    // Implements setstate for pickling.
    PyObject* structure_setstate(StructureData* _self, PyObject *_tuple)
    {
      if(not PyTuple_Check(_tuple))
      {
        LADA_PYERROR(TypeError, "Expected state to be a tuple.");
        return NULL;
      }
      if(PyTuple_Size(_tuple) != 3)
      {
        LADA_PYERROR(TypeError, "Expected state to be a 4-tuple.");
        return NULL;
      }
      // first cell and scale.
      if(structure_setcell(_self, PyTuple_GET_ITEM(_tuple, 0), NULL) < 0) return NULL;
      if(structure_setscale(_self, PyTuple_GET_ITEM(_tuple, 1), NULL) < 0) return NULL;

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
          _self->atoms.push_back(Atom::acquire_((PyObject*)wrapper));
          Py_INCREF(_self);
          return (PyObject*)_self;
        }
      }
      // create new atom and its wrapper.
      Atom atom(_args, _kwargs);
      if(not atom) return NULL;
      // Add it to the container.
      _self->atoms.push_back(atom);
      // Finally, returns this very function for chaining.
      Py_INCREF(_self);
      return (PyObject*)_self;
    }

    static PyObject* structure_transform(StructureData* _self, PyObject* _args)
    {
      Eigen::Matrix<types::t_real, 4, 3> op;
      if(not python::convert_to_matrix(_args, op)) return NULL;
      _self->cell = op.block<3,3>(0,0) * _self->cell;
      std::vector<Atom>::iterator i_atom = _self->atoms.begin();
      std::vector<Atom>::iterator i_atom_end = _self->atoms.end();
      for(; i_atom != i_atom_end; ++i_atom)
        i_atom->pos() = op.block<3,3>(0,0) * i_atom->pos() + ~op.block<1, 3>(3, 0);
      Py_RETURN_NONE;
    }
  }
}
