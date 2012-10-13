#ifndef LADA_VFF_NODE_DATA_H
#define LADA_VFF_NODE_DATA_H

#include "LaDaConfig.h"

#include <vector>

#include <crystal/atom/atom.h>

//! \def PyStructure_Check(object)
//!      Returns true if an object is a node or subtype.
#define PyNodeData_Check(object) PyObject_TypeCheck(object, LaDa::vff::node_type())
//! \def PyNodeData_CheckExact(object)
//!      Returns true if an object is a node.
#define PyNodeData_CheckExact(object) object->ob_type == LaDa::vff::node_type()
      

namespace LaDa 
{
  namespace vff
  {
    extern "C" 
    {
      class EdgeData;
      //! Describes a node in a first neighbor net.
      struct NodeData
      {
        PyObject_HEAD 
        //! Holds list of weak pointers.
        PyObject *weakreflist;
        //! Holds reference to other bonds.
        std::vector<EdgeData*> bonds;
        //! Holds reference to other an atom.
        crystal::Atom center;
        //! Index of the atom in the structure.
        types::t_unsigned index;
      };
      //! Creates a new node.
      NodeData* PyNodeData_New();
      // Creates a new structure with a given type.
      NodeData* PyNode_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Adds an edge between two bonds. 
      bool PyNode_AddEdge(NodeData* _a, NodeData* _b, math::rVector3d &_trans);
      // Returns pointer to node type.
      PyTypeObject* node_type();

      //! Returns pointer to bond iterator type.
      PyTypeObject* bonditerator_type();
    }
  } // namespace vff

} // namespace LaDa

#endif
