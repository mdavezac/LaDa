#ifndef LADA_VFF_EDGE_DATA_H
#define LADA_VFF_EDGE_DATA_H

#include "LaDaConfig.h"

//! \def check_structure(object)
//!      Returns true if an object is an edge or subtype.
#define PyEdgeData_Check(object) PyObject_TypeCheck(object, LaDa::vff::edge_type())
//! \def PyNodeData_CheckExact(object)
//!      Returns true if an object is a edge.
#define PyEdgeData_CheckExact(object) object->ob_type == LaDa::vff::edge_type()
      

namespace LaDa 
{
  namespace vff
  {
    extern "C" 
    {
      struct NodeData;
      //! \brief Describes an edge between first-neighbor nodes.
      //! This object is not meant to be returned directly in python.
      //! Rather, it is exists as a python object only so that python can do
      //! the memory management.
      struct EdgeData
      {
        PyObject_HEAD 
        //! Whether a periodic translation is involved.
        bool do_translate;
        //! Translation.
        math::rVector3d translation;
        //! Node A of bond.
        NodeData* a; 
        //! Node B of bond.
        NodeData* b; 
      };
      //! Creates a new edge.
      EdgeData* PyEdge_New();
      //! Returns pointer to node type.
      PyTypeObject* edge_type();
    }
    //! Returns tuple with both a and b.
    PyObject* edge_to_tuple(EdgeData* _data);
    //! Returns tuple with only one of a and b.
    PyObject* edge_to_tuple(NodeData* _self, EdgeData* _data);

  } // namespace vff

} // namespace LaDa

#endif

