//
//  Version: $Id: confsplit.h 844 2008-11-08 01:22:54Z davezac $
//
#ifndef _LADA_CRYSTAL_NEIGHBORS_H_
#define _LADA_CRYSTAL_NEIGHBORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <utility>

#include <opt/types.h>
#include <opt/debug.h>

#include "structure.h"

namespace LaDa
{
  namespace Crystal
  {
     //! Holds neigbor data.
     struct Neighbor
     {
       //! Index to atom in structure.
       size_t index;
       //! Position with respect to atom.
       atat::rVector3d pos;
       //! Distance from atom.
       types::t_real distance;
     };

     class Neighbors
     {
         //! Type of the list of neighbours.
         typedef std::list<Neighbor> t_Neighbors;
       public:
         //! Type of the crystal structure used.
         typedef Crystal::Structure t_Structure;
         //! Iterator over neighbors.
         typedef t_Neighbors::const_iterator const_iterator;
         //! The number of first neighbors to compute.
         size_t nmax;
         //! The origin from which to compute neighbors.
         atat::rVector3d origin;

         //! Constructor.
         Neighbors   (size_t _nmax = 0, atat::rVector3d const &_vec = atat::rVector3d(0,0,0) )
                   : nmax(_nmax), origin(_vec) {};
         //! returns iterator to first neighbor list.
         const_iterator begin() const { return neighbors_.begin(); }
         //! constructs first neighbor list and returns first iterator.
         const_iterator begin(t_Structure const& _str) 
           { create_neighbors_list_(_str); return begin(); }
         //! returns end of neighbors list.
         const_iterator end() const { return neighbors_.end(); }

       private:
         //! Private copy Constructor.
         Neighbors(Neighbors const &_c) {}

         //! Creates list of atoms.
         void create_neighbors_list_(t_Structure const& _str);
         //! List of neighbors.
         t_Neighbors neighbors_;
     };


  } // end of Crystal namespace.
} // namespace LaDa

#endif
