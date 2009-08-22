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
#include <algorithm>

#include <boost/lambda/lambda.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/fuzzy.h>

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

     inline bool operator<( Neighbor const& _a, Neighbor const& _b )
       { return _a.distance < _b.distance; }

     class Neighbors
     {
         //! Type of the list of neighbours.
         typedef std::list<Neighbor> t_Neighbors;
       public:
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
         template<class T_TYPE> 
           const_iterator begin(Crystal::TStructure<T_TYPE> const& _str) 
             { create_neighbors_list_(_str); return begin(); }
         //! returns end of neighbors list.
         const_iterator end() const { return neighbors_.end(); }

       private:
         //! Private copy Constructor.
         Neighbors(Neighbors const &_c) {}

         //! Creates list of atoms.
         template<class T_TYPE>
           void create_neighbors_list_(Crystal::TStructure<T_TYPE> const& _str);
         //! List of neighbors.
         t_Neighbors neighbors_;
     };

     template<class T_TYPE>
       void Neighbors :: create_neighbors_list_(Crystal::TStructure<T_TYPE> const& _structure)
       {
         const types::t_int N( _structure.atoms.size() );
         const types::t_int umax = nmax / _structure.atoms.size() + 1;
         neighbors_.clear();
         size_t size(0);
         
         atat::rMatrix3d const inv_cell( !_structure.cell );
         typedef Crystal::TStructure<T_TYPE> t_Structure;
         typename t_Structure::t_Atoms::const_iterator i_atom = _structure.atoms.begin();
         typename t_Structure::t_Atoms::const_iterator i_atom_end = _structure.atoms.end();
         Neighbor neighbor;
         neighbor.index = 0;
         for(; i_atom != i_atom_end; ++i_atom, ++neighbor.index ) 
         {
           atat::rVector3d const frac( inv_cell * (i_atom->pos - origin) );
           atat::rVector3d const centered
           ( 
             frac(0) - std::floor( frac(0) + 0.500000001e0 ),
             frac(1) - std::floor( frac(1) + 0.500000001e0 ),
             frac(2) - std::floor( frac(2) + 0.500000001e0 ) 
           );
           for( types::t_int x(-umax); x <= umax; ++x )
             for( types::t_int y(-umax); y <= umax; ++y )
               for( types::t_int z(-umax); z <= umax; ++z )
               {
                  neighbor.pos = _structure.cell * ( frac + atat::rVector3d(x,y,z) );
                  neighbor.distance = atat::norm( neighbor.pos );
                  if( Fuzzy::is_zero( neighbor.distance ) ) continue;
       
                  t_Neighbors :: iterator i_found 
                  (
                    std::find_if
                    (
                      neighbors_.begin(), neighbors_.end(),
                      boost::lambda::constant(neighbor) < boost::lambda::_1  
                    ) 
                  );
       
                  
                  if( i_found != neighbors_.end() ) 
                  {
                    neighbors_.insert(i_found, neighbor);
       
                    if( size < nmax ) ++size;
                    else              neighbors_.pop_back();
                    continue;
                  }
                  
                  if( size == nmax ) continue;
       
                  ++size;
                  neighbors_.push_back( neighbor );
               }
         } // loop over atoms.
       };

  } // end of Crystal namespace.
} // namespace LaDa

#endif
