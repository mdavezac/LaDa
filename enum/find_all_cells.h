//
//  Version: $Id$
//
#ifndef LADA_ENUM_FIND_ALL_CELLS_H_
#define LADA_ENUM_FIND_ALL_CELLS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <iostream>
#include <boost/shared_ptr.hpp>

#include <opt/debug.h>
#include <opt/types.h>
#include <math/fuzzy.h>


namespace LaDa
{
  //! \cond
  namespace Crystal
  {
    class Lattice;
  }
  //! \endcond

  namespace enumeration
  {
    //! \brief Finds all inequivalent supercells of size \a _nmax.
    boost::shared_ptr< std::vector<math::rMatrix3d> >
      find_all_cells( Crystal::Lattice const &_lattice, size_t _nmax);

    struct SmithGroup
    {
      //! left transform and hermite normal form.
      struct Supercell
      {
        //! left tranform;
        math::rMatrix3d transform;
        //! Hermite form;
        math::rMatrix3d hermite;
        //! Constructor.
        Supercell () { transform.zero(); hermite.zero(); }
        //! Constructor.
        Supercell   ( math::rMatrix3d const & _transform, math::rMatrix3d const & _hermite )
                  : transform(_transform), hermite(_hermite) {}
        //! Copy Constructor.
        Supercell   (Supercell const &_c)
                  : transform(_c.transform), hermite(_c.hermite) {}
        //! Comparison operator.
        bool operator==(Supercell const& _a) const
        {
          for( size_t i(0); i < 3; ++i )
            for( size_t j(0); j < 3; ++j )
              if( not math::is_zero(hermite(i,j)-_a.hermite(i,j)) ) return false;
          return true;
        }
      };
      //! Supercells container.
      typedef std::vector<Supercell> t_Supercells;
 
      //! Constructor
      SmithGroup( math::iVector3d const &_id = math::iVector3d(0,0,0) ) : smith(_id) {}
      //! Copy constructor
      SmithGroup( SmithGroup const &_c ) : smith(_c.smith), supercells(_c.supercells) {}
 
      //! Compare smith group by id.
      bool operator==( math::iVector3d const &_vec ) const { return smith == _vec; }
      //! Compare smith group by id.
      bool operator==( SmithGroup const &_sg ) const
      {
        if( not operator==(_sg.smith) ) return false; 
        return supercells == _sg.supercells;
      }
 
      //! Quotient group.
      math::iVector3d smith;
      //! left transform and hermite normal form.
      std::vector<Supercell> supercells;
    };

    inline std::ostream& operator<<(std::ostream &_stream, SmithGroup const& _sg)
    {
      return _stream << "Smith Group (" << _sg.smith << ") of "
                     << _sg.supercells.size() << " Matrices.";
    }
    inline std::ostream& operator<<(std::ostream &_stream, SmithGroup::Supercell const& _sg)
    {
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j <= i; ++j )
          _stream << _sg.hermite(i,j) << " ";
      return _stream;
    }
 
    //! Creates a vector of supercells structure according to their Smith normal form.
    boost::shared_ptr< std::vector<SmithGroup> >
      create_smith_groups( Crystal::Lattice const &_lattice,
                           boost::shared_ptr< std::vector<math::rMatrix3d> > const & _s );
 
    //! Creates a vector of supercells structure according to their Smith normal form.
    inline boost::shared_ptr< std::vector<SmithGroup> >
      create_smith_groups( Crystal::Lattice const &_lattice, size_t _nmax)
       { return create_smith_groups( _lattice, find_all_cells( _lattice, _nmax) ); }
  }
}

#endif
