#ifndef LADA_ATOMIC_POTENTIAL_REPRESENTATION_H_
#define LADA_ATOMIC_POTENTIAL_REPRESENTATION_H_

#include "LaDaConfig.h"

#include <iostream>
#include <utility>
#include <vector>

#include <opt/debug.h>

#include "numeric_types.h"
#include <math/fuzzy.h>

namespace LaDa
{
  // Forward declaration.
  //! \cond
  namespace Crystal
  {
    template<class T_TYPE> class TStructure;
  }
  //! \endcond

 
  //! Contains classes for generalized atomic potentials.
  namespace atomic_potential
  {
    //! Unique set of variables representing one possible symmetry-equivalent structure.
    struct VariableSet
    {
      //! Type of a variable.
      typedef std::pair<numeric_type, specie_type> t_Variable;
      //! Type of a representation.
      typedef std::vector<t_Variable> t_Variables;
      
      //! Weight of this unique variable set in the representation.
      numeric_type weight;
      //! Varibale set.
      t_Variables variables;
    };

    //! Dumps a variable to a stream.
    inline std::ostream& operator<<( std::ostream& _stream, VariableSet::t_Variable const &_var )
      { return _stream << "(" << (math::is_zero(_var.first) ? 0: _var.first)
                       << ", " << _var.second << ")"; }

    //! Compares to variable sets.
    bool operator==(VariableSet const& _a, VariableSet const &_b);
    //! Dumps a variable set to a stream.
    std::ostream& operator<<( std::ostream& _stream, VariableSet const &_rep );


    //! Representation of a structure a set of symmetry-equivalents set of variables.
    class Representation
    {
        //! Type of a set of representations.
        typedef std::vector<VariableSet> t_Sets;
      public:
        //! Coordinate Type.
        enum CoordinateSystem
        {
          cartesian, //! Cartesian coordinates.
          spherical  //! Spherical coordinates.
        };
        //! Switch from one coordinate system to another.
        static CoordinateSystem coord_system;
        //! Iterator over the representations.
        typedef t_Sets::const_iterator const_iterator;
  
        //! Constructor.
        Representation(Crystal::TStructure<std::string> const &_str, size_t _natoms);
        //! Copy Constructor.
        Representation(Representation const& _c): sets_(_c.sets_) {}

        //! Returns a const iterator to the first variable set.
        const_iterator begin() const { return sets_.begin(); }
        //! Returns a const iterator to the past-the-end variable set.
        const_iterator end() const { return sets_.end(); }
        //! Returns  set \a _n.
        t_Sets::value_type& operator[](size_t _n) 
        {
          LADA_ASSERT( _n < sets_.size(), "Index out-of-range.\n");
          return sets_[_n];
        }
        //! Returns  set \a _n.
        t_Sets::value_type const& operator[](size_t _n) const
        {
          LADA_ASSERT( _n < sets_.size(), "Index out-of-range.\n");
          return sets_[_n];
        }
        //! Returns number of coordinates.
        size_t nb_coords() const { return sets_.front().variables.size(); }
        //! Returns number of atoms.
        size_t nb_atoms() const { return (sets_.front().variables.size() + 5)/3; }
        //! Returns number of atoms.
        size_t size() const { return sets_.size(); }

      private:
        //! Create a the representation.
        void create_(Crystal::TStructure<std::string> const &_str, size_t _natoms);
        //! Adds a variable set to this representation.
        void add_( VariableSet const &_rep );
        //! Symmetry-equivalent sets of variables.
        t_Sets sets_;
    };

    //! Dumps a representation to a stream.
    std::ostream& operator<<( std::ostream& _stream, Representation const &_rep );

  } // namespace atomic_potential
} // namespace LaDa
#endif
