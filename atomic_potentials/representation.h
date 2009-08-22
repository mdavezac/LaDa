//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_REPRESENTATION_H_
#define LADA_ATOMIC_POTENTIAL_REPRESENTATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <utility>
#include <vector>

#include "numeric_types.h"

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

    //! Compares to variable sets.
    bool operator==(VariableSet const& _a, VariableSet const &_b) 
      { return _a.variables == _b.variables; }


    //! Representation of a structure a set of symmetry-equivalents set of variables.
    class Representation
    {
        //! Type of a set of representations.
        typedef std::vector<VariableSet> t_Sets;
      public:
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

      private:
        //! Adds a variable set to this representation.
        void add_( VariableSet const &_rep );
        //! Symmetry-equivalent sets of variables.
        t_Sets sets_;
    };


  } // namespace atomic_potential
} // namespace LaDa
#endif
