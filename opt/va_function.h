//
//  Version: $Id$
//
#ifndef _FUNCTION_VA_H_
#define _FUNCTION_VA_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include <lamarck/structure.h>

#ifdef _MPI 
  #include <mpi/mpi_object.h>
#endif

namespace function
{

  //! \brief Implements a stub Virtual Atom functional.
  //! \details In other words, this functional contains a reference to a
  //!          structure which it can convert to a string of atomic occupation
  //!          variables. It should interact quite well with Minimizer::VA and
  //!          Minimizer::Beratan,
  class VirtualAtom 
  {
     //! The type of the atom container
     typedef Ising_CE::Structure::t_Atoms t_Atoms;
     //! The type of the atom
     typedef Ising_CE::Structure::t_Atom  t_Atom;

     public:
       //! see functional::Base::t_Type
       typedef types::t_real t_Type;
       //! see functional::Base::t_Container
       typedef std::vector< t_Type >  t_Container;

     protected:
       Ising_CE::Structure& structure;
       t_Container va_vars;

     public:
       //! Constructor and Initializer
       VirtualAtom   ( Ising_CE::Structure &_str )
                   : structure( _str )
        { va_vars.reserve( _str.atoms.size() ); }
       //! Copy Constructor
       VirtualAtom   ( const VirtualAtom &_c )
                   : structure( _c.structure ), va_vars( _c.va_vars ) {}
        
       // Simple constainer behaviors required by Minimizer::VA and
       // Minimizer::Beratan

       //! Returns the size of VirtualAtom::va_vars.
       types::t_unsigned size() const { return va_vars.size(); }
       //! Returns an iterator to the first \e virtual variable (atomic occupation).
       t_Container::iterator begin() { return va_vars.begin(); }
       //! \brief Returns an iterator to one past the last \e virtual variable
       //!        (atomic occupation).
       t_Container::iterator end() { return va_vars.end(); }
       //! \brief Returns a constant iterator to the first \e virtual variable
       //!        (atomic occupation).
       t_Container::const_iterator begin() const { return va_vars.begin(); }
       //! \brief Returns a constant iterator to one past the last \e virtual
       //!        variable (atomic occupation).
       t_Container::const_iterator end() const { return va_vars.end(); }

       // Now truly "functional" stuff.
       
       //! Initializes the variables with respect to Functional::structure.
       bool init();

     protected:
       //! Transfers occupations from VirtualAtom::va_vars to Functional::structure.
       void unpack_variables();
  };

  inline bool VirtualAtom :: init()
  {
    va_vars.clear();
    va_vars.reserve( structure.atoms.size() );
    t_Atoms :: const_iterator i_atom = structure.atoms.begin();
    t_Atoms :: const_iterator i_atom_end = structure.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
      if( not ( i_atom->freeze & t_Atom::FREEZE_T ) ) 
        va_vars.push_back( i_atom->type );

    return not va_vars.empty();
  }

  inline void VirtualAtom :: unpack_variables()
  {
    t_Atoms :: iterator i_atom = structure.atoms.begin();
    t_Atoms :: iterator i_atom_end = structure.atoms.end();
    t_Container :: const_iterator i_var = va_vars.begin();
    for(; i_atom != i_atom_end; ++i_atom )
    {
      if( i_atom->freeze & t_Atom::FREEZE_T ) continue;

      i_atom->type = *i_var > t_Type(0) ? 1.0: -1.0;
      ++i_var;
    }
  }

} // namespace vff 

#endif // _VFF_FUNCTIONAL_H_
