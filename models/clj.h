//
//  Version: $Id: vff.h 895 2008-12-22 02:04:18Z davezac $
//
#ifndef _CLJ_FUNCTIONAL_H_
#define _CLJ_FUNCTIONAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>

#include <tinyxml/tinyxml.h>

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <opt/types.h>
#include <opt/debug.h>

#include "ewald.h"
#include "lennard-jones.h"

namespace LaDa
{
  //! Empirical models.
  namespace Models
  {
    //! The adapter for Coulomb + Lennard-Jones fortran functional itself.
    class Clj : protected Ewald, protected LennardJones
    {
      public:
        //! Argument type.
        typedef Crystal :: TStructure< std::string > t_Arg;
        //! Return type.
        typedef types::t_real t_Return;
        //! Constructor and Initializer
        Clj() {}
        //! Copy Constructor
        Clj( const Clj &_c ) {}
        //! Destructor.
        ~Clj() {}

        //! \brief computes energy and stress and forces.
        t_Return energy(const t_Arg& _in, t_Arg& _out) const;

        //! Loads parameters from XML.
        bool Load( const TiXmlElement& _node );

      protected:
        //! Checks coherency between Ewald and LennardJones.
        void check_coherency() const;
    };

    //! Creates a CLJ functional using same input as fortran.
    void read_fortran_input( Clj &_clj, boost::filesystem &_path );

  } // namespace CLJ.
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_
