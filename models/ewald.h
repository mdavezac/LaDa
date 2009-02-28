//
//  Version: $Id: vff.h 895 2008-12-22 02:04:18Z davezac $
//
#ifndef _MODELS_EWALD_FUNCTIONAL_H_
#define _MODELS_EWALD_FUNCTIONAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>

#include <tinyxml/tinyxml.h>

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <opt/types.h>


namespace LaDa
{
  namespace Models
  {
    //! \cond
    class Clj;
    void read_fortran_input( Clj&, const boost::filesystem::path& );
    //! \endcond

    //! \brief Interface to fortran ewald sum routine.
    //! \details Ad-hoc implementation for Clj functional.
    class Ewald
    {
      friend  void read_fortran_input( Clj&, const boost::filesystem::path&);
      public:
        //! Type of the key used in the map.
        typedef std::string Key;
        //! Type of the container of atomic charges.
        typedef std::map< Key, types::t_int > t_Charges;
        //! Argument type.
        typedef Crystal :: TStructure< std::string > t_Arg;
        //! Return type.
        typedef types::t_real t_Return;

        //! Contains all atomic species.
        t_Charges charges;

        //! Constructor and Initializer
        Ewald() {}
        //! Copy Constructor
        Ewald( const Ewald &_c ) : charges(_c.charges), cutoff_(_c.cutoff_) {}
        //! Destructor.
        ~Ewald() {}

        //! \brief computes energy and stress and forces.
        //! \param[in] _in input structure, with reduced coordinates.
        //! \param[inout] _out stress and forces. Both are \e added to existing stress and forces.
        t_Return energy(const t_Arg& _in, t_Arg& _out) const;

        //! Loads parameters from XML.
        bool Load( const TiXmlElement& _node );

        //! Sets real-space cutoff.
        types::t_real get_rcutoff() const { return cutoff_; }
        //! Sets real-space cutoff.
        void set_rcutoff(types::t_real _r) { cutoff_ = _r; }

      protected:
        //! Cutoff.
        types::t_real cutoff_;
    };

  } // namespace CLJ.
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_
