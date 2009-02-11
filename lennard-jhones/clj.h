//
//  Version: $Id: vff.h 895 2008-12-22 02:04:18Z davezac $
//
#ifndef _CLJ_FUNCTIONAL_H_
#define _CLJ_FUNCTIONAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <opt/types.h>
#include <opt/debug.h>


namespace LaDa
{
  //! Lennard-Jhones+Coulomb functionals.
  namespace CLJ
  {
    //! The adapter for Coulomb + Lennard-Jhones fortran functional itself.
    class Clj
    {
      //! A structure holding the specie types.
      struct Specie;
      //! Type of the key used in the map.
      typedef std::string Key;
      //! Type of the container of atomic species.
      typedef std::map< Key, Specie > t_Species;
      public:
        //! Argument type.
        typedef Crystal :: TStructure< std::string > t_Arg;
        //! Return type.
        typedef types::t_real t_Return;
        //! Constructor and Initializer
        Clj() {}
        //! Copy Constructor
        Clj( const Clj &_c ) {}
        //! \brief Destructor
        Clj() {}

        //! computes energy and stress and forces.
        t_Return energy(const t_Arg& _in, t_Arg& _out) const;

        //! Add a specie.
        void add_specie( Specie::t_Type _type, Specie::t_Charge _charge, Specie::t_Radius _radius );
        //! Get a specie.
        const Specie& get_specie( size_t n ) const 
          { __ASSERT( n >= species_.size(), "index out-of-range." ) return species_[_n]; }
        //! Get a specie.
        Specie& get_specie( size_t n ) 
          { __ASSERT( n >= species_.size(), "index out-of-range." ) return species_[_n]; }
        //! Get nb of species.
        size_t get_nbspecies() const  { return species_.size(); }

        //! Initializes fortran module with values of this instance.
        void init () const;

      protected:
        //! Computes Lennard-Jhones term only.
        t_Return lennard_jhones( const t_Arg& _in, t_Arg &_out );
       
        //! Contains all atomic species.
        t_Species species_;
        //! Epsilon.
        types::t_real epsilon;
        //! scale for cutoff radius for lennard-jhones potential.
        types::t_real rcut_lj;
        //! scale cutoff radius for Ewald sum.
        types::t_real rcut_ewald;
        //! Scale of hard-sphere atomic radii.
        types::t_real radii_scale;
        //! Percent of r^12 term in lennard-jhones.
        types::t_real PEGS;
    };

    struct Clj :: Specie
    {
      //! Type of the atomic mass.
      typedef std::string t_Type;
      //! Type of the atomic charge.
      typedef types::t_int t_Charge;
      //! Type of the atomic radius.
      typedef types::t_real t_Radius;
      //! Atomic Mass.
      t_Type type;
      //! Atomic Charge.
      t_Charge charge;
      //! Atomic radius.
      t_Radius radius;
      //! Constructor.
      Specie( t_Type _t, t_Charge _c, t_Radius _r ) : type( _t ), charge( _c ), radius( _r ) {}
      //! Copy Constructor.
      Specie( const Specie& _t ) : type( _c.type ), charge( _c.charge ), radius( _c.radius ) {}
    };

  } // namespace CLJ.
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_
