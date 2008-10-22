//
//  Version: $Id$
//
#ifndef _PESCAN_INTERFACE_DIPOLE_ELEMENTS_H_ 
# define _PESCAN_INTERFACE_DIPOLE_ELEMENTS_H_
# ifdef _DIRECTIAGA 

# ifdef HAVE_CONFIG_H
#   include <config.h>
# endif

# include<vector>
# include<complex>

# include<opt/types.h>

  //! \cond
  namespace Crystal 
  {
    class Structure;
  }
  //! \endcond

  namespace Pescan
  {
    //! \cond
    class BandGap;
    //! \endcond

    //! Holds dipole moments.
    struct Dipole
    {
      enum t_SpinTransition
      {
        DOWN2DOWN, //!< Spin down to spin down transitions.
        DOWN2UP, //!< Spin down to spin up transitions.
        UP2UP, //!< Spin up to spin up transitions.
        UP2DOWN //!< Spin up to spin down transitions.
      };
      typedef std::pair<size_t, size_t> t_Band2Band;
      //! The complex dipole moment.
      std::complex<types::t_real> r[3];
      //! Indices of bands in transition.
      t_Band2Band band2band;
      //! Spin transition.
      t_SpinTransition spin2spin;
      //! Serializes a dipole moment.
      template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
        {  _ar & r; _ar & band2band; _ar & spin2spin; }
    };
    //! Prints a dipole moment.
    std::ostream& operator<<( std::ostream &_stream, const Dipole& );
    //! Computes dipole elements between valence and conduction bands.
    void dipole_elements( std::vector< Dipole > &_dipoles,
                          const BandGap& _bandgap,
                          const Crystal::Structure &_structure,
                          types::t_real _degeneracy = types::tolerance );
    //! Returns the valence-conduction the norm of the band dipole elements.
    types::t_real oscillator_strength( const std::vector<Dipole> &_dipoles );
    //! Returns the valence-conduction the norm of the band dipole elements.
    types::t_real oscillator_strength( const BandGap& _bandgap,
                                       const Crystal::Structure &_structure,
                                       types::t_real _degeneracy = types::tolerance,
                                       bool _print = false );

  }

# endif
#endif
