//
//  Version: $Id: interface.h 816 2008-10-17 01:29:20Z davezac $
//
#if undefined( _PESCAN_INTERFACE_DIPOLE_ELEMENTS_H_ ) && defined( _DIRECTIAGA ) 
# define _PESCAN_INTERFACE_DIPOLE_ELEMENTS_H_

# ifdef HAVE_CONFIG_H
#   include <config.h>
# endif

# include<cmath>

# include<opt/types.h>

  // declares fortran interface
  //! \cond
  extern "C"
  {
    void FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( MPI_Fint * );
    void FC_FUNC_(iaga_call_escan, IAGA_CALL_ESCAN)( int* );
    void FC_FUNC_(iaga_get_eigenvalues, IAGA_GET_EIGENVALUES)( double*, int* );
  }
  namespace Pescan { class BandGap; }

  //! \endcond


  namespace Pescan
  {
    struct Dipole
    {
      //! The complex dipole moment.
      std::complex<types::t_real> r[3];
      //! band index i. Negative are conduction bands.
      int i;
      //! band index j. Negative are conduction bands.
      int j;
      //! 0 is up, 1 is down.
      bool ispin;
      //!  0 is up, 1 is down.
      bool jspin;
      //! Serializes a dipole moment.
      template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
        {  _ar & r; _ar & i; _ar & j; _ar & ispin ; _ar & jspin; }
    };
    //! Returns the valence-conduction the norm of the band dipole elements.
    types::t_real dipole_elements( const BandGap& _bandgap, 
                                   types::t_real _degeneracy = types::tolerance );
  }

#endif
