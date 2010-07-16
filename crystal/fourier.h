#ifndef _LADA_CRYSTAL_FOURIER_H_
#define _LADA_CRYSTAL_FOURIER_H_

#include "LaDaConfig.h"

#include <math.h>

#include <opt/types.h>

namespace LaDa 
{
  namespace Crystal {

    //! \brief Defines a %Fourier transform for a structure within a single-site lattice.
    //! \warning There is no reason this should work well in multiple-site lattices....
    struct Fourier
    {
      //! \brief From real to k space
      //! \details The first range is most likely some instance of
      //!          Crystal::Structure::t_Atoms. The second range is similarly an
      //!          instance of Crystal::Structure::t_kAtoms. The third range
      //! \param[in, out] _rfirst iterator to the first real space atom (of a type
      //!                similar to Crystal::Atom_Type)
      //! \param[in, out] _rend iterator to the last real space atom (of a type
      //!              similar to Crystal::Atom_Type)
      //! \param[in] _kfirst iterator to the first real space atom (of a type
      //!                similar to Crystal::Atom_Type< std::complex >)
      //! \param[in] _kend iterator to the last real space atom (of a type
      //!              similar to Crystal::Atom_Type< std::complex >)
      template<class T_R_IT, class T_K_IT>
      Fourier( T_R_IT _rfirst, T_R_IT _rend,
               T_K_IT _kfirst, T_K_IT _kend );
      //! \brief From k space to real space. 
      //! \details The first range is most likely some instance of
      //!          Crystal::Structure::t_Atoms. The second range is similarly an
      //!          instance of Crystal::Structure::t_kAtoms. The third range
      //!          should be iterators to std::complex.
      //! \pre The range [ \a _rout, \a _rout += \a _rfirst - \a _rend  ) should
      //!      be valid.
      //! \param[in] _rfirst iterator to the first real space atom (of a type
      //!                similar to Crystal::Atom_Type)
      //! \param[in] _rend iterator to the last real space atom (of a type
      //!              similar to Crystal::Atom_Type)
      //! \param[in] _kfirst iterator to the first real space atom (of a type
      //!                similar to Crystal::Atom_Type< std::complex >)
      //! \param[in] _kend iterator to the last real space atom (of a type
      //!              similar to Crystal::Atom_Type< std::complex >)
      //! \param[out] _rout iterator to the first complex real-space
      //!              occupation ( of std::complex type )
      template<class T_R_IT, class T_K_IT, class T_O_IT >
      Fourier( T_R_IT _rfirst, T_R_IT _rend,
               T_K_IT _kfirst, T_K_IT _kend,
               T_O_IT _rout ); // sets rvector values from kspace values
    };

    template<class T_R_IT, class T_K_IT>
    Fourier :: Fourier( T_R_IT _rfirst, T_R_IT _rend,
                        T_K_IT _kfirst, T_K_IT _kend )
    {
      const types::t_complex
         imath(0, 2*3.1415926535897932384626433832795028841971693993751058208);
      
      for (; _kfirst != _kend; ++_kfirst)
      {
        _kfirst->type = std::complex<types::t_real>(0);
        for(T_R_IT i_r( _rfirst ); i_r != _rend; ++i_r )
        {
          _kfirst->type +=   exp( imath * ( i_r->pos[0] * _kfirst->pos[0] +
                                            i_r->pos[1] * _kfirst->pos[1] +
                                            i_r->pos[2] * _kfirst->pos[2] ) )
                           * i_r->type;
        }
      }
    }
    template<class T_R_IT, class T_K_IT, class T_O_IT >
    Fourier :: Fourier( T_R_IT _rfirst, T_R_IT _rend,
                        T_K_IT _kfirst, T_K_IT _kend,
                        T_O_IT _rout ) // sets rvector values from kspace values
    {
      const types::t_complex
         imath(0, 2*3.1415926535897932384626433832795028841971693993751058208);
      for (; _rfirst != _rend; ++_rfirst, ++_rout)
      {
        *_rout = types::t_complex(0,0);
        for(T_K_IT i_k=_kfirst; i_k != _kend; ++i_k)
        {
          *_rout +=   exp( imath * ( _rfirst->pos[0] * i_k->pos[0] +
                                     _rfirst->pos[1] * i_k->pos[1] +
                                     _rfirst->pos[2] * i_k->pos[2] ) )
                    * i_k->type;
        }
      }
    }
  } // namespace Crystal

} // namespace LaDa

#endif
