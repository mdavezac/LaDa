//
//  Version: $Id$
//
#ifndef _LADA_GA_PURELAYERS_FOURIER_H_
#define _LADA_GA_PURELAYERS_FOURIER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <crystal/structure.h>
#include <opt/types.h>

#include "concentration.h"

namespace LaDa
{
  namespace GA
  {
    namespace PureLayers
    {
      //! \brief Defines fourier transforms for two-site lattices
      //! \details The fourier transform assign +/-1 to the
      //!          first lattice site and +/- i to the second lattice site.
      struct Fourier
      {
        //! \brief From real to k space
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
        const std::complex<types::t_real>
           imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
        
        for (; _kfirst != _kend; ++_kfirst)
        {
          _kfirst->type = std::complex<types::t_real>(0);
          for(T_R_IT i_r( _rfirst ); i_r != _rend; ++i_r )
          {
            // used to be directly in std::complex constructor below...
            // but then it hit me: which constructor argument does c++ look at first?
            // this is safe now.
            types::t_real a = i_r->type;
            types::t_real b = (++i_r)->type;
            _kfirst->type +=   exp( imath * ( i_r->pos[0] * _kfirst->pos[0] +
                                              i_r->pos[1] * _kfirst->pos[1] +
                                              i_r->pos[2] * _kfirst->pos[2] ) )
                             * std::complex<types::t_real>(a, b);
          }
        }
      }
      template<class T_R_IT, class T_K_IT, class T_O_IT >
      Fourier :: Fourier( T_R_IT _rfirst, T_R_IT _rend,
                          T_K_IT _kfirst, T_K_IT _kend,
                          T_O_IT _rout ) // sets rvector values from kspace values
      {
        const std::complex<types::t_real>
           imath(0, 2*3.1415926535897932384626433832795028841971693993751058208);
        for (; _rfirst != _rend; _rfirst+=2, ++_rout)
        {
          *_rout = 0.0;
          for(T_K_IT i_k=_kfirst; i_k != _kend; ++i_k)
          {
            *_rout +=   exp( imath * ( _rfirst->pos[0] * i_k->pos[0] +
                                       _rfirst->pos[1] * i_k->pos[1] +
                                       _rfirst->pos[2] * i_k->pos[2] ) )
                      * i_k->type;
          }
        }
      }
    } // namespace PureLayers
  } // namespace GA
} // namespace LaDa


#endif // _TWOSITES_OBJECT_H_
