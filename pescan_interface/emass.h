//
//  Version: $Id$
//
#ifndef _PESCAN_EMASS_H_
#define _PESCAN_EMASS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <utility>
#include <string>

#include <tinyxml/tinyxml.h>

#include <minimizer/cgs.h>
#include <crystal/structure.h>
#include <opt/types.h>

#include "interface.h"

namespace LaDa
{
  namespace Pescan 
  {
    class eMass
    {
      public:
        //! \brief Type of the ouput argument.
        //! \details First of pair is the eigenvalue, second of pair is the effective mass.
        typedef std::vector< std::pair<types::t_real, types::t_real> > t_Output;
        //! Conjugate gradient minimizer.
        typedef Fitting :: Cgs t_Minimizer;
        //! Conjugate gradient minimizer.
        t_Minimizer cgs;
        //! Order of interpolation.
        size_t order;
        //! Number of interpolation points.
        size_t npoints;
        //! Stepsize for interpolation points.
        types::t_real stepsize;
        //! k-vector for which to compute effective mass.
        Eigen::Vector3d kpoint;
        //! Direction for which to compute effective mass.
        Eigen::Vector3d direction;
        //! Number of states for which to compute effective mass.
        size_t nbstates;

        //! Constructor.
        eMass() : order( 2 ), npoints(2), stepsize( 1e-2 ),
                  kpoint( 0,0,0 ), direction(1,1,1), nbstates( 2 ) {}
        //! Copy Constructor.
        eMass   ( const eMass &_c )
              : cgs( _c.cgs ), order( _c.order ), 
                npoints( _c.npoints ), stepsize( _c.stepsize ),
                kpoint( _c.kpoint ), direction( _c.direction ),
                nbstates( _c.nbstates ) {}

        //! \brief Computes effective mass.
        //! \param[in] _interface the escan functional.
        //! \param[in] _ocell cell vectors of the structure prior to relaxation.
        //! \param[in] _structure \e relaxed structure for which to compute emasses.
        //! \param[in] _eref reference energy at which to compute bands.
        //! \param[out] _out resulting effective masses and band eigenvalues.
        void operator()
        (
          const Interface& _interface,
          const Eigen::Matrix3d &_ocell, 
          const Crystal::Structure &_structure,
          const types::t_real &_eref,
          t_Output &_out 
        ) const;

        //! Loads parameters from XML node.
        bool load( const TiXmlElement& _node, const std::string &_name = "emass" );

      protected:
        //! Computes kpoint in distorted structure from original lattice point.
        Eigen::Vector3d distort_kpoint( const Eigen::Matrix3d &_ocell, 
                                        const Eigen::Matrix3d &_dcell,
                                        const Eigen::Vector3d &_kpoint ) const
          { return ( !(~_dcell) ) * (~_ocell) * _kpoint; }
        //! Computes eigenvalues for single kpoint.
        void compute_( Interface& _interface,
                       const Eigen::Vector3d &_kpoints ) const;

        //! Adds a weight to interpolation points.
        types::t_real weight_( size_t _i, size_t j ) const;

    };

  } // namespace Pescan
} // namespace LaDa

#endif // _PESCAN_EMASS_H_
