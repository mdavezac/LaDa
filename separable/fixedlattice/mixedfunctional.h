//
//  Version: $Id$
//
#ifndef _CE_MIXED_SEPARABLE_H_
#define _CE_MIXED_SEPARABLE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>
#include<string>

#include <opt/types.h>
#include <opt/debug.h>
#include <ce/clusters.h>
#include "functional.h"

namespace CE
{
  //! \cond 
  template< class T_TRAITS, class T_HARMONICS > class MixedSeparables;
  //! \endcond

  //! Computes predictions for all structures in PI file.
  template< class T_TRAITS, class T_HARMONICS >
    void evaluate_pifile( const std::string &_file,
                          MixedSeparables< T_TRAITS, T_HARMOMICS > &_harmonics );

  //! \brief A mixed CE/Sum of sep. functional.
  //! \details Because of its mixed status, the functional takes a complete
  //!          structure on input. 
  template< class T_TRAITS, class T_HARMONICS >
    class MixedSeparables :: protected Separable< T_TRAITS >,
                             protected PosToConfs,
                             protected ConstituentStrain::Functional< T_HARMONICS>
    {
      template<class TT_TRAITS> friend class MixedSeparables;

      public:
        //! Traits of the separable function.
        typedef T_TRAITS t_Traits;
        //! Type of the separable base class.
        typedef Separable<T_TRAITS> t_SepBase;
        //! Type of the position base class.
        typedef PosToConf t_PosBase;
        //! Type of the position base class.
        typedef ConstituentStrain :: Functional<T_HARMONIC> t_CSBase;

        //! Copy constructor.
        MixedSeparables () {}
        //! Copy constructor.
        template<class TT_TRAITS >
          MixedSeparables   ( const MixedSeparables<TT_TRAITS> &_c ) 
                          : t_SepBase( _c ), t_PosBase( _c ),
                            positions_(_c.positions_),
                            clusterclasses_(_c.clusterclasses_) {}

        //! Copies the separable part only.
        template< class TT_TRAITS > 
          void operator=( const Separables<TT_TRAITS > &_c );
        //! Copies the PosToConf part only.
        void operator=( const t_PosBase &_c )
          {  t_PosBase :: syms = _c.syms; }

        //! Returns a reference to the cluster classes.
        t_ClusterClasses& clusterclasses() { return clusterclasses_; }
        //! Returns a constant reference to the cluster classes.
        const t_ClusterClasses& clusterclasses() const { return clusterclasses_; }
        //! Returns a reference to the cluster interaction energies.
        t_Ecis& ecis() { return ecis_; }
        //! Returns a constant reference to the cluster interaction energies.
        const t_Ecis& ecis() const { return ecis_; }

        //! Returns the energy.
        typename t_SepBase::t_Matrix::value_type
          operator()( const Crystal::Structure &_str );

        //! Returns the number of degrees of liberty.
        size_t dof() const { return ecis().size() + t_SepBase::dof(); }

        //! Initializes from input parameters.
        void init( const std::string &_csxml,
                   const std::string &_desc, 
                   const types::t_unsigned _nbpairs = 0,
                   const std::string &_jtypes = "",
                   bool _rmpairs = false,
                   bool _addJO = false,
                   bool _addJ1 = false );


      protected:
        //! Type of the positions basis for the separable function.
        typedef PosToConfs::t_Positions t_Positions;
        //! Type of the clusters.
        typedef Cluster t_Cluster;
        //! Type of an equivalent class of clusters.
        typedef std::vector< t_Cluster > t_ClusterClass;
        //! Type of a all classes of clusters.
        typedef std::list< t_ClusterClass > t_ClusterClasses;
        //! Type of the ecis.
        typedef typename t_SepBase::t_Vector t_Ecis;

        //! All cluster classes.
        t_ClusterClasses clusterclasses_;
        //! The cluster interaction energies;
        t_Ecis ecis_;
    };

  //! Copies norms, coefs, and ecis.
  template<class T_SEPTRAITS, class T_COLTRAITS > 
  void operator=( MixedSeparables<T_SEPTRAITS> &_s, 
                  const MixedApproach<T_COLTRAITS> &_col );

  //! Prints out the separable to a stream.
  template<class T_TRAITS >
  std::ostream& operator<<( std::ostream& _stream,
                            const MixedSeparables<T_TRAITS> &_sep );

}

#include "functional.impl.h"

#endif
