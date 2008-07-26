//
//  Version: $Id: functional_builder.impl.h 685 2008-07-22 17:20:39Z davezac $
//

#ifndef _CE_REGULARIZATION_H_
#define _CE_REGULARIZATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include <opt/debug.h>
#include <crystal/structures.h>

namespace CE
{
  //! \brief Computes pis of \a _str for \a _clusters.
  //! \param[in] _cluster is a vector of containers of symmetrically equivalent
  //!                     cluster, centered on the origin.
  //! \param[in] _str the structure for which to compute the pis.
  //! \param[out] _pis the computed pis, one per class of symmetrically
  //!                  equivalent clusters.
  template< class T_CLUSTERS, class T_PIS >
  void find_pis( const T_CLUSTERS &_clusters,
                 const Crystal::Structure & _str,
                 T_PIS &_pis );
  //! \brief Computes pis of \a _str for \a _clusters.
  //! \see[in] _cluster is a vector of containers of symmetrically equivalent
  //!                     cluster, centered on the origin.
  //! \param[in] _str structures for which to compute the pis.
  //! \param[out] _pis the computed pis, one per class of symmetrically
  //!                  equivalent clusters.
  template< class T_CLUSTERS, class T_PIS >
  void find_pis( const T_CLUSTERS &_clusters,
                 const std::vector< Crystal::Structure > & _str,
                 T_PIS &_pis );

  //! \brief reads lda energies and structures from NREL input files.
  //! \param[in] _path is the full or relative path to the "LDAs.dat" file.
  //!                  Structure files are expected to be found in the same
  //!                  directory as the LDAs.dat and is deduced from \a _path.
  //! \param[inout] _structures are added to this container.
  template< template<class> T_CONTAINER >
  void read_ce_structures( const std::string &_path,
                           T_CONTAINER<Crystal::Structure> &_structures );
                           

  //! \brief Regulated Cluster-Expansion.
  //! \see <A HREF="http://dx.doi.org/10.1103/PhysRevB.73.224207"> Ralf Drautz
  //!      and Alejandro Diaz-Ortiz, PRB \bf 73, 224207 (2007)</A>.
  class Regulated
  {
    public:
      //! Type of the class representing a cluster.
      typedef CE:Cluster t_Cluster;
      //! Container of equivalent clusters.
      typedef std::vector< t_Cluster > t_EquivClusters;
      //! Container of classes of equivalent xlusters.
      typedef std::vector< t_EquivCluster > t_Clusters;
      //! A container of structures.
      typedef std::vector< Crystal::Structure > t_Structures;
      //! A container of Pis for a single structure.
      typedef std::vector< types::t_int > t_StructPis;
      //! A container of Pis for a single structure.
      typedef std::vector< t_StructPis > t_Pis;
      //! A container of weights.
      typedef std::vector< types::t_real > t_Weights;
      //! Type of the Jacobian.
      typedef std::vector<types::t_real> t_Jacobian;
      //! Type of the return.
      typedef std::vector<types::t_real> t_Return;
      //! Type of the input variables.
      typedef std::vector< t_Return > t_Arg;


    public:
      //! Constructor.
      Regulated() {};
      //! Destructor.
      ~Regulated() {};

      //! Evaluates the cv score for the weights on input.
      void operator( const t_Arg &_arg, t_Return & _return );
      //! Evaluates the jacobian.
      void jacobian( const t_Arg &_arg, t_Jacobian & _jacobian );

    protected:
      //! Initializes Regulated::sums.
      void init_sums();

    protected:
      //! A container of symmetrically equivalent clusters.
      t_Clusters clusters;
      //! A container of pis for all structures.
      t_Pis pis;
      //! Type of Regulated::numsums.
      typedef std::vector< types::t_real > t_NumSums;
      //! Type of Regulated::denumsums.
      typedef std::vector< std::vector< types::t_real > > t_DenumSums;
      //! \f$=\sum_s w_s \phi_{\alpha, s}\phi_{\beta, s}\f$.
      t_DenumSums denumsums;
      //! \f$=\sum_s w_s E_s\phi_{\beta, s}\f$.
      t_NumSums numsums;
      //! A container of weights.
      t_Weights weights;
  }


} // end of namespace CE

#endif 
