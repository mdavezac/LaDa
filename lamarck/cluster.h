//
//  Version: $Id$
//
#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <tinyxml/tinyxml.h>

#include "opt/types.h"

#include "atat/machdep.h" 
#include "atat/vectmac.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace VA_CE{ class Functional_Builder; };
namespace Ising_CE
{
  class Cluster 
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<Ising_CE::Cluster> ( Ising_CE::Cluster& );
#endif
    friend class VA_CE::Functional_Builder;
    
    protected:
      std::vector<atat::rVector3d> vectors;
      Real eci;

    public:
      Cluster () { eci = 0; };
      Cluster (const Cluster &_cluster){ copy( _cluster ); }
      ~Cluster (){ vectors.clear(); };

      Real set_eci(const Real _eci) { eci = _eci; return eci; };
      void copy(const Cluster &_cluster) { eci = _cluster.eci; vectors = _cluster.vectors; }
      void operator=( const Cluster &_cluster) { copy (_cluster); };
      void apply_symmetry(const atat::rMatrix3d &point_op, const atat::rVector3d &trans);
      bool equivalent_mod_cell( Cluster &equiv, const atat::rMatrix3d &inv_cell);

      const atat::rVector3d& access(types::t_int i) const 
        { return *(vectors.begin()+i); }
      atat::rVector3d& access(types::t_int i)
        { return *(vectors.begin()+i); }
      const atat::rVector3d& operator[](types::t_int i) const 
        { return *(vectors.begin()+i); }
      atat::rVector3d& operator[](types::t_int i) 
        { return *(vectors.begin()+i); }

      types::t_int size() const { return vectors.size(); };

      void print_out( std::ostream &stream ) const; 
      bool Load(const TiXmlElement &_node);
  };
  
} // namespace Ising_CE
#endif
