//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include "ce.hpp"
#include "clusters.hpp"
#include "create_pairs.hpp"
#include "create_clusters.hpp"
#include "find_pis.hpp"
#include "mlcluster.hpp"
#include "mlclusters.hpp"

BOOST_PYTHON_MODULE(ce)
{
  LaDa::Python::expose_ce();
  LaDa::Python::expose_clusters();
  LaDa::Python::expose_create_pairs();
  LaDa::Python::expose_create_clusters();
  LaDa::Python::expose_find_pis();
  LaDa::Python::expose_mlcluster();
}
