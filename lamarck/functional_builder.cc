#include <iostream>
#include <sstream>
#include <string>
#include <exception>
#include <stdexcept>

#include <opt/types.h>
#include <analysis/analyze_code.h>
#include <analysis/call_log.h>

#include "functional_builder.h"


namespace VA_CE 
{
   Ising_CE::Lattice*              Functional_Builder::lattice     = NULL;
   std::vector<Ising_CE::Cluster>* Functional_Builder::clusters    = NULL;
   Constituent_Strain*             Functional_Builder::harmonics   = NULL;

   Functional_Builder :: ~Functional_Builder()
   {
     if ( lattice ) 
       delete lattice; 
     if ( clusters ) 
       delete clusters; 
     if ( harmonics ) 
       delete harmonics; 
     lattice = NULL; clusters = NULL;
     harmonics = NULL;
   }


   bool Functional_Builder :: Load (const TiXmlElement &_node)
   {
     const TiXmlElement *child, *parent;
     atat::rVector3d vec;
     
     // This whole section tries to find a <Functional type="vff"> tag
     // in _node or its child
     std::string str = _node.Value();
     if ( str.compare("Functional" ) != 0 )
       parent = _node.FirstChildElement("Functional");
     else
       parent = &_node;
     
     
     while (parent)
     {
       str = "";
       if ( parent->Attribute( "type" )  )
         str = parent->Attribute("type");
       if ( str.compare("CE" ) == 0 )
         break;
       parent = parent->NextSiblingElement("Functional");
     }
     if ( not parent )
     {
       std::cerr << "Could not find an <Functional type=\"CE\"> tag in input file" 
                 << std::endl;
       return false;
     } 

     // creates static quantities
     try
     {
       lattice    = new Ising_CE::Lattice;
       clusters   = new std::vector<Ising_CE::Cluster>;
       harmonics  = new Constituent_Strain;
     }
     catch( std::exception &e )
     {
       std::ostringstream s(" Memory allocation failure in Functional_Builder :: Load (...) ");
       s << std::endl << e.what() << std::endl;
       throw std::runtime_error(s.str());
     }
     
     // then loads lattice
     child = parent->FirstChildElement( "Lattice" );
     if ( not child )
     {
       const TiXmlNode *grandparent = parent->Parent();
       if ( parent )
         child = grandparent->FirstChildElement( "Lattice" );
       if ( not child )
       {
         std::cerr << "Could not find Lattice in input" << std::endl;
         return false;
       }
     }
     if ( not lattice->Load(*child) )
     {
       std::cerr << "Error while reading Lattice from input" << std::endl;
       return false;
     }
     Ising_CE :: Structure :: lattice = lattice;
     
     // then load clusters
     child = parent->FirstChildElement( "Clusters" );
     if (not child )
     {
       std::cerr << "Could not find Clusters in input" << std::endl;
       return false;
     }
     child = child->FirstChildElement( "Cluster" );
     Ising_CE::Cluster cluster;
     for (  ; child; child = child->NextSiblingElement("Cluster") )
     {
       if( !cluster.Load( *child ) )
         return false;
       clusters->push_back(cluster);
     }
     
     // reads in Constituent Strain (coefficients are static)
     child = parent->FirstChildElement( "CS" );
     if ( not child )
     {
       std::cerr << "Could not find CS in input" << std::endl;
       return false;
     }
     if( not harmonics->Load_Harmonics( *child ) )
     {
       std::cerr << "Error while loading harmonics from input" << std::endl;
       return false;
     }

     return true;
   }

   // finds all clusters, including symmetric equivalents
   // starting from cluster included in Lamarck::clusters
   // results are stored in a new Lamarck::clusters
   void  Functional_Builder :: add_equivalent_clusters() 
   {
     START_ANALYSIS("Lamarck :: add_equivalent_clusters");
     const atat::rMatrix3d &cell            = lattice->cell;
     const atat::rMatrix3d inv_cell         = !cell;
     const atat::Array<atat::rMatrix3d> &point_op = lattice->space_group.point_op;
     const atat::Array<atat::rVector3d> &trans    = lattice->space_group.trans;
     // new clusters will be added directly to Lamarck::clusters
     std::vector<Ising_CE::Cluster> old_cluster_list = *clusters;  
 
     std::vector<Ising_CE::Cluster> :: iterator i_old_list_last = old_cluster_list.end();
     std::vector<Ising_CE::Cluster> :: iterator i_old_list = old_cluster_list.begin();
     Ising_CE::Cluster *transfo_cluster = new Ising_CE::Cluster;
 
     for (; i_old_list != i_old_list_last ; ++i_old_list )
     {
       for (types::t_int op=0; op<point_op.get_size(); op++)
       {
         // initialize a new cluster to object pointed by i_old_list
         transfo_cluster->copy( *i_old_list );
         
         // transforms cluster according to symmetry group
         transfo_cluster->apply_symmetry(point_op(op),trans(op));
 
         // checks wether transformed cluster is in Lamarck::clusters
         std :: vector<Ising_CE::Cluster> :: iterator i_cluster  = clusters->begin();
         std :: vector<Ising_CE::Cluster> :: iterator i_last = clusters->end();
         for ( ; i_cluster != i_last ; ++i_cluster)
           if ( transfo_cluster->equivalent_mod_cell(*i_cluster, inv_cell) ) 
             break;
 
         // if it isn't, adds cluster to clusters
         if (i_cluster == i_last ) 
           clusters->push_back( *transfo_cluster );
       }
     }
     delete transfo_cluster;
     END_ANALYSIS;
   }

   // converts the clusters into a polynomial forme
   // converts each cluster for each atom to a monome
   // One function per structure
   // each function is made up of a polynomial and constituent strain
   // part
   bool Functional_Builder :: generate_functional( const Ising_CE::Structure &str,
                                                   t_VA_Functional * const functional )
   {
     Constituent_Strain *strain;
     Polynome *polynome;
 
     // finally, creates polynomials
     START_ANALYSIS( "Functional_Builder :: generate_functionals" )
 
     atat::rMatrix3d inv_cell = !str.cell;
     polynome = new Polynome();
     std::vector<Ising_CE::Atom> :: const_iterator i_atom = str.atoms.begin();
     std::vector<Ising_CE::Atom> :: const_iterator i_atom_last = str.atoms.end();
 
     for(; i_atom != i_atom_last; ++i_atom) // loop over atoms
     {
       std::vector<Ising_CE::Cluster> :: iterator i_cluster = clusters->begin();
       std::vector<Ising_CE::Cluster> :: iterator i_cluster_last = clusters->end();
       atat::rVector3d atom_pos = i_atom->pos;
       
       for( ; i_cluster != i_cluster_last; ++i_cluster ) // loop over clusters
       {
         std::vector<atat::rVector3d> :: iterator i_cpos_begin = i_cluster->vectors.begin();
         std::vector<atat::rVector3d> :: iterator i_cpos_center = i_cluster->vectors.begin();
         std::vector<atat::rVector3d> :: iterator i_cpos_last = i_cluster->vectors.end();
         std::vector<atat::rVector3d> :: iterator i_cpos;
         if ( i_cluster->vectors.size() == 0 )
         {
           i_cpos_center = i_cpos_last;
           polynome->add( function::Monome<>(i_cluster->eci) ); 
         }
 
         for (; i_cpos_center != i_cpos_last; ++i_cpos_center ) 
         {   
           // sets up a monome with the correct coefficient
           function::Monome<> monome( i_cluster->eci / ( (Real) i_cluster->vectors.size() ) );
 
           i_cpos = i_cpos_begin;
           for (; i_cpos != i_cpos_last; ++i_cpos ) // loop over cluster points
           {
             atat::rVector3d shift = atom_pos - *i_cpos_center;
             
             if (is_int( (!lattice->cell)*shift)) 
             {
               // finds atom to which "point" is equivalent
               types::t_unsigned index;
               for (index = 0; index<str.atoms.size(); ++index)  
                 if ( atat::equivalent_mod_cell(*i_cpos + shift, str.atoms[index].pos,inv_cell) ) 
                   break;
               
               #ifdef _DEBUG_LADA_
                 if ( index == str.atoms.size() )
                 {
                   std::cerr << "Could not find equivalent atom"
                             << " of cluster in Lamarck::generate_functionals" 
                             << std::endl;
                   exit(0);
                 }
               #endif 
               monome.add_term(index, true);
             }
             
           }  // end of loop over cluster points
 
           polynome->add(monome);
           
         } // end of rotation
 
         
       }  // end of loop over clusters
     }  // end of loop over atoms 
 
     (*polynome) *= 1.0 /( (types::t_real) str.atoms.size() );
 
     // now computes constituent strain 
     strain = new Constituent_Strain(str);
     
     // finally creates FUNCTION object and adds it to the list
     functional->set_functional2( strain );
     functional->set_functional1( polynome );
 
     END_ANALYSIS;
 
     return true;
   }




} // namespace VA_CE




#ifdef _MPI

namespace mpi
{
  template<>
  bool BroadCast :: serialize<VA_CE::Functional_Builder>( VA_CE::Functional_Builder &_b )
  {
    if ( not _b.lattice ) _b.lattice = new Ising_CE::Lattice;
    if ( not _b.harmonics ) _b.harmonics = new Ising_CE::Constituent_Strain;
    if ( not _b.clusters ) _b.clusters = new std::vector<Ising_CE::Cluster>;

    if ( not serialize( *_b.lattice ) ) return false;
    if ( not serialize( *_b.harmonics ) ) return false;
    
    // then serializes rvecs and kvecs
    types::t_int n = _b.clusters->size();
    if( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      _b.clusters->resize(n);
    std::vector<Ising_CE::Cluster> :: iterator i_cluster = _b.clusters->begin();
    std::vector<Ising_CE::Cluster> :: iterator i_cluster_end = _b.clusters->end();
    for(; i_cluster != i_cluster_end; ++i_cluster )
      if ( not serialize( *i_cluster ) ) return false;

    return true;
  }
}

#endif
