//
//  Version: $Id$
//
#include <iostream>
#include <sstream>
#include <string>

#include <opt/types.h>
#include <opt/debug.h>

#include "functional_builder.h"


namespace CE 
{
   template< class T_HARMONIC > Crystal::Lattice*
     Builder<T_HARMONIC>::lattice     = NULL;
   template< class T_HARMONIC > typename Builder<T_HARMONIC> :: t_Clusters*
     Builder<T_HARMONIC>::clusters    = NULL;
   template< class T_HARMONIC > typename Builder<T_HARMONIC> :: t_CS*
     Builder<T_HARMONIC>::harmonics   = NULL;

   template< class T_HARMONIC >
   Builder<T_HARMONIC> :: ~Builder<T_HARMONIC>()
   {
     if ( lattice )   delete lattice; 
     if ( clusters )  delete clusters; 
     if ( harmonics ) delete harmonics; 
     lattice = NULL; clusters = NULL; harmonics = NULL;
   }


   template< class T_HARMONIC >
   bool Builder<T_HARMONIC> :: Load (const TiXmlElement &_node)
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
     __DOASSERT( not parent,
                 "Could not find an <Functional type=\"CE\"> tag in input file.\n" )

     // creates static quantities
     __TRYCODE(
       lattice    = new Crystal::Lattice;
       clusters   = new t_Clusters;
       harmonics  = new t_CS;,
       " Memory allocation failure in Builder :: Load (...) \n" 
     )

     // then loads lattice
     child = parent->FirstChildElement( "Lattice" );
     if ( not child )
     {
       const TiXmlNode *grandparent = parent->Parent();
       if ( parent )
         child = grandparent->FirstChildElement( "Lattice" );
       __DOASSERT( not child, "Could not find Lattice in input.\n" )
     }
     __DOASSERT( not lattice->Load(*child),
                 "Error while reading Lattice from input.\n" )
     lattice->find_space_group();
     Crystal :: Structure :: lattice = lattice;
     
     // then load clusters
     child = parent->FirstChildElement( "Clusters" );
     __DOASSERT( not child,
                 "Could not find Clusters in input.\n" )
     child = child->FirstChildElement( "Cluster" );
     t_Cluster cluster;
     for (  ; child; child = child->NextSiblingElement("Cluster") )
     {
       __DOASSERT( not cluster.Load( *child ),
                   "Error while loading cluster from input.\n" )
       clusters->push_back(cluster);
     }
     
     // reads in Constituent Strain (coefficients are static)
     child = parent->FirstChildElement( "CS" );
     __DOASSERT( not child,
                 "Could not find CS in input.\n" )
     __DOASSERT( not harmonics->Load_Harmonics( *child ),
                 "Error while loading harmonics from input.\n" )

     return true;
   }

   // finds all clusters, including symmetric equivalents
   // starting from cluster included in Lamarck::clusters
   // results are stored in a new Lamarck::clusters
   template< class T_HARMONIC >
   void  Builder<T_HARMONIC> :: add_equivalent_clusters() 
   {
     const atat::rMatrix3d &cell            = lattice->cell;
     const atat::rMatrix3d inv_cell         = !cell;
     const atat::Array<atat::rMatrix3d> &point_op = lattice->space_group.point_op;
     const atat::Array<atat::rVector3d> &trans    = lattice->space_group.trans;
     // new clusters will be added directly to Lamarck::clusters
     t_Clusters old_cluster_list = *clusters;  
 
     t_Clusters :: iterator i_old_list_last = old_cluster_list.end();
     t_Clusters :: iterator i_old_list = old_cluster_list.begin();
     t_Cluster *transfo_cluster = new t_Cluster;
 
     for (; i_old_list != i_old_list_last ; ++i_old_list )
     {
       for (types::t_int op=0; op<point_op.get_size(); op++)
       {
         // initialize a new cluster to object pointed by i_old_list
         *transfo_cluster = *i_old_list;
         
         // transforms cluster according to symmetry group
         transfo_cluster->apply_symmetry(point_op(op),trans(op));
 
         // checks wether transformed cluster is in Lamarck::clusters
         t_Clusters :: iterator i_cluster  = clusters->begin();
         t_Clusters :: iterator i_last = clusters->end();
         for ( ; i_cluster != i_last ; ++i_cluster)
           if ( transfo_cluster->equivalent_mod_cell(*i_cluster, inv_cell) ) 
             break;
 
         // if it isn't, adds cluster to clusters
         if (i_cluster == i_last ) 
           clusters->push_back( *transfo_cluster );
       }
     }
     delete transfo_cluster;
   }

   // converts the clusters into a polynomial forme
   // converts each cluster for each atom to a monome
   // One function per structure
   // each function is made up of a polynomial and constituent strain
   // part
   template< class T_HARMONIC >
   bool Builder<T_HARMONIC> :: generate_functional( const Crystal::Structure &str,
                                                   t_VA_Functional * const functional )
   {
     std::pair<t_Chemical*, t_CS*> pair( generate_functional( str ) );
     if( (not pair.first) or (not pair.second) ) 
     {
       if( pair.first ) delete pair.first;
       if( pair.second ) delete pair.second;
       return false;
     }
     functional->set_functional2( pair.second );
     functional->set_functional1( pair.first );
 
     return true;
   }

   template< class T_HARMONIC >
   std::pair< typename Builder<T_HARMONIC> :: t_Chemical*,
              typename Builder<T_HARMONIC> :: t_CS* >
     Builder<T_HARMONIC> :: generate_functional( const Crystal::Structure &str )
     {
       t_Chemical *polynome;
  
       // finally, creates polynomials
  
       atat::rMatrix3d inv_cell = !(~str.cell);
       polynome = new t_Chemical();
       Crystal :: Structure :: t_Atoms :: const_iterator i_atom = str.atoms.begin();
       Crystal :: Structure :: t_Atoms :: const_iterator i_atom_last = str.atoms.end();
  
       for(; i_atom != i_atom_last; ++i_atom) // loop over atoms
       {
         t_Clusters :: iterator i_cluster = clusters->begin();
         t_Clusters :: iterator i_cluster_last = clusters->end();
         atat::rVector3d atom_pos = i_atom->pos;
         
         for( ; i_cluster != i_cluster_last; ++i_cluster ) // loop over clusters
         {
           typedef std::vector<atat::rVector3d> :: iterator vec_iterator;
           vec_iterator i_cpos_begin = i_cluster->vectors.begin();
           vec_iterator i_cpos_center = i_cluster->vectors.begin();
           vec_iterator i_cpos_last = i_cluster->vectors.end();
           vec_iterator i_cpos;
           if ( i_cluster->vectors.size() == 0 )
           {
             i_cpos_center = i_cpos_last;
             polynome->add( function::Monome<>(i_cluster->eci) ); 
           }
  
           for (; i_cpos_center != i_cpos_last; ++i_cpos_center ) 
           {   
             // sets up a monome with the correct coefficient
             function::Monome<> monome(   i_cluster->eci
                                        / ( (Real) i_cluster->vectors.size() ) );
  
             i_cpos = i_cpos_begin;
             for (; i_cpos != i_cpos_last; ++i_cpos ) // loop over cluster points
             {
               atat::rVector3d shift = atom_pos - *i_cpos_center;
               
               if (is_int( (!lattice->cell)*shift)) 
               {
                 // finds atom to which "point" is equivalent
                 types::t_unsigned index;
                 for (index = 0; index<str.atoms.size(); ++index)  
                   if ( atat::equivalent_mod_cell(*i_cpos + shift,
                                                  str.atoms[index].pos,inv_cell) ) 
                     break;
                 
                 __ASSERT( index == str.atoms.size(),
                             "Could not find equivalent of atom "
                           << (*i_cpos + shift) << " in Lamarck::generate_functionals\n" )
                 monome.add_term(index, true);
               }
               
             }  // end of loop over cluster points
  
             polynome->add(monome);
             
           } // end of rotation
  
           
         }  // end of loop over clusters
       }  // end of loop over atoms 
  
       (*polynome) *= 1.0 /( (types::t_real) str.atoms.size() );
  
       // now computes constituent strain 
       return std::pair< t_Chemical*, t_CS* >( polynome, new t_CS(str) );
     }
} // namespace CE



