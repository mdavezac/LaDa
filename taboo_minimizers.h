#ifndef _TABOO_MINIMIZERS_H_
#define _TABOO_MINIMIZERS_H_

#include <algorithm>
#include <vector>

#include <eo/eoOp.h>

#include "taboo.h"

#include <opt/types.h>
using namespace types;

namespace LaDa
{
  template <class t_Object, class t_Call_Back> 
  class SA_TabooOp : public eoMonOp<t_Object>
  {
    protected:
      t_Call_Back &call_back;
      Taboo_Base< t_Object > &taboo;
      std::vector<t_unsigned> directions;
      t_unsigned max_directions_checked;
      static t_unsigned nb_evals;

    public:
      SA_TabooOp   ( t_Call_Back &_call_back,
                     Taboo_Base<t_Object> &_taboo,
                     t_unsigned &_max )
                 : call_back( _call_back ), taboo(_taboo),
                   max_directions_checked(_max) {}
      SA_TabooOp   ( const SA_TabooOp<t_Object, t_Call_Back> &_sa )
                 : call_back( _sa.call_back ), taboo(_sa.taboo),
                   max_directions_checked(_sa.max_directions_checked) {}
      virtual ~SA_TabooOp() {}

      virtual std::string className() const { return "LaDa::SA_TabooOp"; }

      // minimizes a function obtained from call_back, starting from
      // position _object. The result is stored in _object
      virtual bool operator() ( t_Object &_object )
      {
        // first gets the function to minimize
        typename t_Call_Back :: t_Functional &functional = call_back.get_functional( _object );
        functional.set_variables( _object.get_variables() );

        // then gets problem size of object
        t_unsigned size = _object.size();

        // creates a vector of ordered directions
        directions.clear();
        { // implements sgl iota
          // -- i.e. fills which_direction with increasing intergers
          for( t_unsigned i =0; i < size; ++i)
            directions.push_back(i);
        }

        // now for the simmulated annealing part
        t_unsigned count = 0; // number of calls
        std::vector<t_unsigned> :: iterator i_begin = directions.begin();
        std::vector<t_unsigned> :: iterator i_last = directions.end();
        std::vector<t_unsigned> :: iterator i_dir;
        bool is_not_converged;
        typename t_Object::iterator i_var_begin = functional.begin();
        typename t_Object::iterator i_var;
        t_real current_e, next_e;

        do 
        {
          // shuffle directions
          std::random_shuffle(i_begin, i_last);
          is_not_converged = false;

          // evaluates first position
          current_e = functional.evaluate();
          ++nb_evals;

          for( i_dir = i_begin; i_dir != i_last; ++i_dir )
          {
            ++count; // one more flip 
            i_var = i_var_begin + *i_dir; // sets direction
            *i_var = ( *i_var > 0 ) ? -1.0 : 1.0; // flips spins

            // we only want to flip if structure not in taboo lists
            if ( taboo(_object) )  
              *i_var = ( *i_var > 0 ) ? -1.0 : 1.0; // flips spin back
            else
            {
              next_e = functional.evaluate(); // evaluates
               ++nb_evals;
              if ( current_e > next_e ) // we do wanna flip!
              {
                is_not_converged = true;
                current_e = next_e;
              }
              else
                *i_var = ( *i_var > 0 ) ? -1.0 : 1.0; // flips spin back
            }

            if ( count >= max_directions_checked )  // breaks out after n calls to functional
              { is_not_converged = false; break; }

          } // loop over shuffled directions

        } while ( is_not_converged );

        // we have values, so lets use them
        next_e = call_back.evaluate( _object.get_concentration() );
        _object.set_quantity( current_e + next_e );
        _object.set_baseline( next_e );
        _object.set_fitness();
        // and add to Convex_Hull if necessary
        if ( current_e < 0 )
          call_back.add_to_convex_hull( _object );
        return false; // fitness is set, returns false
      } // end of operator()(t_Object &_object)

  };
  
  template <class t_Object, class t_Call_Back> 
  t_unsigned SA_TabooOp<t_Object, t_Call_Back> :: nb_evals = 0;

  template <class t_Object, class t_Call_Back> 
  class GradientSA_TabooOp : public eoMonOp<t_Object>
  {
    protected:
      t_Call_Back &call_back;
      Taboo_Base< t_Object > &taboo;
      std::vector<t_unsigned> directions;
      t_unsigned max_directions_checked;
      static t_unsigned nb_evals, nb_grad_evals;

    public:
      GradientSA_TabooOp   ( t_Call_Back &_call_back,
                             Taboo_Base<t_Object> &_taboo,
                             t_unsigned &_max )
                         : call_back( _call_back ), taboo(_taboo),
                           max_directions_checked(_max) {}
      GradientSA_TabooOp   ( const GradientSA_TabooOp<t_Object, t_Call_Back> &_sa )
                         : call_back( _sa.call_back ), taboo(_sa.taboo),
                           max_directions_checked(_sa.max_directions_checked) {}

      virtual ~GradientSA_TabooOp() {}

      virtual std::string className() const { return "LaDa::GradientSA_TabooOp"; }

      // minimizes a function obtained from call_back, starting from
      // position _object. The result is stored in _object
      virtual bool operator() ( t_Object &_object )
      {
        // first gets the function to minimize
        typename t_Call_Back :: t_Functional &functional = call_back.get_functional( _object );
        functional.set_variables( _object.get_variables() );

        // then gets problem size of object
        t_unsigned size = _object.size();

        // creates a vector of ordered directions
        directions.clear();
        { // implements sgl iota
          // -- i.e. fills which_direction with increasing intergers
          for( t_unsigned i =0; i < size; ++i)
            directions.push_back(i);
        }

        // now for the simmulated annealing part
        t_unsigned count = 0; // number of calls
        std::vector<t_unsigned> :: iterator i_begin = directions.begin();
        std::vector<t_unsigned> :: iterator i_last = directions.end();
        std::vector<t_unsigned> :: iterator i_dir;
        bool is_not_converged;
        typename t_Object::iterator i_var_begin = functional.begin();
        typename t_Object::iterator i_var;
        t_real current_e, next_e;

        do 
        {
          // shuffle directions
          std::random_shuffle(i_begin, i_last);
          is_not_converged = false;

          // evaluates first position
          current_e = functional.evaluate();
          ++nb_evals;

          for( i_dir = i_begin; i_dir != i_last; ++i_dir )
          {
            ++count; // one more flip 
            i_var = i_var_begin + *i_dir; // sets direction
            *i_var = ( *i_var > 0 ) ? -1.0 : 1.0; // flips spins

            // we only want to flip if structure not in taboo lists
            if ( taboo(_object) ) 
              *i_var = ( *i_var > 0 ) ? -1.0 : 1.0; // flips back
            else 
            {
              *i_var = ( *i_var > 0 ) ? -0.998 : 0.998; // flips to derivative position
              next_e = functional.evaluate(); // evaluates
              ++nb_grad_evals;
              if ( current_e > next_e ) // gradient says to flip
              {
                *i_var = ( *i_var > 0.0 ) ? -1.0 : 1.0; // so we flip
                next_e = functional.evaluate();
                ++nb_evals;
                if ( current_e > next_e )
                { // yeah! gradient was right
                  is_not_converged = true;
                  current_e = next_e;
                }
                else // flips back
                  *i_var = ( *i_var > 0.0 ) ? -1.0 : 1.0; 
              }
            } // end of check for taboo

            if ( count >= max_directions_checked )  // breaks out after n calls to functional
              { is_not_converged = false; break; }

          } // loop over shuffled directions

        } while ( is_not_converged );

        // we have values, so lets use them
        next_e = call_back.evaluate( _object.get_concentration() );
        _object.set_quantity( current_e + next_e );
        _object.set_baseline( next_e );
        _object.set_fitness();
        // and add to Convex_Hull if necessary
        if ( current_e < 0 )
          call_back.add_to_convex_hull( _object );
        return false; // fitness is set, returns false
      } // end of operator()(t_Object &_object)

  };

  template <class t_Object, class t_Call_Back> 
  t_unsigned GradientSA_TabooOp<t_Object, t_Call_Back> :: nb_evals = 0;
  template <class t_Object, class t_Call_Back> 
  t_unsigned GradientSA_TabooOp<t_Object, t_Call_Back> :: nb_grad_evals = 0;

} // namespace LaDa
#endif
