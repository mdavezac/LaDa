#ifndef _TABOO_MINIMIZERS_H_
#define _TABOO_MINIMIZERS_H_

#include <algorithm>
#include <vector>

#include <eo/eoOp.h>
#include <eo/utils/eoRNG.h>

#include "taboo.h"

#include "evaluation.h"
#include <opt/types.h>
using namespace types;

namespace LaDa
{
  template <class t_Call_Back, class t_Object = typename t_Call_Back :: t_Object> 
  class SA_TabooOp : public eoMonOp<t_Object>
  {
    protected:
      t_Call_Back &call_back;
      Taboo_Base< t_Object > *taboo;
      std::vector<t_unsigned> directions;
      t_unsigned max_directions_checked;
      Taboo< t_Object, std::list<t_Object> > *path_taboo;
      Evaluation< t_Call_Back, t_Object > &eval;
      
    public:
      SA_TabooOp   ( t_Call_Back &_call_back,
                     t_unsigned &_max,
                     Evaluation< t_Call_Back, t_Object > &_eval,
                     Taboo_Base<t_Object> *_taboo = NULL,
                     Taboo<t_Object, std::list<t_Object> > *_pt = NULL )
                 : call_back( _call_back ), taboo(_taboo),
                   max_directions_checked(_max), path_taboo(_pt), eval(_eval) {}
      SA_TabooOp   ( const SA_TabooOp<t_Object, t_Call_Back> &_sa )
                 : call_back( _sa.call_back ), taboo(_sa.taboo),
                   max_directions_checked(_sa.max_directions_checked),
                   path_taboo(_sa.path_taboo), eval(_sa.eval) {}
      virtual ~SA_TabooOp() {}

      virtual std::string className() const { return "LaDa::SA_TabooOp"; }

      // minimizes a function obtained from call_back, starting from
      // position _object. The result is stored in _object
      virtual bool operator() ( t_Object &_object )
      {
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
        t_Object next_o(_object);
        typename t_Object::iterator i_var_begin = next_o.begin();
        typename t_Object::iterator i_var;
        t_real current_e, next_e;

        // evaluates first position
        current_e = eval.evaluate(_object);
        do 
        {
          // shuffle directions
          std::random_shuffle(i_begin, i_last, eo::random<t_unsigned>);
          is_not_converged = false;

          for( i_dir = i_begin; i_dir != i_last; ++i_dir )
          {
            ++count; // one more flip 
            i_var = i_var_begin + *i_dir; // sets direction
            *i_var = ( *i_var > 0 ) ? -1.0 : 1.0; // flips spins

            // we only want to flip if structure not in taboo lists
            if ( taboo and (*taboo)(next_o) )  
              *i_var = ( *i_var > 0 ) ? -1.0 : 1.0; // flips spin back
            else
            {
              next_o.invalidate();
              next_e = eval.evaluate(next_o); // evaluates
              if ( path_taboo ) // add to path_taboo
                path_taboo->add(next_o, false);
              if ( current_e > next_e ) // we do wanna flip!
              {
                is_not_converged = true;
                current_e = next_e;
                _object = next_o;
              }
              else
                *i_var = ( *i_var > 0 ) ? -1.0 : 1.0; // flips spin back
            }

            if ( count >= max_directions_checked )  // breaks out after n calls to functional
              { is_not_converged = false; break; }

          } // loop over shuffled directions

        } while ( is_not_converged );

        return false; // fitness is set, returns false
      } // end of operator()(t_Object &_object)

  };
  
  template <class t_Call_Back, class t_Object = typename t_Call_Back::t_Object> 
  class GradientSA_TabooOp : public eoMonOp<t_Object>
  {
    protected:
      t_Call_Back &call_back;
      Taboo_Base< t_Object > *taboo;
      std::vector<t_unsigned> directions;
      t_unsigned max_directions_checked;
      Taboo<t_Object, std::list<t_Object> > *path_taboo;
      Evaluation<t_Call_Back, t_Object > &eval;

    public:
      GradientSA_TabooOp   ( t_Call_Back &_call_back,
                             t_unsigned &_max,
                             Evaluation<t_Call_Back, t_Object > &_eval,
                             Taboo_Base<t_Object> *_taboo = NULL,
                             Taboo<t_Object, std::list<t_Object> > *_pt = NULL )
                         : call_back( _call_back ), taboo(_taboo),
                           max_directions_checked(_max), path_taboo(_pt), eval(_eval) {}
      GradientSA_TabooOp   ( const GradientSA_TabooOp<t_Object, t_Call_Back> &_sa )
                         : call_back( _sa.call_back ), taboo(_sa.taboo),
                           max_directions_checked(_sa.max_directions_checked),
                           path_taboo( _sa.path_taboo ), eval(_sa.eval) {}

      virtual ~GradientSA_TabooOp() {}

      virtual std::string className() const { return "LaDa::GradientSA_TabooOp"; }

      // minimizes a function obtained from call_back, starting from
      // position _object. The result is stored in _object
      virtual bool operator() ( t_Object &_object )
      {
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
        t_Object next_o = _object;
        t_Object copy = _object;
        typename t_Object::iterator i_var_begin = next_o.begin();
        typename t_Object::iterator i_var;
        t_real current_e, next_e;
  
        // evaluates first position
        current_e = eval.evaluate(_object);

        do 
        {
          // shuffle directions
          std::random_shuffle(i_begin, i_last, eo::random<t_unsigned>);
          is_not_converged = false;

          for( i_dir = i_begin; i_dir != i_last; ++i_dir )
          {
            ++count; // one more flip 
            i_var = i_var_begin + *i_dir; // sets direction
            *i_var = ( *i_var > 0 ) ? -1.0 : 1.0; // flips spins

            // we only want to flip if structure not in taboo lists
            if ( taboo and (*taboo)(next_o) ) 
              *i_var = ( *i_var > 0 ) ? -1.0 : 1.0; // flips back
            else if ( eval.is_known( next_o ) )
            {
              next_e = next_o.value();
              if ( current_e > next_e )
              {
                _object = next_o;
                current_e = next_e; 
                is_not_converged = true;
              }
              else
                *i_var = ( *i_var > 0 ) ? -1.0 : 1.0; // flips back
            }
            else
            {
              // flips to derivative position 
              // ** we've flipped one already for taboo, 
              // ** hence negative sign
              *i_var = ( *i_var > 0 ) ? -0.998 : 0.998; 
              next_e = eval.fast_eval(next_o); // evaluates fast!! no history, etc..
              if ( current_e > next_e ) // gradient says to flip
              {
                *i_var = ( *i_var > 0.0 ) ? -1.0 : 1.0; // so we flip
                next_o.invalidate();
                next_e = eval.evaluate(next_o);
                current_e = eval.evaluate(_object); // baseline may have changed
                if ( path_taboo ) // add to path_taboo
                  path_taboo->add(next_o, false);
                if ( current_e > next_e )
                { // yeah! gradient was right
                  is_not_converged = true;
                  current_e = next_e;
                  _object = next_o;
                }
                else // flips back -- gradient was wrong
                  *i_var = ( *i_var > 0.0 ) ? -1.0 : 1.0; 
              }
              else // flips back -- gradient says dont investigate
                *i_var = ( *i_var > 0.0 ) ? 1.0 : -1.0; 
            } // end of check for taboo

            if ( count >= max_directions_checked )  // breaks out after n calls to functional
              { is_not_converged = false; break; }

          } // loop over shuffled directions

        } while ( is_not_converged );

//       _object = copy;
        eval(_object);
        return false; // fitness is set, returns false
      } // end of operator()(t_Object &_object)

  };

} // namespace LaDa
#endif
