//
//  Version: $Id: main.extras.h 868 2008-11-22 02:02:34Z davezac $
//


#if ! defined( _LADA_GA_MAIN_EXTRAS_H_ )
#  define _LADA_GA_MAIN_EXTRAS_H_ 0
#  include<boost/lexical_cast.hpp> 
   template< class T_TYPE > T_TYPE call_lexcast( const std::string& _str )
    { return boost::lexical_cast<T_TYPE>( _str ); }
    inline LaDa::types::t_unsigned call_lexcast_unsigned( const std::string &_str )
      { return boost::lexical_cast<LaDa::types::t_unsigned>( _str ); }
    inline LaDa::types::t_real call_lexcast_real( const std::string &_str )
      { return boost::lexical_cast<LaDa::types::t_real>( _str ); }

#elif _LADA_GA_MAIN_EXTRAS_H_ == 0 || _LADA_GA_MAIN_EXTRAS_H_ == 1
    // Creates factories.
#  if _LADA_GA_MAIN_EXTRAS_H_ == 0
#    undef _LADA_GA_MAIN_EXTRAS_H_
#    define _LADA_GA_MAIN_EXTRAS_H_ 1
#  else
#    undef _LADA_GA_MAIN_EXTRAS_H_
#    define _LADA_GA_MAIN_EXTRAS_H_ 2
#  endif

  typedef LaDa::Factory::Factory
          < 
            void( boost::function< bool( const t_Individual& ) >&, const TiXmlElement& ), 
            std::string
          > t_TabooFactory;
  typedef t_Evaluator :: t_GATraits :: t_Population t_Population;
  typedef t_Evaluator :: t_GATraits :: t_Islands t_Islands;

  ga.taboo_factory.connect
    (
      "Population", "No two individuals can be the same in the population.",
      boost::bind( &LaDa::GA::Taboo::Factory::islands<t_Islands>, 
                   _1, _2, boost::cref( ga.islands ) )
    )
    (
      "Offspring", "No two individuals can be the same in the offspring population.",
      boost::bind( &LaDa::GA::Taboo::Factory::offspring<t_Population>, 
                   _1, _2, boost::cref( ga.offspring ) )
    );

  typedef LaDa :: GA :: Factory :: XmlOperators< t_Individual, eoPopulator<t_Individual> >
    t_OpFactory;

  ga.operator_factory.connect_attribute
    ( "period", "operators with this attribute will be"
      "called only every n=\"period\" generations.",
       boost::bind( &LaDa::GA::Factory::periodic< t_OpFactory >,
                     _1, _2, _3, ga.get_counter() ) );
  ga.operator_factory.connect
    ( "Operators", "container of operators."
      "Accepts attribute type=\"and\" or type=\"or\".",
      boost::bind( &LaDa::GA::Factory::containers<t_OpFactory>, _1, _2, _3 ) ) 
    ( "And", "Sequential container."
      "same as <Operators type=\"and\"> ... <Operators/>.",
      boost::bind( &LaDa::GA::Factory::sequential<t_OpFactory>, _1, _2, _3 ) ) 
    ( "Or", "Proportional container."
      "same as <Operators type=\"or\"> ... <Operators/>.",
      boost::bind( &LaDa::GA::Factory::sequential<t_OpFactory>, _1, _2, _3 ) )
    ( "TabooOp", "An offspring created by the operators within this one"
      " will be rejected if it is in the tabooed individuals"
      " defined in <Taboos> .. </Taboos>. Offspring creation"
      " will repeat for at most n (attribute max=\"n\") trials,"
      " at which point a random individual is created.",
      boost::bind( &LaDa::GA::Operator::taboo_factory<t_OpFactory>, 
                   _1, _2, _3, boost::ref(ga.taboos), "Random" ) );

    ga.tournament_size = 2;
    ga.replacement_rate = 0.5;
    ga.pop_size = 1;
    ga.max_generations = 0;
    ga.population_creator
      = LaDa::GA::Operators::populator
        ( 
          bl::bind( &t_Evaluator::initialize, boost::cref(ga.evaluator), bl::_1 ),
          ga.taboos 
        );
    ga.att_factory.connect
      (
        "tournament", "size of the deterministic tournaments"
        " for each parent selection. Default = 2.",
        bl::var(ga.tournament_size)
           = bl::bind<LaDa::types::t_unsigned> 
                     ( &call_lexcast_unsigned, bl::_1 )
      )
      (
        "rate", "replacement rate (eg offspring creation rate)."
        " Default = 0.5.",
        bl::var(ga.replacement_rate) = bl::bind<LaDa::types::t_real>( &call_lexcast_real, bl::_1 ) 
      )
      (
        "popsize", "population size. Default = 1.",
        bl::var(ga.pop_size) = bl::bind( &call_lexcast_unsigned, bl::_1 ) 
      )
      (
        "maxgen", "Max number of generations before quitting."
        "Default = 0.",
        (
          bl::bind
          (
            LaDa::GA::CheckPoint::AddressOf::maxgenerations( ga.checkpoints ),
            bl::var(ga.checkpoints), bl::_1, bl::constant(ga.counter) 
          ),
          bl::var(ga.max_generations) = bl::bind( &call_lexcast_unsigned, bl::_1 )
        ) 
      )
      (
        "maxeval", "Max number of evaluations before quitting."
        "Default = 0.",
        boost::bind
        (
          LaDa::GA::CheckPoint::AddressOf::maxevaluations( ga.checkpoints, ga ),
          boost::ref(ga.checkpoints), _1, boost::cref( ga )
        )
      )
      (
        "islands", "Number of independent islands. Default = 1.",
        bl::var(ga.nb_islands) = bl::bind( &call_lexcast_unsigned, bl::_1 ) 
      )
      (
        "starting_population_only", "Quits after evaluating starting population."
        "Default = false.",
        boost::lambda::var( ga.do_starting_population_only ) = true
      )
      (
        "seed", "random seed. Default = 0."
        " In MPI, seed0, seed1, ... are accepted, if not necessarily used.",
        boost::bind( &LaDa::GA::Topology::seed_n, ga.topology, 0, _1 )
      );
#    ifdef _MPI
        for( size_t ii(0); ii < world->size(); ++ii )
          ga.att_factory.connect
          (
               "seed" + boost::lexical_cast<std::string>(ii),
               "random seed for process " + boost::lexical_cast<std::string>(ii) 
               + ". Default = 0.",
            boost::bind( &LaDa::GA::Topology::seed_n, ga.topology, 0, _1 )
          );
#     endif

       ga.restart_filename = input;
       ga.evaluator_filename = input;
       ga.checkpoint_factory.connect_value
         (
           LaDa::Factory::XmlKey( "Print", "type", "population" ),
           "Prints population at end of generation.",
           std::ptr_fun( 
             LaDa::GA::CheckPoint::AddressOf::print_population( ga.checkpoints )
           )
         )
         (
           LaDa::Factory::XmlKey( "Print", "type", "new" ),
           "Prints new individuals in population at end of generation.",
           boost::bind
           (
             LaDa::GA::CheckPoint::AddressOf::print_newindividuals( ga.checkpoints ),
             _1, boost::cref( ga.counter )
           )
         )
         (
           LaDa::Factory::XmlKey( "Print", "type", "average" ),
           "Prints average (scalar) fitness of individuals in population at end of generation.",
           std::ptr_fun( LaDa::GA::CheckPoint
                             ::AddressOf
                             ::print_averagefitness( ga.checkpoints, ga.offspring ) )
         )
         (
           LaDa::Factory::XmlKey( "Print", "type", "census" ),
           "Prints the number of different individuals in the population.",
           std::ptr_fun( LaDa::GA::CheckPoint
                             ::AddressOf
                             ::print_truecensus( ga.checkpoints, ga.offspring ) )
         );
       ga.checkpoint_factory.connect_attribute
       (
         LaDa::Factory::XmlKey( "Print", "xyz" ),
         "Prints best current result(s) structure in xyz format to file.",
         boost::bind
         ( 
           LaDa::GA::CheckPoint::AddressOf::xyz_anim( ga.checkpoints, ga ),
           _1, _2, boost::cref( ga ), boost::cref( ga.evaluator.get_structure() )
         )
       )
       (
         LaDa::Factory::XmlKey( "Print", "xcrysden" ),
         "Prints best current result(s) structure in xcrysden format to file.",
         boost::bind
         ( 
           LaDa::GA::CheckPoint::AddressOf::xcrysden_anim( ga.checkpoints, ga ),
           _1, _2, boost::cref( ga ), boost::cref( ga.evaluator.get_structure() )
         )
       )
       (
         LaDa::Factory::XmlKey( "Filenames", "stop" ),
         "GA stops when this file is found.",
         std::ptr_fun( LaDa::GA::CheckPoint::AddressOf::stop_onfile( ga.checkpoints ) )
       ) 
       (
         LaDa::Factory::XmlKey( "Filenames", "evaluator" ),
         "Reads evaluators from this file.",
         LaDa::GA::CheckPoint::Factory::Dummy()
       )
       (
         LaDa::Factory::XmlKey( "Filenames", "restart" ),
         "Reads restart data from this file.",
         boost::bind
         (
           LaDa::GA::CheckPoint::AddressOf::filename( ga.checkpoints ),
           _1, _2, "Will restart from data in ", boost::ref(ga.restart_filename)
         )
       )
       (
         LaDa::Factory::XmlKey( "Filenames", "xmgrace" ),
         "Xmgrace formatted output goes to this file: ",
         LaDa::GA::CheckPoint::Factory::Dummy()
       )
       (
         LaDa::Factory::XmlKey( "Filenames", "out" ),
         "Writes output to this file.",
         LaDa::GA::CheckPoint::Factory::Dummy()
       );
       ga.checkpoint_factory.connect_node
       (
         LaDa::Factory::XmlKey( "Filenames" ),
         "when attributes save=\"filename\" and what=\"all or history or "
         "population or results\" are set, saves the relevant information "
         "every n generations, with n specified by  attribute every=\"n\""
         "(Default n=0, eg each and every generation).",
         boost::bind
         (
           LaDa::GA::CheckPoint::AddressOf::save_every( ga.checkpoints, ga ),
           _1, _2, boost::cref(ga) 
         )
       );



#     if _LADA_GA_MAIN_EXTRAS_H_ == 2 
        LaDa::GA::CheckPoint::Factory::print_at_iteration( ga.checkpoints, true, ga );
#     endif
//     ga.checpoint_factory.connect_
#else 
   // makes sure this file cannot be (meaningfully) included anymore.
#  undef _LADA_GA_MAIN_EXTRAS_H_
#  define _LADA_GA_MAIN_EXTRAS_H_ nofurtherinclude
#endif

