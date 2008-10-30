//
//  Version: $Id: debug.h 754 2008-08-20 22:56:21Z davezac $
//
#ifndef __OPT_BPO_MACROS__
#define __OPT_BPO_MACROS__

#include <boost/program_options.hpp>

# define __BPO_START__ \
  namespace po = boost::program_options; \
  po::options_description generic("Generic Options"); \
  generic.add_options() \
         ("help,h", "produces this help message.") \
         ("version,v", "prints version string.")
# define __BPO_SPECIFICS__( section_string ) \
  po::options_description specific( section_string ); \
  specific.add_options() \
      ("input,i", po::value<std::string>()->default_value("input.xml"), \
       "input filename." ) 
# define __BPO_RERUNS__ \
      ("reruns,r", po::value<unsigned>()->default_value(1), \
                   "number of times to run the algorithm.\n" \
                   "Is equivalent to manually re-launching the program.\n") 

# define __BPO_GENERATE__( others ) \
  po::options_description all; \
  all.add(generic).add(specific) others ; \
  po::positional_options_description p; \
  p.add("input", 1); 

# define __BPO_MAP__ \
  po::variables_map vm; \
  po::store(po::command_line_parser(argc, argv). \
            options(all).positional(p).run(), vm); \
  po::notify(vm); \

# define __BPO_HELP_N_VERSION__ \
  std::cout << "\n" << __PROGNAME__ \
            << " from the " << PACKAGE_STRING << " package.\n" \
            << "Subversion Revision: " << SVN::Revision << "\n\n"; \
  if ( vm.count("help") ) \
  { \
    __ROOTCODE( (*::mpi::main), \
      std::cout << "Usage: " << argv[0] << " [options] file.xml\n" \
                << "  file.xml is an optional filename for XML input.\n" \
                << "  Default input is input.xml.\n\n" \
                << all << "\n";  \
    ) \
    return 1; \
  } \
  if ( vm.count("version") ) \
  { \
    __ROOTCODE(  (*::mpi::main), \
      std::cout << "\n" << __PROGNAME__ \
                << " from the " << PACKAGE_STRING << " package\n" \
                << "Subversion Revision: " << SVN::Revision << "\n\n";  \
    ) \
    return 1; \
  }

#define __BPO_CATCH__(code) \
  } \
  catch ( boost::program_options::invalid_command_line_syntax &_b)\
  {\
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"\
              << "Something wrong with the command-line input.\n"\
              << _b.what() << std::endl;\
    code; \
    return 0;\
  }\
  catch ( boost::program_options::invalid_option_value &_i )\
  {\
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"\
              << "Argument of option in command-line is invalid.\n"\
              << _i.what() << std::endl;\
    code; \
    return 0;\
  }\
  catch ( boost::program_options::unknown_option &_u)\
  {\
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"\
              << "Unknown option in command-line.\n"\
              << _u.what() << std::endl;\
    code; \
    return 0;\
  }\
  catch (  boost::program_options::too_many_positional_options_error &_e )\
  {\
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"\
              << "Too many arguments in command-line.\n"\
              << _e.what() << std::endl;\
    code; \
    return 0;\
  }\
  catch ( std::exception &e )\
  {\
    std::cout << "Caught error while running " << __PROGNAME__ \
              << "\n" << e.what() << std::endl;\
    code; \
    return 0;\
  }


#else

# ifdef __BPO_CATCH__
#   undef __BPO_CATCH__
# endif
# ifdef __BPO_HELP_N_VERSION__
#   undef __BPO_HELP_N_VERSION__
# endif

#endif
