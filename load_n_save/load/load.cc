#include "LaDaConfig.h"

#include <boost/exception/get_error_info.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include "load.h"
#include "../exceptions.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace load
    {
      bool Load::operator()( tree::Base const& _tree,
                             xpr::Section const& _sec,
                             version_type _version ) const
      {
        tree::Base base(_tree);
        tree::Section text("");
        text.tree::Base::swap(base);
        try { return Section(text, verbose, _version)(_sec); }
        catch(error::internal &_e)
        {
          std::cerr << "Caught internal error:\n" << boost::diagnostic_information(_e);
        }
        catch(error::too_few_sections &_e)
        {
          std::string const * option_name = boost::get_error_info<error::option_name>(_e);
          std::string const * section_name = boost::get_error_info<error::section_name>(_e);
          std::cerr << boost::diagnostic_information(_e) << "\n" 
                    << "Too few sections of following kind:\n";
          if(section_name) std::cerr << "section: " << *section_name << "\n";
          if(option_name) std::cerr << "option: " << *option_name << "\n";
        }
        catch(error::required_section_not_found &_e)
        {
          std::string const * section_name = boost::get_error_info<error::section_name>(_e);
          std::cerr << boost::diagnostic_information(_e) << "\n" 
                    << "Following required section not found.\n";
          if(section_name) std::cerr << "section: " << *section_name << "\n";
        }
        catch(error::section_parse_error &_e)
        {
          std::string const * section_name = boost::get_error_info<error::section_name>(_e);
          std::cerr << boost::diagnostic_information(_e) << "\n" 
                    << "Following section did not parse correctly.\n";
          if(section_name) std::cerr << "section: " << *section_name << "\n";
          return false;
        }
        catch(error::section &_e)
        {
          std::string const * section_name = boost::get_error_info<error::section_name>(_e);
          std::cerr << boost::diagnostic_information(_e) << "\n" 
                    << "Unknown error caught when parsing following section.\n";
          if(section_name) std::cerr << "section: " << *section_name << "\n";
          return false;
        }
        catch(error::required_option_not_found &_e)
        {
          std::string const * option_name = boost::get_error_info<error::option_name>(_e);
          std::string const * section_name = boost::get_error_info<error::section_name>(_e);
          std::cerr << boost::diagnostic_information(_e) << "\n" 
                    << "Following required option not found.\n";
          if(section_name) std::cerr << "section: " << *section_name << "\n";
          if(option_name) std::cerr << "option: " << *option_name << "\n";
        }
        catch(error::option_parse_error &_e)
        {
          std::string const * option_name = boost::get_error_info<error::option_name>(_e);
          std::string const * section_name = boost::get_error_info<error::section_name>(_e);
          std::cerr << boost::diagnostic_information(_e) << "\n" 
                    << "Following option did not parse correctly.\n";
          if(section_name) std::cerr << "section: " << *section_name << "\n";
          if(option_name) std::cerr << "option: " << *option_name << "\n";
          return false;
        }
        catch(error::option &_e)
        {
          std::string const * option_name = boost::get_error_info<error::option_name>(_e);
          std::string const * section_name = boost::get_error_info<error::section_name>(_e);
          std::cerr << boost::diagnostic_information(_e) << "\n" 
                    << "Unknown error caught when parsing following option.\n";
          if(section_name) std::cerr << "section: " << *section_name << "\n";
          if(option_name) std::cerr <<  "option: " << *option_name << "\n";
          return false;
        }
        return false;
      }
      bool Load::operator()( tree::Section const& _tree,
                             xpr::Section const& _sec,
                             version_type _version ) const
        { return Section(_tree, verbose, _version)(_sec); }

    } // namespace initializer.
  } // namespace load_n_save
} // namespace LaDa
