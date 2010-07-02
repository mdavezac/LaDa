#include <stdlib.h>

#include <boost/xpressive/basic_regex.hpp>
#include <boost/xpressive/regex_algorithms.hpp>
#include <boost/xpressive/regex_primitive.hpp>
#include <boost/xpressive/regex_actions.hpp>


#include "path.h"

namespace LaDa
{
  namespace opt
  {
    struct GetEnv
    {
      typedef std::string const result_type;

      result_type operator()(std::string const &val) const
        { return getenv(val.c_str()); }
    };
    struct GetHome
    {
      typedef std::string const result_type;

      result_type operator()(std::string const &val) const
      {
        char *orig_ = get_current_dir_name();
        chdir(val.c_str());
        char *result_ = get_current_dir_name();
        std::string const result = result_;
        free(result_);
        chdir(orig_.c_str());
        free(orig_);
        return result;
      }
    };


    boost::filesystem::path const expandpath(boost::filesystem::path const &_path) 
    {
      namespace bx = boost::xpressive
      std::string path = _path;   
      bx::sregex home = (bx::beos|bx::as_xpr('~')) >> +(bx::eos|~bx::as_xpr('/'));
      bx::sregex var1 = bx::as_xpr('$') >> (bx::s1_ = +(~bx::as_xpr('/')));
      bx::sregex var2 = bx::as_xpr("$(") >> (bx::s1_ = +(~bx::as_xpr('/'))) >> bx::as_xpr(')');
      bx::cmatch what;
      bx::function<GetEnv>::type const rep_env = {{}};
      bx::function<GetEnv>::type const rep_home = {{}};


      path = bx::regex_replace(path, var1, env(bx::as<std::string const>(bx::s1_)));
      path = bx::regex_replace(path, var2, env(bx::as<std::string const>(bx::s1_)));
      path = bx::regex_replace(path, home, rep_home(bx::as<std::string const>(bx::s1_)));
    }
     

  }
}
#endif
