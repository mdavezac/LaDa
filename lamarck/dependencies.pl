#! /usr/bin/perl
#

my $computer = "lester";
my %params;

my $HOME = `cd; pwd`; chomp $HOME;

@{$params{"defs"}} = ( "_MPI" );

@{$params{"Includes"}} = (".", "$HOME/usr/src/atat/src",
                          "$HOME/usr/include/opt", "$HOME/usr/include/analysis");

if ( $computer =~ /home/ )
{
  @{$params{"make include"}} = ( "atat", "$HOME/usr/include");
  @{$params{"make lib"}} = ( ); 
  $params{"CC"}  = "gcc";
  $params{"CXX"} = "g++";
  $params{"LD"}  = "g++";
  $params{"F77"}  = "g77";
  $params{"CXXFLAGS"}  = "-mtune=athlon64 -ffriend-injection";
}
elsif ( $computer =~ /office/ )
{
  @{$params{"make include"}} = ( "$HOME/usr/include"
                                 ,"/opt/mpich/include"
                               );
  @{$params{"make lib"}} = ( "-lm", "-lstdc++", "-L $HOME/usr/lib/", "-llamarck", "-latat", 
                             "-L /opt/mpich/ch-p4/lib/", "-lpmpich++", "-lpmpich", "-lmpiobject", 
                             "-ltinyxml" );
  $params{"CC"}  = "gcc";
  $params{"CXX"} = "gcc";
  $params{"LD"}  = "gcc";
  $params{"F77"}  = "g77";
  $params{"CXXFLAGS"}  = "-malign-double -ffriend-injection";
}
elsif ( $computer =~ /lester/ )
{
  @{$params{"make include"}} = ( "atat", 
                                 "/opt/mpich.gcc/include",
                                 "$HOME/usr/include/" );
  @{$params{"make lib"}} = ( "-lm", "-lstdc++", "-L $HOME/usr/lib/", "-llamarck", "-latat", 
                             "-L /opt/mpich.gcc/lib/", "-lpmpich++", "-lpmpich", "-lmpiobject", 
                             "-ltinyxml" );
  $params{"CC"}  = "gcc";
  $params{"CXX"} = "g++";
  $params{"LD"}  = "g++";
  $params{"F77"}  = "f77";
  $params{"CXXFLAGS"}  = "-mtune=opteron";
}
elsif ( $computer =~ /super/ )
{
  @{$params{"make include"}} = ( "atat", "$HOME/usr/include");
  @{$params{"make lib"}} = ();
  $params{"CC"}  = "icc";
  $params{"CXX"} = "icc";
  $params{"LD"}  = "icc";
  $params{"F77"}  = "icc";
  $params{"CXXFLAGS"}  = "-mtune=itanium2";
}



# find source files in present directory
my %dependencies;
foreach $file ( <*.cc> )
{
  if ( $file =~ /lada/ )
    { next; }
  if ( $file =~ /fftw_interface/ )
    { next; }

  
  my $key = $file; $key =~ s/\.cc//;
  if ( $key =~ /\// )
  {
    $_ = $key; /(\S+)\/(\S+)/; $key = $2;
    $dependencies{$key}{"location"} = $1;
  }
  else
  {
    $dependencies{$key}{"location"} = ".";
  }
  if ( -e "$dependencies{$key}{'location'}/$key" . ".cc" )
    { $dependencies{$key}{"source"} = 1; }
  if ( -e "$dependencies{$key}{'location'}/$key" . ".h" )
  {
    $dependencies{$key}{"header"} = 1;
    push @{$dependencies{$key}{"depends on"}}, ( $key );
  }
}

while ( get_dependencies() == 0 ) {};

copy_dependencies();
write_make_file();
write_dependencies();


exit;

# now finds dependencies for each key
# if a new key is added, then its dependencies are also investigated
sub get_dependencies()
{
  my $key;
  my $new_stuff = 1;

  foreach $key ( keys %dependencies )
  {
    if ( exists $dependencies{$key}{"header"} ) 
    {
      open IN, "<$dependencies{$key}{'location'}/$key.h"
        or die "file $dependencies{$key}{'location'}/$key.h not found.\n";
    
      while ( <IN> )
      {
        if (/\#include(\s+|)(\"|\<)(\S+\/|)(\S+)\.h(\"|\>)/)
        {
          my $new_key = $4;
          if ( !( exists $dependencies{$new_key} ) )
          {
            foreach $location (  @{$params{"Includes"}} ) 
            { 
              if ( -e "$location/$new_key.h" )
              { 
                $dependencies{$new_key}{"location"} = $location;
                $dependencies{$new_key}{"header"} = 1;
                if ( -e "$location/$new_key.cc" )
                  { $dependencies{$new_key}{"source"} = 1; }
              }
            }
            if( exists $dependencies{$new_key} ) 
              { $new_stuff = 0; }

#           if( not exists $dependencies{$new_key} ) 
#           { print " location of $new_key not found\n"; }
          } # end if ( !( exists $dependencies{$new_key} ) )
          
          add_to_dependencies( $key, $new_key ); 

        } # end if (/\#include(\s+|)\"(\S+)\.h\"/)
      } # end while ( <IN> )
      close IN;
    } # end if ( exists $dependencies{$key}{"header"} ) 
    if ( exists $dependencies{$key}{"source"} ) 
    {
      open IN, "<$dependencies{$key}{'location'}/$key.cc"
        or die "file $dependencies{$key}{'location'}/$key.cc not found.\n";
      while ( <IN> )
      {
        if (/\#include(\s+|)(\"|\<)(\S+\/|)(\S+)\.h(\"|\>)/)
        {
          my $new_key = $4;
          if ( !( exists $dependencies{$new_key} ) )
          {
            foreach $location (  @{$params{"Includes"}} ) 
            { 
              if ( -e "$location/$new_key.h" )
              { 
                $dependencies{$new_key}{"location"} = $location;
                $dependencies{$new_key}{"header"} = 1;
                if ( -e "$location/$new_key.cc" )
                  { $dependencies{$new_key}{"source"} = 1; }
              }
            }
            if( exists $dependencies{$new_key} ) 
              { $new_stuff = 0; }
#           if( not exists $dependencies{$new_key} ) 
#           { print " location of $new_key not found\n"; }
           
          } # end if ( !( exists $dependencies{$new_key} ) )

          add_to_dependencies( $key, $new_key );

        } # end if (/\#include(\s+|)\"(\S+)\.h\"/)
      } # end while (<IN>)
      close IN;
    } # end  ( exists $dependencies{$key}{"source"} ) 
  }

  return $new_stuff;
}

sub template()
{
  open OUT, ">make_file_template" or die;
  printf OUT "CC     := $params{'CC'}\nCXX    := $params{'CXX'}\nF77    := $params{'F77'}\nLD     := $params{'LD'}\n";
  printf OUT "AR     := ar rvu\nRANLIB := ranlib\nDEBUG = NO\n";
  printf OUT "\nFFLAGS     := \n";
  printf OUT  "CFLAGS     := \n";
  printf OUT  "CXXFLAGS   := $params{'CXXFLAGS'} \n";
  printf OUT "DEFS       := ";
  foreach $def ( @{$params{"defs"}} )
    { print OUT "-D ", $def, " "; }
  print OUT "\n\n";


  printf OUT "\nDEBUG_CFLAGS     := -Wall -Wno-format -g -O0 -Wno-unknown-pragmas -fbounds-check\n";
  printf OUT "DEBUG_CXXFLAGS    := \${DEBUG_CFLAGS} -D_DEBUG_LADA_ \n";
  printf OUT "DEBUG_LDFLAGS     := -g \n";
  printf OUT "\nRELEASE_CXXFLAGS := \${RELEASE_CFLAGS}\n";
  printf OUT "RELEASE_LDFLAGS   := \n";
  printf OUT "RELEASE_CFLAGS    := -Wall -Wno-unknown-pragmas -Wno-format -O3\n";
  printf OUT "\nLIBS:=  ";
  foreach $include (@{$params{"make lib"}})
    { print OUT " $include "; } 
  print OUT "\n";
  printf OUT "INCS := ";
  foreach $lib (@{$params{"make include"}})
    { print OUT "-I ", $lib, " "; } 
  print OUT "\n";
  printf OUT "\n\nifeq (YES, \${DEBUG})\n";
  printf OUT "   CFLAGS      := \${CFLAGS} \${DEBUG_CFLAGS} \n";
  printf OUT "   CXXFLAGS    := \${CXXFLAGS} \${DEBUG_CXXFLAGS}\n";
  printf OUT "   LDFLAGS     := \${LDFLAGS} \${DEBUG_LDFLAGS}\n";
  printf OUT "else\n";
  printf OUT "   CFLAGS      := \${CFLAGS} \${RELEASE_CFLAGS} -O3\n";
  printf OUT "   CXXFLAGS    := \${CXXFLAGS} \${RELEASE_CXXFLAGS} -O3 \n";
  printf OUT "   LDFLAGS     := \${LD_FLAGS} \${RELEASE_LDFLAGS} -O3\n";
  printf OUT "endif\n";
  printf OUT "\nCFLAGS   := \${CFLAGS}   \${DEFS}\n";
  printf OUT "CXXFLAGS := \${CXXFLAGS} \${DEFS}\n";
  printf OUT "\nOUTPUT := liblamarck.a\n";
  printf OUT "\nall: atat \$(OUTPUT) test \n";
  printf OUT "\nSRCS := main.cc\n";
  printf OUT "\n\?LIBSRCS :=\?\n";
  printf OUT "\n\?ATATSRCS :=\?\n";
  printf OUT "\nOBJS := \$(addsuffix .o,\$(basename \${SRCS}))\n";
  printf OUT "\nLIBOBJS := \$(addsuffix .o,\$(basename \${LIBSRCS}))\n";
  printf OUT "\nATATOBJS := \$(addsuffix .o,\$(basename \${ATATSRCS}))\n";
  printf OUT "\n.PHONY: clean cleanall\n";
  printf OUT "\n\${OUTPUT}: lib atat \n\n";
# printf OUT "\t\${LD} \${LDFLAGS} -o \$@ \${OBJS} ";
# printf OUT "-L . -L atat -llamarck -latat \${LIBS} \${EXTRALIBS}\n";
  printf OUT "\n\${OBJS} : \n";
  printf OUT "\t\${CXX} -c \${CXXFLAGS} \${INCS} \$< -o \$@\n\n";
  printf OUT "\n\${ATATOBJS} : \n";
  printf OUT "\t\${CXX} -c \${CXXFLAGS} $params{'cxx'}{'atat flags'} \${INCS} \$< -o \$@\n\n";
  printf OUT "\n\${LIBOBJS} : \n";
  printf OUT "\t\${CXX} -c \${CXXFLAGS} $params{'cxx'}{'lamarck flags'} \${INCS} \$< -o \$@\n\n";
  printf OUT "\?dependencies\?";
  printf OUT "\ntest: \${OBJS} lib atat\n";
  printf OUT "\t\${LD} \${LDFLAGS} -o lamarck \${OBJS}  \${LIBS} \${EXTRALIBS}\n";

            

  printf OUT "\n\nlib: \${LIBOBJS}\n";
  if ( exists $params{'ar'}{'lamarck flags'} ) 
  {
    printf OUT "\t%s \n", $params{'ar'}{'lamarck flags'}[0];
    printf OUT "\t\${AR} \${OUTPUT} \${LIBOBJS} %s/*.o\n",
               $params{'ar'}{'lamarck flags'}[1];
  }
  else 
  {
    printf OUT "\t\${AR} \${OUTPUT} \${LIBOBJS} \n";
  }
  printf OUT "\t\${RANLIB} \${OUTPUT}\n";

  printf OUT "\n\natat: \${ATATOBJS}\n";
  if ( exists $params{'ar'}{'atat flags'} ) 
  {
    printf OUT "\t%s \n", $params{'ar'}{'atat flags'}[0];
    printf OUT "\t\${AR} atat/lib\$@.a \${ATATOBJS} %s/*.o\n",
               $params{'ar'}{'atat flags'}[1];
  }
  else 
  {
    printf OUT "\t\$\{AR} atat/lib\$@.a \${ATATOBJS} \n";
  }
  printf OUT "\t\${RANLIB} atat/lib\$@.a\n";

  printf OUT "\n\nclean:\n\t- rm -f \${OBJS} \${LIBOBJS}\n";
  printf OUT "\t- rm -f \${OUTPUT}\n\n";

  printf OUT "\n\ncleanall:\n\t- rm -f \${OBJS} \${LIBOBJS} \${ATATOBJS} \${TIXMLOBJS}\n";
  printf OUT "\t- rm -f \${OUTPUT}\n";
  if (exists $params{'cleanall'} )
  { 
    foreach $clean ( (@{$params{"cleanall"}}) )
    {
      printf OUT "\t%s\n", $clean;
    }
  }
  printf OUT "\t- rm -f atat/libatat.a \${OUTPUT} \n\n";

  printf OUT "\n\ninclude .dependencies\n\n"; 


  close OUT;
}



sub copy_dependencies()
{
  system "rm -f atat/*";
  foreach $key ( keys %dependencies )
  {
    if ( $dependencies{$key}{"location"} =~ /atat/ )
    {
      if ( exists $dependencies{$key}{"header"} ) 
        { copy_atat($key,"header") }
      if ( exists $dependencies{$key}{"source"} ) 
        { copy_atat($key,"source") }
    }
  }
}

sub copy_atat($$)
{
  my $key = $_[0];
  my $type = $_[1];
  my $suffix = ".cc";
  if ($type =~ /header/ )
    { $suffix = ".h"; }


  open IN, "<$dependencies{$key}{'location'}/$key$suffix" 
    or die "Could not open $dependencies{$key}{'location'}/$key$suffix\n";

  my $last_include = 0;
  my $last_endif = 0;
  my $counter = 0;
  my $first_other = 100;

  if ( $type =~ /header/ and $counter == 0 )
  {
    my $def;
    while ( ($_=<IN>) && ($_!~/\#ifndef(\S+|)(\s+)/ ) ) 
      { $counter++; }
    $def = $2; chomp $def;
    do
      { $counter++; }
    while ( ($_=<IN>) && ($_!~/\#define(\S+|)$def/ ) );
    $last_include = $counter + 1;
  }
  while( <IN> )
  {
    if( /\#include/ )
      { $last_include = $counter + 2; }
    elsif( /RESULT_FORALLi/ 
           and $first_other > $counter)
      { $first_other = $counter; }
    elsif( /(class|template|void|int|Real)/ 
           and $first_other > $counter)
      { $first_other = $counter; }
    elsif( /\#endif/ )
      { $last_endif = $counter; }
    $counter++;
  }

  if ( $last_endif < $counter - 5 )
   { $last_endif = $counter -1;}
  if ( $last_endif == 0 )
   { $last_endif = $counter -1;}
  if ( $first_other < $last_include )
   { $last_include = $first_other; }

  open OUT, ">atat/$key$suffix" 
    or die "Could not open $dependencies{$key}{'location'}/$key$suffix\n";

  seek (IN, 0, 0);
  $counter = 0;
  while( <IN> )
  {
    if ( $counter == $last_include )
    {
      printf OUT "\n#include <opt/types.h>\n\nnamespace atat\n{ \n";
      if ($key =~ /xtalutil/ )
        { printf OUT "\nusing ::ceil;using ::floor;\n"; }
    }
    elsif ( $counter == $last_endif )
    {
      printf OUT "\n\n} \/\/ namespace atat\n";
    }
    elsif( $key=~/arraylist/ and $suffix =~/\.h/ and
           $computer =~/lester/ and 
           /operator void \* \(\) {return \(void \*\)valid;}/ )
    {
      $_ = "  operator void * () { return valid ? (void*)(1) : (void*)(NULL); }\n";
    }
    if( $_=~/int/ and $_ !~ /typedef+\s+(unsigned\s+|)(long\s+|short\s+|)int/ )
    {
      $_ =~ s/\bint\b/types::t_int/cg;
    }
    if( $_=~/operator(\s+|)\+\+\((\s+|)types::t_int(\s+|)\)/ )
    {
      $_ =~
      s/operator(\s+|)\+\+\((\s+|)types::t_int(\s+|)\)/operator\+\+(int)/cg;
    }
    if( $_=~/unsigned/ and $_ !~/long/ and $_ !~ /typedef\s+unsigned/)
    {
      $_ =~ s/\bunsigned\b/types::t_unsigned/cg;
    }
    if( $_=~/Real/ and  $_!~/#define Real double/)
    {
      $_ =~ s/\bReal\b/types::t_real/cg;
    }
    if( $_=~/double/ and  $_!~/#define Real double/)
    {
      $_ =~ s/\bdouble\b/types::t_real/cg;
    }
    if ( $key =~ /vectmac/  and
         $_ =~
         /\#define\s+(\S)(Vector|Matrix|BoundingBox)(\S\S)\s(Vector|Matrix|BoundingBox)(\S\S)(\S*)/  )
    {   
       $_ = "typedef $4$5$6 $1$2$3;\n";
    }
    print OUT $_;
    $counter++;
  }
 
  close IN;
  close OUT;
}

sub  write_make_file()
{
  my @sorted_keys = sort{ sort_hash($a, $b) } keys %dependencies;
  template();

  open IN, "make_file_template" or die;
  open OUT, ">makefile" or die;
  
  while ( ($_=<IN> ) )
  {
    if ( /\?LIBSRCS :=\?/ )
    {
      print OUT 'LIBSRCS := ';
      my $i = 0;
      foreach $key ( @sorted_keys )
      {
        if ( $key =~ /(main|lamarck)/ )
          { next; }
        if ( exists $dependencies{$key}{"source"} )
        {
          if ( $dependencies{$key}{"location"} =~ /atat/ )
            { next; }

          if ( $i % 5 == 0 and $i != 0 )
            { print OUT "\\\n\t"; }

          if ( $dependencies{$key}{"location"} eq "." )
            { print OUT $key, ".cc "; $i++; }
          else
            { print OUT $dependencies{$key}{"location"},"/", $key, ".cc "; $i++; }
        }
      }
      print OUT "\n";
      next;
    }
    elsif ( /\?ATATSRCS :=\?/ )
    {
      print OUT 'ATATSRCS := ';
      my $i = 0;
      foreach $key ( @sorted_keys )
      {
        if ( exists $dependencies{$key}{"source"} )
        {
          if ( $i % 5 == 0 and $i != 0 )
            { print OUT "\\\n\t"; }
          if ( $dependencies{$key}{"location"} !~ /atat/ )
            { next; } 
          print OUT "atat/", $key, ".cc "; $i++;
        }
      }
      print OUT "\n";
      next;
    }
    elsif (/\?dependencies\?/)
    { next; }
    print OUT $_;
  }
  
  close IN;
  close OUT;
  system "rm make_file_template";

}

sub add_to_dependencies(\$\$)
{
  my $key = $_[0];
  my $new_key = $_[1];

  if ( ! exists $dependencies{$key}{"depends on"} )
  {
    push @{$dependencies{$key}{"depends on"}}, ( "$new_key" ); 
    return;
  }

  $w = 0;
  foreach $dep ( @{$dependencies{$key}{"depends on"}} )
  {
     if ( $dep eq $new_key )
       { last; }
     $w ++ ;
  }
  
  if ( $w == scalar( @{$dependencies{$key}{"depends on"}} ) )
    { push @{$dependencies{$key}{"depends on"}}, ( "$new_key" ); }
}


sub write_dependencies()
{
  open OUT, ">.dependencies"
    or die " could not open .dependencies for writing \n";

  my @sorted_keys = sort{ sort_hash($a, $b) } keys %dependencies;
  foreach $key ( @sorted_keys )
  {    
    delete $params{"already"};
    if ( not exists $dependencies{$key}{"source"} )
      { next; }
    if ( $dependencies{$key}{"location"} =~ /atat/ )
      { print OUT "atat/", $key, ".o: atat/", $key, ".cc ";}
    elsif ( $dependencies{$key}{"location"} ne "." )
      { print OUT $dependencies{$key}{"location"}, "/", $key, ".o: ",
                  $dependencies{$key}{"location"}, "/", $key, ".cc ";}
    else
      { print OUT $key, ".o: ", $key, ".cc ";}
    if ( exists $dependencies{$key}{"header"} )
    {
      if ( $dependencies{$key}{"location"} =~ /atat/ )
        { print OUT "atat/", $key, ".h ";}
      elsif ( $dependencies{$key}{"location"} ne "." )
        { print OUT $dependencies{$key}{"location"}, "/", $key, ".h ";}
      else
        { print OUT $key, ".h ";}
    }
    print_recurrent_deps($key);
    print OUT "\n\n";
  }    

  close OUT;

}

sub print_recurrent_deps($)
{
  my $key = $_[0];
  push @{$params{"already"}}, ( $key );
  foreach $dep ( sort { sort_hash($a,$b) } @{$dependencies{$key}{"depends on"}} )
  {
    if ( already_included($dep) )
      { next; }
    if(     ( $dependencies{$key}{"location"} !~ /atat/ )
        and ( $dependencies{$dep}{"location"} =~ /atat/ ) )
      { next; }
    if ( (exists $dependencies{$dep}{"source"}) and
         (exists $dependencies{$dep}{"header"})     )
    {
      print OUT " \\\n\t";
      if ( $dependencies{$dep}{"location"} =~ /atat/ )
        { print OUT "atat/"; }
      elsif ( $dependencies{$dep}{"location"} ne "." )
        { print OUT $dependencies{$dep}{"location"}, "/"; }
      print OUT $dep, ".o";
    }
    elsif ( (exists $dependencies{$dep}{"source"}) )
    {
      print OUT " \\\n\t";
      if ( $dependencies{$dep}{"location"} =~ /atat/ )
        { print OUT "atat/"; }
      elsif ( $dependencies{$dep}{"location"} ne "." )
        { print OUT $dependencies{$dep}{"location"}, "/"; }
      print OUT $dep, ".cc";
    }
    elsif ( exists $dependencies{$dep}{"header"} )
    {
      print OUT " \\\n\t";
      if ( $dependencies{$dep}{"location"} =~ /atat/ )
        { print OUT "atat/"; }
      elsif ( $dependencies{$dep}{"location"} ne "." )
        { print OUT $dependencies{$dep}{"location"}, "/"; }
      print OUT $dep, ".h";
      print_recurrent_deps($dep);
      next;
    }
    push @{$params{"already"}}, ( $dep );
    
  }
}

sub already_included($)
{
  my $is_included = $_[0];
  foreach $dep ( @{$params{"already"}} )
  {
    if ( $dep eq $is_included )
    {  return 1; }
  }
  return 0;
}

sub sort_hash ()
{
  my $a = $_[0], $b = $_[1];

  if(     ($dependencies{$a}{"location"} eq ".")
      and ($dependencies{$b}{"location"} ne ".") )
    { return 0; }
  if(     ($dependencies{$a}{"location"} ne ".")
      and ($dependencies{$b}{"location"} eq ".") )
    { return 1; }
  if( $dependencies{$a}{"location"} ne $dependencies{$b}{"location"} )
    { return ($dependencies{$a}{"location"} > $dependencies{$b}{"location"} ); } 
  if ($a > $b)
    { return 0; }

  return 1;
}

