#! /usr/bin/perl

my $max_revision = -1;
my $min_revision = -1;
my $cur_revision = -1;
my $year;
my $month;
my $day;
my $hour;
my $minute;
my $second;
my $user;
my @alldires;

# get_all_dirs();
get_current_revision();
get_max_revision();
dodir();


if ( $cur_revision < $max_revision )
{
  open OUT, ">revision.h";

  printf OUT "#ifndef _REVISION_H_\n#define _REVISION_H_\n\n";
  printf OUT "#include <opt/types.h> \n\n";
  printf OUT "//! \\brief Contains SVN related data.\n";
  printf OUT "//! \\details Everything here is determined at compile time by\n";
  printf OUT "//!          the script max_revision.pl\n";
  printf OUT "//! \\warning Don't trust the documentation, get your numbers \n";
  printf OUT "//!          directly from revision.h. Or rerun ./max_revision.pl \n";
  printf OUT "//!          followed by doxygen doxy.conf.\n";
  printf OUT "namespace SVN \n{\n";
  printf OUT "  //! Revision number of this version of LaDa. \n";
  printf OUT "  const types::t_unsigned Revision = %i; \n", $max_revision;
  printf OUT "  //! Year of last revision. \n";
  printf OUT "  const types::t_unsigned Year = %i; \n", $year;
  printf OUT "  //! Month of last revision. \n";
  printf OUT "  const types::t_unsigned Month = %i; \n", $month;
  printf OUT "  //! Day of last revision. \n";
  printf OUT "  const types::t_unsigned Day = %i; \n", $day;
  printf OUT "  //! Last revision committed by user %s \n", $user;
  printf OUT "  const std::string User = \"%s\"; \n", $user;
  printf OUT "}\n\n";
  printf OUT "#endif\n";

  close OUT;
}

sub get_current_revision()
{
  open IN, "revision.h" or return;
  while ( $_=<IN> )
  {
    next if( not /types::t_unsigned Revision = (\d+)/ );
    
    $cur_revision = $1;
    last;
  }
}

sub dodir()
{
  foreach $dir ( <*/> )
  {
    next if (not -d $dir );
    next if ($dir =~ /LaDa/);
    chdir $dir;
    get_max_revision();
    chdir "..";
  }
}

sub get_max_revision()
{
  foreach $file ( <*.h>, <*.cc> ) 
  {
    next if ( $file =~ /config\.h/ );

    open IN, "$file";
    my $n = 0;
    my $found = 0;
    while ( ($_=<IN>) and ($n < 5) )
    {
      if( /Version:+(\s+|)\$Id:\s+(\S+)\s+(\d+)\s+(\d+)-(\d+)-(\d\d) (\d\d):(\d\d):(\d\d)\S (\S+)/ )
      {
        $year = $4 if ( $3 > $max_revision );
        $month = $5 if ( $3 > $max_revision );
        $day = $6 if ( $3 > $max_revision );
        $hour = $7 if ( $3 > $max_revision );
        $minute = $8 if ( $3 > $max_revision );
        $second = $9 if ( $3 > $max_revision );
        $user = $10 if ( $3 > $max_revision );
        $max_revision = $3 if ( $3 > $max_revision );
        $min_revision = $3 if ( $min_revision == -1 );
        $min_revision = $3 if ( $3 < $min_revision and $3 != -1);
        $found = 1;
        last;
      }
      $n ++;
    }
    close IN; 
    print "Revision String not found in $file\n" if ($found != 1 and $file !~ /revision.h/ );
  }
}

