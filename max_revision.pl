#! /usr/bin/perl

my $max_revision = -1;
my $min_revision = -1;
my $cur_revision = -1;
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
  printf OUT "namespace SVN \n{\n  const types::t_unsigned Revision = %i; \n}\n\n", $max_revision;
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
      if( /Version:+(\s+|)\$Id:\s+(\S+)\s+(\d+)/ )
      {
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

