#! /usr/bin/perl


my $argument = shift;
my $np = shift;
my $prefix = "";
$prefix = "valgrind -v --leak-check=full --show-reachable=yes "  if( $argument  =~ /valgrind/ );
$prefix = "mpirun -n $np "  if( $argument  =~ /mpirun/ and $np > 1);
system "qsub sendme > jobid"  if( $argument  =~ /pbs/ );
system "qdel `cat jobid`; rm jobid"  if( $argument  =~ /kill/  and -e "jobid");
system "cd ../../; make"  if( $argument  =~ /make/ );
exit if( $argument  =~ /(pbs|kill|make)/ );

if( not -e "../../darwin/bandgap_opt" )
{
  print "../../darwin/bandgap_opt was not found.\n";
  exit;
}

system "rm -rf calc/*";

chdir "calc";
system "ln -s ../../../darwin/bandgap_opt " if ( not -e "bandgap_opt" );

print "$prefix ./bandgap_opt -i ../input.xml >& stdout \n";
system "$prefix ./bandgap_opt -i ../input.xml>& stdout";

chdir "..";

