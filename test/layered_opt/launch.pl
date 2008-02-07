#! /usr/bin/perl


my $argument = shift;
my $np = shift;
my $prefix = "";
$prefix = "valgrind -v --leak-check=full --show-reachable=yes "  if( $argument  =~ /valgrind/ );
$prefix = "mpirun -n $np "  if( $argument  =~ /mpirun/ and $np > 1);
$prefix = "kdbg "  if( $argument  =~ /kdbg/ );
$prefix = "mpirun ", shift if( $argument =~/mpipbs/ );
system "qsub sendme > jobid"  if( $argument  =~ /\bpbs/ );
system "qdel `cat jobid`; rm jobid"  if( $argument  =~ /kill/  and -e "jobid");
system "cd ../../; make"  if( $argument  =~ /make/ );
exit if( $argument  =~ /(\bpbs|kill|make)/ );

if( not -e "../../darwin/layered_opt" )
{
  print "../../darwin/layered_opt was not found.\n";
  exit;
}

system "rm -rf calc/*";

chdir "calc";
system "ln -s ../../../darwin/layered_opt ." if ( not -e "layered_opt" );
system "ln -s ../input.xml";

print "$prefix ./layered_opt >& stdout\n" if( $argument !~ /kdbg/ );
system "$prefix ./layered_opt >& stdout" if( $argument !~ /kdbg/ );
system "$prefix ./layered_opt" if( $argument =~ /kdbg/ );

chdir "..";

