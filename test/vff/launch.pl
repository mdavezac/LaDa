#! /usr/bin/perl


my $argument = shift;
my $prefix = "";
my $program = "vff";
my $program2 = "layered_vff";
my $directory = "../../vff";
$prefix = "valgrind -v --leak-check=full --show-reachable=yes "  if( $argument  =~ /valgrind/ );
system "qsub sendme > jobid"  if( $argument  =~ /pbs/ );
system "qdel `cat jobid`; rm jobid"  if( $argument  =~ /kill/  and -e "jobid");
system "cd ../../; make"  if( $argument  =~ /make/ );

if( not -e "$directory/$program" )
{
  print "$directory/$program was not found.\n";
  exit;
}
if( not -e "$directory/$program2" )
{
  print "$directory/$program2 was not found.\n";
  exit;
}

system "rm -rf calc/*";

chdir "calc";
system "ln -s ../$directory/$program " if ( not -e "$program" );
system "ln -s ../$directory/$program2 " if ( not -e "$program2" );

system "$prefix ./$program -i ../quaternary.xml | grep \"Energy\" ";
system "$prefix ./$program -i ../ternary.xml | grep \"Energy\" ";
system "$prefix ./$program2 -i ../quaternary.xml | grep \"Energy\" ";
system "$prefix ./$program2 -i ../ternary.xml | grep \"Energy\" ";

chdir "..";

