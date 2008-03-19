#! /usr/bin/perl


if( not -e "../../pescan_interface/escan" )
{
  print "../../pescan_interface/escan was not found.\n";
  exit;
}
system "rm -rf calc/*";

my $argument = shift;
my $np = shift;
my $prefix = "";
$prefix = "valgrind -v --leak-check=full --show-reachable=yes "  if( $argument  =~ /valgrind/ );
$prefix = "mpirun -n 2 "  if( $argument  =~ /mpirun/ );

chdir "calc";
system "ln -s ../../../pescan_interface/escan " if ( not -e "escan" );

print "$prefix ./escan -i ../input.xml >& stdout \n";

system "$prefix ./escan -i ../input.xml >& stdout" if ( -e "escan" );
chdir "..";

