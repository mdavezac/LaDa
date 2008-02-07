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

chdir "calc";
system "ln -s ../../../pescan_interface/escan " if ( not -e "escan" );

print "$prefix ./escan -i ../input.xml >& stdout \n";

system "bandgap_opt" if ( -e "escan" );
chdir "..";

