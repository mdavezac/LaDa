#! /usr/bin/perl
use Getopt::Long;

exit if( not -e "../../darwin/ce_opt" );
my %options = ();
$options{"reruns"} = 1;
GetOptions (\%options, 'mpirun|n=i', 'reruns=i', 'make', 'kdbg' );
die "Cannot perform $options{'reruns'} runs.\n" if ( $options{"reruns"} <= 0 );

my $cmdl = "./ce_opt -i ../input.xml";


chdir "calc";
system "rm krossover_VA_gen:20.agr" if ( -e "krossover_VA_gen:20.agr" );
system "ln -s ../../../darwin/ce_opt " unless ( -e "ce_opt" );

system "cd ../../../; make" if ( exists $options{"make"} );

if( exists $options{"kdbg"} ) 
{
  system "ln -s ../input.xml" if( not -e "input.xml" );
  system "kdbg ce_opt";
  exit;
}

if( exists $options{"mpirun"} )
{
  die "Cannot run on $options{'mpirun'} procs.\n" if ( $options{"mpirun"} <= 0 );
  $cmdl = sprintf "mpirun -n %i %s", $options{"mpirun"}, $cmdl;
} 
$cmdl = sprintf "%s --reruns %i ", $cmdl, $options{"reruns"} unless( $options{"reruns"} == 1 );
print "$cmdl\n";
system "$cmdl";
system "rm ce_opt" if ( -e "ce_opt" );
chdir "../";

