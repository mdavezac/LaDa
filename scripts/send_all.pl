#! /usr/bin/perl
#
#  Version: $Id$
#

my @job = ( 
  ["crossover mutation:0.2_0.02 terminator:10000 CH popsize:200 poptaboo ", 6],
  ["krossover terminator:30000 CH popsize:800 poptaboo ", 6] );

for( my $i = 0; $i < scalar(@job) ; $i++)
{
  my $prefix = sprintf "atoms:2_2_%i", $job[$i][1];
  my $work = "work_$prefix";

  system "rm -f sendme_$prefix%i",$i; # $work$i/present_job";
  if ( not -e "$work$i" ) 
    { system "mkdir $work$i"; }
  else
    { system "rm -f $work$i/*"; }
  open OUT, ">sendme_$prefix$i" or die;
  
  printf OUT "#! /bin/tcsh\n#\n";
  printf OUT "#PBS -l ncpus=1,walltime=05:00:0\n";
  printf OUT "#PBS -q Std\n";
  printf OUT "#PBS -m n \n";
# printf OUT "#PBS -e $work$i/std.err\n";
# printf OUT "#PBS -o $work$i/std.out\n\n";
  printf OUT "set WORK_DIR=(\$PBS_O_WORKDIR/$work$i)\n";
  printf OUT "set LOAD_DIR=(\$PBS_O_WORKDIR)\n";
  printf OUT "set RESULT_DIR=(\$PBS_O_WORKDIR/results)\n\n";
  printf OUT "module load intel\n\n";
  printf OUT "cd \"\$LOAD_DIR\"\n";
  printf OUT "rm sendme_$prefix$i.*\n";
  printf OUT "cp ../darwin/.libs/ce_opt ../scripts/nxnxn.pl ../scripts/input_template.xml \"\$WORK_DIR\"\n\n";
  printf OUT "cd \"\$WORK_DIR\"\n\n";
  printf OUT "cat >COMMANDS<<EOF\n";
  printf OUT "%s\n",$prefix;
  printf OUT "\$RESULT_DIR\n";
  printf OUT "\$WORK_DIR\n";
  printf OUT "%s\n",
	 $job[$i][0];
  printf OUT "EOF\n";
  printf OUT "chmod u+x nxnxn.pl\n\n";
  printf OUT "./nxnxn.pl > status\n\n";
  printf OUT "./ce_opt input.xml -r 400 \n";
  printf OUT "exit\n";

  close OUT;

  system "chmode u+x sendme_$prefix$i";

  system "qsub sendme_$prefix$i >> $work$i/present_job"; 

} 

system "chmod u-x send_all.pl";

