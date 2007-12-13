#!/usr/bin/perl
#
#  Version: $Id$
#

# 1.  run type (multistart|lamarck + true )
# 2.  miniizer type (none|sa|linear|wang|ssquare)
# 3.  CH (one point|all points)
# 4.  max calls
# 5.  pop size
# 6.  replacement rate
# 7.  minimize best rate

my @job = (#["popsize:20", "CH", "Krossover VA",         "pop", 400, "" ],
#           ["popsize:60", "CH", "Krossover VA",         "pop", 400, "" ], 
#           ["popsize:100","CH", "Krossover VA",         "pop", 400, "" ], 
            ["popsize:20 rate:0.1", "CH", "Krossover",           "pop", 400, "" ],
            ["popsize:20 rate:0.2", "CH", "Krossover",           "pop", 400, "" ],
            ["popsize:20 rate:0.4", "CH", "Krossover",           "pop", 400, "" ],
            ["popsize:20 rate:0.6", "CH", "Krossover",           "pop", 400, "" ],
            ["popsize:20 rate:0.8", "CH", "Krossover",           "pop", 400, "" ] );
#           ["popsize:60", "CH", "Krossover",            "pop", 400, "" ], 
#           ["popsize:100","CH", "Krossover",            "pop", 400, "" ], 
#           ["popsize:1 rate=1",  "CH", "UtterRandom VA",       "",    400, "" ],
#           ["popsize:1 rate=1",  "CH", "UtterRandom Beratan",  "",    400, "" ] );
#           ["popsize:60", "CH", "Krossover VA",         "pop", 400, " true " ] );
my @prefixes = ("atoms:2_2_3");
# 

foreach $prefix ( @prefixes )
{
  for( my $i = 0; $i < scalar(@job) ; $i++)
  {
    system "rm -f sendme_$prefix$i ";
    open OUT, ">sendme_$prefix$i" or die;
    my $work = "work_$i";
system "mkdir $work";
    
    printf OUT "#! /bin/tcsh\n#\n";
    printf OUT "#PBS -l nodes=1:ppn=1,walltime=24:00:0\n";
    printf OUT "#PBS -q Std\n";
#   printf OUT "#PBS -l ncpus=1,walltime=48:00:0\n";
#   printf OUT "#PBS -q single\n";
    printf OUT "#PBS -m n \n";
    printf OUT "#PBS -e std.err\n";
    printf OUT "#PBS -o std.out\n\n";
    printf OUT "set WORK_DIR=(\$PBS_O_WORKDIR/$work)\n";
    printf OUT "set LOAD_DIR=(\$PBS_O_WORKDIR)\n";
    printf OUT "set RESULT_DIR=(\$PBS_O_WORKDIR/results)\n\n";
    printf OUT "module load intel\n\n";
    printf OUT "cd \"\$LOAD_DIR\"\n";
    printf OUT "rm sendme_$prefix$i.*\n";
    printf OUT "cp ../lada ../n.pl input_template.xml \"\$WORK_DIR\"\n\n";
    if ( $job[$i][0] =~ /local/ )
    {
      $_ = $prefixes[$i]; /atoms:(\d+)_(\d+)_(\d+)/;
      my $chname = sprintf "true_ch_%i_%i_%i.xml", $1, $2, $3;
      printf OUT "cp $chname \"\$WORK_DIR\"\n";
    }
    printf OUT "cd \"\$WORK_DIR\"\n\n";
    printf OUT "cat >COMMANDS<<EOF\n";
    printf OUT "%s\n",$prefix;
    printf OUT "\$RESULT_DIR\n";
    printf OUT "%s\n%s\n%s\n%s\n%s\n%s\n",
  	 $job[$i][0], $job[$i][1], $job[$i][2], $job[$i][3],
  	 $job[$i][4], $job[$i][5];
    printf OUT "EOF\n";
    printf OUT "chmod u+x n.pl\n\n";
    printf OUT "n.pl > status\n\n";
    printf OUT "exit\n";
  
    close OUT;
  
    system "chmode u+x sendme_$prefix$i";
  
    system "qsub sendme_$prefix$i > $prefix$i\_present_job"; 
  
  } 
  
  system "chmod u-x send_all.pl";

}
