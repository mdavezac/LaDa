#!/usr/bin/perl

# 1.  run type (multistart|lamarck + true )
# 2.  miniizer type (none|sa|linear|wang|ssquare)
# 3.  CH (one point|all points)
# 4.  max calls
# 5.  pop size
# 6.  replacement rate
# 7.  minimize best rate

my @job = ( #["darwin krossover true", "none",       "none", 0, 50, 0.1, 10000, 0.0, "pop"] );
            #["darwin",           "none",       "none", 0, 200, 0.1, 10000, 0.0, "pop"],
             ["darwin krossover full:0-0 ", "none",       "none", 0,  200,  0.1, 10000000, 0.0, "pop"] );
            #["lamarck krossover terminator:3000", "GradientSA",       "one point", 0, 10,  0.1, 10000000, 0.0, "pop"],
            #["lamarck krossover terminator:3000", "GradientSA",       "one point", 0, 20,  0.1, 10000000, 0.0, "pop"],
            #["lamarck krossover terminator:3000", "GradientSA",       "one point", 0, 40,  0.1, 10000000, 0.0, "pop"],
            #["lamarck krossover terminator:3000", "GradientSA",       "one point", 0, 60,  0.1, 10000000, 0.0, "pop"],
            #["lamarck krossover terminator:3000", "GradientSA",       "one point", 0, 80,  0.1, 10000000, 0.0, "pop"],
            #["darwin krossover terminator:4000", "none",       "one point", 0, 20,  0.1, 10000000, 0.0, "pop"],
            #["darwin krossover terminator:4000", "none",       "one point", 0, 40,  0.1, 10000000, 0.0, "pop"],
            #["darwin krossover terminator:4000", "none",       "one point", 0, 80,  0.1, 10000000, 0.0, "pop"] );
            #["multistart terminator:200000",  "GradientSA", "one point", 0, 100, 0.1, 10000000, 0.0, "pop"],
            #["multistart terminator:200000",  "TabooSA",    "one point", 0, 100, 0.1, 10000000, 0.0, "pop"] );
            #["multistart true",       "GradientSA", "none", 0, 200, 0.1, 10000, 0.0, "pop"] );:w
            #["darwin",           "none",       "one point", 0, 200, 0.1, 10000, 0.0, "pop"],
            #["darwin krossover", "none",       "one point", 0, 60, 0.1, 10000, 0.0, "pop"],
            #["darwin krossover", "none",       "one point", 0, 70, 0.1, 10000, 0.0, "pop"],
            #["darwin krossover", "none",       "one point", 0, 80, 0.1, 10000, 0.0, "pop"],
            #["darwin krossover", "none",       "one point", 0, 90, 0.1, 10000, 0.0, "pop"],
            #["darwin krossover", "none",       "one point", 0, 100, 0.1, 10000, 0.0, "pop"] );
            #["multistart",       "GradientSA", "one point", 0, 200, 0.1, 10000, 0.0, "pop"] );
my @prefixes = ("full");#, "atoms:2_2_4" );
# 

foreach $prefix ( @prefixes )
{
  my $work = "work_$prefix";
  
  for( my $i = 0; $i < scalar(@job) ; $i++)
  {
    system "rm -f sendme_$prefix$i $work$i/present_job";
    if ( not -e "$work$i" ) 
      { system "mkdir $work$i"; }
    else
      { system "rm -f $work$i/*"; }
    open OUT, ">sendme_$prefix$i" or die;
    
    printf OUT "#! /bin/tcsh\n#\n";
    printf OUT "#PBS -l ncpus=1,walltime=48:00:0\n";
    printf OUT "#PBS -q single\n";
    printf OUT "#PBS -m n \n";
  # printf OUT "#PBS -e $work$i/std.err\n";
  # printf OUT "#PBS -o $work$i/std.out\n\n";
    printf OUT "set WORK_DIR=(\$PBS_O_WORKDIR/$work$i)\n";
    printf OUT "set LOAD_DIR=(\$PBS_O_WORKDIR)\n";
    printf OUT "set RESULT_DIR=(\$PBS_O_WORKDIR/results)\n\n";
    printf OUT "module load intel\n\n";
    printf OUT "cd \"\$LOAD_DIR\"\n";
    printf OUT "rm sendme_$prefix$i.*\n";
    printf OUT "cp ../lada ../vadsds.pl input_template.xml \"\$WORK_DIR\"\n\n";
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
    printf OUT "\$WORK_DIR\n";
    printf OUT "%s\n%s\n%s\n%i\n",
  	 $job[$i][0], $job[$i][1], $job[$i][2], $job[$i][3];
    printf OUT "popsize: %i  replacement: %f gen: %i \n", 
  	     $job[$i][4], $job[$i][5], $job[$i][6];
    printf OUT "minimize best: %.4f  every: %i \n", 
  	     $job[$i][7], 1000;
    printf OUT "taboos: %s \n", 
  	     $job[$i][8];
    printf OUT "EOF\n";
    printf OUT "chmod u+x vadsds.pl\n\n";
    printf OUT "vadsds.pl > status\n\n";
    printf OUT "exit\n";
  
    close OUT;
  
    system "chmode u+x sendme_$prefix$i";
  
    system "qsub sendme_$prefix$i >> $work$i/present_job"; 
  
  } 
  
  system "chmod u-x send_full.pl";

}
