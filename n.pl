#!/usr/bin/perl

my %params;
open IN, "COMMANDS";
my @COMMANDS;
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
close IN;

$params{"home"} = `cd; pwd`; chomp $params{"home"};
$params{'iaga call'}           =  "lada > out";
$params{'max GA iters'} = 400;

$params{'GA'}{'max generations'} = 200;
$params{'agr'}{'filename'} =  $COMMANDS[0]; chomp $params{'agr'}{'filename'};
$params{'directory'}{'result'} =  $COMMANDS[1]; chomp $params{'directory'}{'result'};
$params{'directory'}{'work'}   =  $COMMANDS[2]; chomp $params{'directory'}{'work'};
$params{'GA style'}            =  $COMMANDS[3]; chomp $params{'GA style'};
$params{'GA style'}            =~ s/^\s+//;
$params{'GA style'}            =~ s/\s+$//;
$params{'minimizer'}           =  $COMMANDS[4]; chomp $params{'minimizer'};
$params{'minimizer'}            =~ s/\s+$//;
$params{'minimizer'}            =~ s/^\s+//;
$params{'CH'}                  =  $COMMANDS[5]; chomp $params{'CH'};
$_ = $COMMANDS[6]; /(\d+)/; $params{'max calls'} = $1;
$_ = $COMMANDS[7]; /popsize:\s+(\d+)\s+replacement:\s+(\S+) gen:\s+(\d+)/;
$params{"GA"}{"population"} = $1;
$params{"GA"}{"replace per generation"} = $2;
$params{"GA"}{"max generations"} = $3;
$_ = $COMMANDS[8]; /(\d+)/;
$params{"every"} = $1;
$params{'taboos'} = $COMMANDS[9]; 
if ( $params{'taboos'} !~ /(age|pop|path|hst)/ )
  { delete $params{'taboos'}; }


$params{'ssc'} = 0;
$params{'nb atoms'} = 20;
if ( $params{'agr'}{'filename'} =~ /atoms:(\d+)_(\d+)_(\d+)/ )
{
  $params{'nb atoms'} = $1*$2*$3;
  $params{'x'} = $1;
  $params{'y'} = $2;
  $params{'z'} = $3;
  $params{'agr'}{'filename'} =~ s/atoms:(\d+)_(\d+)_(\d+)(\_|)//;
  $params{'agr'}{'filename'} =~ s/^\s+//;
  $params{'agr'}{'filename'} =~ s/\s+$//;
  $params{'agr'}{'filename'} = sprintf "atoms:%i_%i_%i_%s",
                                       $params{'x'}, 
                                       $params{'y'}, 
                                       $params{'z'}, 
                                       $params{'agr'}{'filename'};
}
if ( $params{'agr'}{'filename'} =~ /cosa:(\d+)/ )
{
  $params{'ssc'} = $1;
  $params{'agr'}{'filename'} =~ s/cosa:(\d+)(\_|)//;
  $params{'agr'}{'filename'} =~ s/^\s+//;
  $params{'agr'}{'filename'} =~ s/\s+$//;
  $params{'agr'}{'filename'} = sprintf "cosa:%i_%s", $params{'ssc'},
                                       $params{'agr'}{'filename'};
}

if ( $params{'agr'}{'filename'} =~ /full/ )
{
  $params{'ssc'} = $1;
  $params{'agr'}{'filename'} =~ s/full(\_|)//;
  $params{'agr'}{'filename'} =~ s/^\s+//;
  $params{'agr'}{'filename'} =~ s/\s+$//;
  $params{'agr'}{'filename'} = sprintf "full_%s", $params{'ssc'},
                                       $params{'agr'}{'filename'};
  $params{'max GA iters'} = 20;
  $params{'population'} = 50 if ( $params{'nb atoms'} < 12 );
}

if ( $params{'GA style'} =~ /phenotype/ )
{
  $params{'phenotype'} = "true";
  $params{'GA style'} =~ s/phenotype//;
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  $params{'GA style'} = sprintf "pheno_%s", $params{'GA style'};
}




if ( $params{'GA style'} =~ /true/ )
{ 
  $params{'GA style'} =~ s/true//; 
  $params{'GA style'}            =~ s/^\s+//;
  $params{'GA style'}            =~ s/\s+$//;
  $params{'GA style'} = "true_$params{'GA style'}"; chomp $params{'GA style'};
}
if ( $params{'CH'} =~ /one point/i )
{
  $params{'GA style'} = "one_point_$params{'GA style'}";
}
if ( $params{'GA style'} =~ /restart:(\d+)_(\S+)/
     and $params{'GA style'} =~/PopAlgo:(\S+)/ )
{
  @{$params{'restart'}{'PopAlgo'}} = (split /\_/, $1);
  $params{'GA style'} =~ s/PopAlgo:(\S+)//; 
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  $_ = $params{'GA style'}; /restart:(\d+)\_(\S+)/;
  $params{'restart'}{'every'} = $1;
  $params{'restart'}{'rate'} = $2;
  $params{'GA style'} =~ s/restart:(\d+)\_(\S+)//; 
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  my $prefix =  sprintf "restart:%i\_%.2f\_PA:",
                         $params{'restart'}{'every'},
                         $params{'restart'}{'rate'};
  foreach $op ( @{$params{'restart'}{'PopAlgo'}} )
  { 
    if ( $op =~ /UtterRandom/ )
      { $prefix .= "UR_"; }
    elsif ( $op =~ /GradientSA/ )
      { $prefix .= "gSA_"; }
    elsif ( $op =~ /TabooSA/ )
      { $prefix .= "tSA_"; }
    elsif ( $op =~ /linear/ )
      { $prefix .= "linear_"; }
    elsif ( $op =~ /SA/ )
      { $prefix .= "SA_"; }
  } 
  if ( scalar (@{$params{'restart'}{'PopAlgo'}} ) > 0 )
  {
    $params{'GA style'} = sprintf "%s%s", $prefix, $params{'GA style'};
                                 
  }
  else
    { delete $params{'every'}; }

}
else 
  { delete $params{'restart'};  }
if ( $params{'GA style'} =~ /krossover/)
{
  $params{'GA style'} =~ s/krossover//; 
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  $params{'GA style'} =  sprintf "krossover_%s",
                         $params{'GA style'};
}
if ( $params{'GA style'} =~ /range/)
{
  $params{'GA style'} =~ s/range//; 
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  $params{'GA style'} =  sprintf "range_%s",
                         $params{'GA style'};
}
if ( $params{'GA style'} =~ /hst/)
{
  $params{'GA style'} =~ s/hst//; 
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  $params{'GA style'} =  sprintf "hst_%s",
                         $params{'GA style'};
}
if ( $params{'GA style'} =~ /islands:(\d+)/)
{
  my $islands = $1;
  $params{'GA style'} =~ s/islands:(\d+)//; 
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  $params{'GA style'} =  sprintf "islands:%i_%s",
                         $islands,
                         $params{'GA style'};
}

if ( $params{'GA style'} =~ /colonize:(\d+)/)
{
  my $colonize = $1;
  $params{'GA style'} =~ s/colonize:(\d+)//; 
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  $params{'GA style'} =  sprintf "colonize:%i_%s",
                         $colonize,
                         $params{'GA style'};
}

if ( $params{'GA style'} =~ /popgrowth:(\d+)_(\d+)_(\d+)/)
{
  my $every = $1, $maxpop = $2, $rate = $3;
  $params{'GA style'} =~ s/popgrowth:(\d+)_(\d+)_(\d+)//; 
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  $params{'GA style'} =  sprintf "popgrowth:%i_%i_%i_%s",
                         $every, $maxpop, $rate,
                         $params{'GA style'};
}




$params{'agr'}{"filename"} .= sprintf "%s_%s",
                              $params{'GA style'}, 
                              $params{'minimizer'};

if ( $params{'max calls'} != 0 )
{
  $params{'agr'}{"filename"} = "$params{'GA style'}_$params{'minimizer'}_n=$params{'max calls'}";
}

if ( $params{"GA"}{"population"} != 100 )
{
  $params{'agr'}{"filename"} .= "_gen:$params{'GA'}{'population'}";
}
if ( $params{"GA"}{"replace per generation"} != 0.1 )
{
  $params{'agr'}{"filename"} .= "_rep:$params{'GA'}{'replace per generation'}";
}
if ( exists $params{"minimize best"} )
{
  $params{'agr'}{"filename"} = sprintf "mb:%.3f_$params{'agr'}{'filename'}",
				       $params{'every'};
}

if ( exists $params{'taboos'} )
{
  if( $params{'taboos'} =~ /age/ )
  {
    $params{'agr'}{"filename"} = sprintf "%s_%s", "age",
                                 $params{'agr'}{'filename'};
  } 
  if( $params{'taboos'} =~ /pop/ )
  {
    $params{'agr'}{"filename"} = sprintf "%s_%s", "pop",
                                 $params{'agr'}{'filename'};
  } 
  if( $params{'taboos'} =~ /path/ )
  {
    $params{'agr'}{"filename"} = sprintf "%s_%s", "path",
                                 $params{'agr'}{'filename'};
  } 
  if( $params{'taboos'} =~ /hst/ )
  {
    if ( $params{'GA style'} !~ /hst/ )
    {
      $params{'taboos'} =~ s/hst//;
      if ( $params{'taboos'} !~ /(age|pop|path)/ )
        { delete $params{'taboo'}; } 
    }
    $params{'agr'}{"filename"} = sprintf "%s_%s", "thst",
                                 $params{'agr'}{'filename'};
  } 
}

$params{'xml'}{"filename"} = $params{'agr'}{'filename'};
$params{'agr'}{"filename"} .= ".agr";
$params{'xml'}{"filename"} .= ".xml";


# begin work
my %structure;

my $cosa = 0;
system "rm -f $params{'agr'}{'filename'}";
read_structure(); # data passes from PIfile to %structure hash
launch_iaga(); 


exit;


sub launch_iaga()
{
  for( my $count =0; $count < $params{'max GA iters'}; $count++)
  {
    write_lamarck_input();
    if ( $params{'GA style'} !~ /true/ and
         $params{'agr'}{'filename'} !~ /full/ )
      { system "rm -f convex_hull.xml"; }
    system "$params{'iaga call'}";
    system "echo '# new GA ' >> $params{'agr'}{'filename'}";
    system "cat convex_hull.agr >> $params{'agr'}{'filename'}";
    system "cp convex_hull.xml $params{'xml'}{'filename'}  ";
  }
  system "cp $params{'xml'}{'filename'} $params{'agr'}{'filename'}  $params{'directory'}{'result'} ";
}

sub read_structure()
{
  my $n = $params{'nb atoms'}; 
  @{$structure{'cell'}} = ( [ 0.0, 0.5*$params{'x'}, 0.5*$params{'x'} ],
                            [ 0.5*$params{'y'}, 0.0, 0.5*$params{'y'} ],
                            [ 0.5*$params{'z'}, 0.5*$params{'z'}, 0.0 ] );
  @{$structure{'atomic positions'}}  = ();
  
  my $index = 0;
  for($i=0;$i<$params{'x'};$i++)
  { 
    for($j=0;$j<$params{'y'};$j++)
    {
      for($k=0;$k<$params{'z'};$k++)
      {
        $structure{'atomic positions'}[$index] = 0.5*($j+$k);
        $index++;
        $structure{'atomic positions'}[$index] = 0.5*($i+$k);
        $index++;
        $structure{'atomic positions'}[$index] = 0.5*($i+$j);
        $index++;
      }
    }
  }

}


sub write_lamarck_input()
{
  open IN, "<input_template.xml"
    or die "No input_template.xml file found";
  open OUT, ">input.xml"
    or die "Could not open input.xml";

  while( <IN> )
  {

    if( /\?structure\?/ )
    {
      printf OUT "  <Structure>\n";
      printf OUT "    <Cell>\n";
      printf OUT "      <column x=\"%.6f\" y=\"%.6f\" z=\"%.6f\" />\n",
                 @{$structure{'cell'}[0]}[0..2];
      printf OUT "      <column x=\"%.6f\" y=\"%.6f\" z=\"%.6f\" />\n",
                 @{$structure{'cell'}[1]}[0..2];
      printf OUT "      <column x=\"%.6f\" y=\"%.6f\" z=\"%.6f\" />\n",
                 @{$structure{'cell'}[2]}[0..2];
      printf OUT "    </Cell>\n";

      for( my $i = 0; $i < scalar( @{$structure{'atomic positions'}} ); 
           $i +=3 )
      {
        printf OUT "    <Atom x=\"%.6f\" y=\"%.6f\" z=\"%.6f\" type=\"%.6f\" />\n",
                 $structure{'atomic positions'}[$i], 
                 $structure{'atomic positions'}[$i+1], 
                 $structure{'atomic positions'}[$i+2]; 
      }
      printf OUT "  </Structure>\n";
    }
    elsif( /\?GA\?/ )
    {
      printf OUT "  <GA";
      if ( $params{'GA style'} =~ /lamarck_evolution/i )
        { printf OUT " method=\"lamarck evolve from start\""; }
      elsif ( $params{'GA style'} =~ /lamarck/i )
        { printf OUT " method=\"lamarck\""; }
      elsif ( $params{'GA style'} =~ /multistart/i )
        { printf OUT " method=\"lamarck multistart evolve from start\""; }
      elsif ( $params{'GA style'} =~ /darwin/i )
        { printf OUT " method=\"darwin\""; }
      if ( $params{'GA style'} =~ /islands:(\d+)_/i )
        { printf OUT " islands=\"$1\""; }
      printf OUT " maxgen=\"%i\">\n",
                 $params{'GA'}{'max generations'};

      if ( exists $params{'phenotype'} )
        { printf OUT "    <Phenotype/> \n"; }

      my $minizertype;
      if ( $params{'minimizer'} =~ /linear/i )
      {
       $minimizertype = sprintf "<Minimizer type=\"linear\" maxeval=\"%i\ notaboo=\"true\" "; 
                        $params{'max calls'};
      }
      elsif ( $params{'minimizer'} =~ /GradientSA/ ) 
      {
        $minimizertype = sprintf "<Minimizer type=\"GradientSA\" maxeval=%i ",
                                 $params{'max calls'}; 
      }
      elsif ( $params{'minimizer'} =~ /TabooSA/ ) 
      {
        $minimizertype = sprintf "<Minimizer type=\"SA\" maxeval=%i ",
                                 $params{'max calls'}; 
      }
      elsif ( $params{'minimizer'} =~ /sa/i )
      {
       $minimizertype = sprintf "<Minimizer type=\"SA\" maxeval=\"%i\" notaboo=\"true\" ",
                                $params{'max calls'};
      }

      # creates the operators
      if ( $params{'GA style'} =~ /multistart/i )
      {
        printf OUT "    <Operators type=\"and\" >\n";
        printf OUT "      <UtterRandom prob=1.0 />\n";
        printf OUT "      %s prob=1.0 />\n", $minimizertype;
        printf OUT "    </Operators>\n";
      }
      elsif($params{'GA style'} =~ /lamarck/i ) 
      {
        printf OUT "      <Operators type=\"and\" >\n";
        printf OUT "        <TabooOp> \n";
        if ( $params{'GA style'} =~ /krossover/ and  $params{'GA style'} =~ /range/) 
          { printf OUT "          <Krossover type=\"range\" />\n"; }
        elsif ( $params{'GA style'} =~ /krossover/ )
          { printf OUT "          <Krossover/>\n"; }
        else
          { printf OUT "          <Crossover />\n"; }
        printf OUT "        </TabooOp> \n";
        printf OUT "        %s />\n", $minimizertype;
        printf OUT "      </Operators>\n";
      }
      elsif($params{'GA style'} =~ /darwin/i ) 
      {
        printf OUT "    <TabooOp> \n";
        if ( $params{'GA style'} =~ /krossover/ and  $params{'GA style'} =~ /range/) 
          { printf OUT "          <Krossover type=\"range\" />\n"; }
        elsif ( $params{'GA style'} =~ /krossover/ )
          { printf OUT "          <Krossover/>\n"; }
        else
          { printf OUT "      <Crossover value=\"0.5\"/>\n"; }
        printf OUT "    </TabooOp> \n";
      }


      printf OUT "    <Selection size=2 />\n";

      printf OUT "    <Population size=\"%i\"/>\n",
                 $params{'GA'}{'population'};
      printf OUT "    <Offsprings rate=\"%.4f\"/>\n",
                 $params{'GA'}{'replace per generation'};
      if( $params{'CH'} =~ /one point/i )
        { printf OUT "    <OnePointHull/>\n"; }
      if( exists $params{'restart'}  )
      {
        printf OUT "    <PopAlgo every=\"%i\" rate=\"%.2f\" >\n", 
                   $params{"restart"}{'every'},
                   $params{"restart"}{'rate'};
        printf OUT "      <Operators type=\"and\">\n"; 
        foreach $op ( @{$params{'restart'}{'PopAlgo'}}  )
        {
          if ( $op =~ /UtterRandom/ )
            { printf OUT "        <UtterRandom/>\n"; }
          elsif ( $op =~ /GradientSA/ )
            { printf OUT "        <TabooMinimizer type=\"GradientSA\"/>\n"; }
          elsif ( $op =~ /TabooSA/ )
            { printf OUT "        <TabooMinimizer type=\"TabooSA\"/>\n"; }
          elsif ( $op =~ /linear/ )
            { printf OUT "        <Minimizer type=\"linear\"/>\n"; }
          elsif ( $op =~ /SA/ )
            { printf OUT "        <Minimizer type=\"sa\"/>\n"; }
        }
        printf OUT "      </Operators>\n";
        printf OUT "    </PopAlgo>\n";
      }
      if( exists $params{'taboos'} )
#         and scalar( @{$structure{'atomic positions'}} ) >= 36 )
      {
        printf OUT "    <Taboos>\n"; 
        if ( $params{'taboos'} =~ /pop/ )
          { printf OUT "      <PopTaboo/>\n      <OffspringTaboo/>\n";  }
        if ( $params{'taboos'} =~ /path/ )
          { printf OUT "      <PathTaboo/>\n";  }
        if ( $params{'taboos'} =~ /hst/ )
        {
          printf OUT "      <HistoryTaboo/> \n";  
        }
        elsif ( $params{'taboos'} =~ /age/ )
        {
          printf OUT "      <AgeTaboo lifespan=%i />\n",
                     $params{'every'};  
          printf OUT "      <NuclearWinter >\n";  
          printf OUT "        <Mutation value=\"0.5\" />\n";  
          printf OUT "      </NuclearWinter>\n";  
        }
        printf OUT "    </Taboos>\n"; 
      }
      if( $params{'GA style'} =~ /colonize:(\d+)/ )
      {  
        printf OUT "    <Colonize every=\"%i\">\n", $1;
          printf OUT "    <TabooOp>\n";
          printf OUT "      <Krossover/>\n";
          printf OUT "    </TabooOp>\n";
        printf OUT "    </Colonize>\n";
      }
      if( $params{'GA style'} =~ /popgrowth:(\d+)_(\d+)_(\d+)/ )
      {  
        printf OUT "    <PopGrowth every=\"%i\" maxpop=\"%i\" rate=\"%i\" >\n", $1, $2, $3;
          printf OUT "    <TabooOp>\n"; 
          printf OUT "      <Krossover/>\n";
          printf OUT "    </TabooOp>\n";
        printf OUT "    </PopGrowth>\n";
      }
      if( $params{'GA style'} =~ /hst/ )
        { printf OUT "    <History/>\n"; }


      if ( $params{"minimizer"} =~ /GradientSA/ )
      {
        if ( $params{"nb atoms"} > 20 )
           { printf OUT "    <Terminator ref=\"gradient\" value=250000/>\n"; }
        else
           { printf OUT "    <Terminator ref=\"gradient\" value=50000/>\n"; }
      }
      else
      {
        if ( $params{"nb atoms"} > 20 )
           { printf OUT "    <Terminator ref=\"evaluation\" value=25000/>\n"; }
        else
           { printf OUT "    <Terminator ref=\"evaluation\" value=5000/>\n"; }
      }

#     printf OUT "    <tSeed n=\"5\"/>\n";
      printf OUT "    <Statistics type=\"diversity\"/>\n";
      printf OUT "    <Statistics type=\"true census\"/>\n";
      printf OUT "  </GA>\n";
    }
    else
    {
      print OUT $_;
    }
  }
  

  close IN;
  close OUT;
}

