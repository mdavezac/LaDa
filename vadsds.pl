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
if ( $params{'taboos'} !~ /(age|pop|unique)/ )
  { delete $params{'taboos'}; }

if ( $params{'GA style'} =~ /PopAlgo/ )
{
  $params{'GA style'} =~ s/PopAlgo//;
}
else 
{ delete $params{'minimize best'}; }
$params{'iaga call'}           =  "lada > out";

$params{'nb atoms'} = 20;
$params{"file"}{"Pi"} = "$params{'home'}/nanopse/cell_shapes/fcc_7-32";
$params{'max GA iters'} = 200;

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
if ( $params{'GA style'} =~ /restart/)
{
  $params{'GA style'} =~ s/restart//; 
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  $params{'GA style'} =  sprintf "restart:%i\_%s",
                         $params{'every'},
                         $params{'GA style'};
}
if ( $params{'GA style'} =~ /krossover/)
{
 print $params{'GA style'}, "\n";
  $params{'GA style'} =~ s/krossover//; 
 print $params{'GA style'}, "\n";
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  $params{'GA style'} =  sprintf "krossover_%s",
                         $params{'GA style'};
 print $params{'GA style'}, "\n";
}
if ( $params{'GA style'} =~ /best:(\d+\.\d+)/)
{
  $params{'rate'} = $1;
  $params{'GA style'} =~ s/best:(\d+\.\d+)//; 
  $params{'GA style'} =~ s/^\s+//;
  $params{'GA style'} =~ s/\s+$//;
  $params{'GA style'} =  sprintf "best:%.2f\_every:%i_%s",
                         $params{'rate'},
                         $params{'every'},
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




$params{'agr'}{"filename"} .= sprintf "%s_%s",
                              $params{'GA style'}, 
                              $params{'minimizer'};
if ( $params{'nb atoms'} != 20 )
{
  $params{'agr'}{"filename"} = sprintf "%i_%s",
                               $params{'nb atoms'},
                               $params{'agr'}{'filename'};
}

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
}

$params{'xml'}{"filename"} = $params{'agr'}{'filename'};
$params{'agr'}{"filename"} .= ".agr";
$params{'xml'}{"filename"} .= ".xml";


# begin work
my %structure;

open (PI, "$params{'file'}{'Pi'}") 
  or die "Did not find PI file $params{'file'}{'Pi'}";

while(<PI>)
{
  if(/NT/)
  {
    read_structure(); # data passes from PIfile to %structure hash
	
    if ( scalar( @{$structure{'atomic positions'}} ) < $params{'nb atoms'}*3 ) 
      { next; }
    else 
      { launch_iaga(); };
  exit;

  }
}

close PI;



exit;


sub launch_iaga()
{
  system "rm -f $params{'agr'}{'filename'}";
  for( my $count =0; $count < $params{'max GA iters'}; $count++)
  {
    write_lamarck_input();
    if ( $params{'GA style'} !~ /true/ )
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
  my $dummy;
  ($dummy,$structure{'nb sites'},$structure{'nb triads'},
   @{$structure{'cell'}[0]}[1..3],@{$structure{'cell'}[1]}[1..3],
   @{$structure{'cell'}[2]}[1..3]) = (split);
  
  for($i=1;$i<=3;$i++)
  { 
    $structure{'cell'}[0][$i] /= 2;
    $structure{'cell'}[1][$i] /= 2;
    $structure{'cell'}[2][$i] /= 2;
  }

  # Skip the PIs
  while( ($_=<PI>) and ($_ !~ /BASIS/) ){}

  # atoms
  s/BASIS//; # Get the basis atom positions
  @{$structure{'atomic positions'}} = (split);
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
                 @{$structure{'cell'}[0]}[1..3];
      printf OUT "      <column x=\"%.6f\" y=\"%.6f\" z=\"%.6f\" />\n",
                 @{$structure{'cell'}[1]}[1..3];
      printf OUT "      <column x=\"%.6f\" y=\"%.6f\" z=\"%.6f\" />\n",
                 @{$structure{'cell'}[2]}[1..3];
      printf OUT "    </Cell>\n";

      for( my $i = 0; $i < scalar( @{$structure{'atomic positions'}} ); 
           $i +=3 )
      {
        printf OUT "    <Atom x=\"%.6f\" y=\"%.6f\" z=\"%.6f\" type=\"%.6f\" />\n",
                 $structure{'atomic positions'}[$i]/2, 
                 $structure{'atomic positions'}[$i+1]/2, 
                 $structure{'atomic positions'}[$i+2]/2; 
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
#     printf OUT "  <Phenotype/> \n";

      my $minizertype;
      if ( $params{'minimizer'} =~ /linear/i )
      {
       $minimizertype = sprintf "<Minimizer type=\"linear\" maxeval=\"%i\" "; 
                        $params{'max calls'};
      }
      elsif ( $params{'minimizer'} =~ /GradientSA/ ) 
      {
        $minimizertype = sprintf "<TabooMinimizer type=\"GradientSA\" maxeval=%i ",
                                 $params{'max calls'}; 
      }
      elsif ( $params{'minimizer'} =~ /TabooSA/ ) 
      {
        $minimizertype = sprintf "<TabooMinimizer type=\"SA\" maxeval=%i ",
                                 $params{'max calls'}; 
      }
      elsif ( $params{'minimizer'} =~ /sa/i )
      {
       $minimizertype = sprintf "<Minimizer type=\"SA\" maxeval=\"%i\" ",
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
        printf OUT "            <Operators type=\"or\">\n";
        if ( $params{'GA style'} =~ /krossover/ ) 
          { printf OUT "              <kCrossover prob=\"0.75\" />\n"; }
        else
          { printf OUT "              <Crossover value=\"0.5\" prob=\"0.75\" />\n"; }
        printf OUT "              <Mutation value=\"0.05\" prob=\"0.25\" />\n";
        printf OUT "            </Operators>\n";
        printf OUT "        </TabooOp> \n";
        printf OUT "        %s />\n", $minimizertype;
        printf OUT "      </Operators>\n";
      }
      elsif($params{'GA style'} =~ /darwin/i ) 
      {
        printf OUT "    <TabooOp> \n";
        printf OUT "      <Operators type=\"or\">\n";
        if ( $params{'GA style'} =~ /krossover/ ) 
          { printf OUT "        <kCrossover prob=\"0.75\" />\n"; }
        else
          { printf OUT "        <Crossover value=\"0.5\" prob=\"0.75\" />\n"; }
        printf OUT "        <Mutation value=\"0.05\" prob=\"0.25\" />\n";
        printf OUT "      </Operators>\n";
        printf OUT "    </TabooOp> \n";
      }


      printf OUT "    <Selection size=2 />\n";

      printf OUT "    <Population size=\"%i\"/>\n",
                 $params{'GA'}{'population'};
      printf OUT "    <Offsprings rate=\"%.4f\"/>\n",
                 $params{'GA'}{'replace per generation'};
      if( $params{'CH'} =~ /one point/i )
        { printf OUT "    <OnePointHull/>\n"; }
      if( $params{'GA style'} =~ /restart/ )
      {
        printf OUT "    <PopAlgo rate=\"1.0\" every=\"%i\" >\n", 
                   $params{"every"};
          printf OUT "      <Operators type=\"and\">\n";
        printf OUT "          <UtterRandom/>\n";
          printf OUT "      </Operators>\n";
        printf OUT "    </PopAlgo>\n";
      }
      if( $params{'GA style'} =~ /best/ )
      {
        printf OUT "    <PopAlgo rate=\"%.2f\" every=\"%i\" >\n", 
                   $params{'rate'}, $params{"every"};
        printf OUT "      <Operators type=\"and\">\n";
        printf OUT "        %s \>\n", $minimizertype;
        printf OUT "        <UtterRandom/>\n";
        printf OUT "      </Operators>\n";
        printf OUT "    </PopAlgo>\n";
      }
      if( exists $params{'taboos'} )
      {
        printf OUT "    <Taboos>\n"; 
        if ( $params{'taboos'} =~ /pop/ )
          { printf OUT "      <PopTaboo/>\n";  }
        if ( $params{'taboos'} =~ /unique/ )
        {
          printf OUT "      <AgeTaboo lifespan=0 />\n";  
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
#     printf OUT "    <Statistics type=\"population\"/>\n";
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

