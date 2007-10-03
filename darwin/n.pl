#!/usr/bin/perl

my %params;
open IN, "COMMANDS";
$params{'iaga call'}            =  <IN>;
$params{'agr'}{'filename'}      =  <IN>;
$params{'directory'}{'result'}  =  <IN>;
$params{'GA string'}            =  <IN>;
$params{'method string'}        =  <IN>;
$params{'breeding string'}      =  <IN>;
$params{'taboos string'}        =  <IN>;
$params{'terminator string'}    =  <IN>;
$params{'other string'}         =  <IN>;
close IN;

$params{"home"} = `cd; pwd`; chomp $params{"home"};
$params{'max GA iters'} = 400;

$params{'GA'}{'maxgen'} = 0;
$params{'GA'}{'popsize'} = 100;
$params{'GA'}{'rate'} = 0.1;

$params{'method'}{'type'} = "convexhull";
$params{'method'}{'target'} = "0";
$params{'method'}{'delta'} = "0";

$params{'breeding'} = "";

$params{'taboos'} = "";

$params{'other'} = "";

$params{'terminator'} = 1;




if ( $params{'GA string'} =~ /maxgen:(\d+)/ )
  { $params{'GA'}{'maxgen'} = $1; }
if ( $params{'GA string'} =~ /popsize:(\d+)/ )
  { $params{'GA'}{'popsize'} = $1; }
if ( $params{'GA string'} =~ /rate:(\S+)/ )
  { $params{'GA'}{'rate'} = $1; }
if ( $params{'GA string'} =~ /x=(\S+)/ )
  { $params{'GA'}{'x'} = $1; }

if ( $params{'method string'} =~ /optimum/ )
  { $params{'method'}{'type'} = "optimum"; }
elsif ( $params{'method string'} =~ /BestOf:(\S+)/ )
{
  $params{'method'}{'type'} =  "BestOf";
  $params{'method'}{'delta'} = $1;
}
elsif ( $params{'method string'} =~ /target:(\S+)/ )
{
  $params{'method'}{'type'} = "Target";
  $params{'method'}{'target'} = $1;
  if( $params{'method'}{'string'} !~ /delta:(\S+)/ )
    { print " Forgot to include delta in method target \n"; exit; }
  $params{'method'}{'delta'} = $1;
}
elsif ( $params{'method string'} =~ /CH/ )
  { $params{'method'}{'type'} = "CH"; }
else
  { print "Couldn't determine method on input!!\n"; exit; }


$params{'breeding'} .= " UtterRandom" if ( $params{'breeding string'} =~ /UtterRandom/ );
$params{'breeding'} .= " Crossover"   if ( $params{'breeding string'} =~ /Crossover/ );
$params{'breeding'} .= " Krossover"   if ( $params{'breeding string'} =~ /Krossover/ );
$params{'breeding'} .= " SA" if ( $params{'breeding string'} =~ /SA/ );
$params{'breeding'} .= " VA" if ( $params{'breeding string'} =~ /VA/ );
$params{'breeding'} .= " Beratan" if ( $params{'breeding string'} =~ /Beratan/ );
  
$params{'taboos'} .= " pop" if ( $params{'taboos string'} =~ /pop/ );
$params{'taboos'} .= " history" if ( $params{'taboos string'} =~ /history/ );
if ( $params{'taboos string'} =~ /c:(\S+)-(\S+)/ )
{
  $params{'taboo'} .= " concentration:$1-$2";
}

if ( $params{'terminator string'} =~ /(\d+)/ )
  { $params{'terminator'} = $1; }

$params{'other'} = " history " if ( $params{'other string'} =~ /history/ );
$params{'restart'} = " results " if ( $params{'other string'} =~ /true/ );

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
}


$params{'agr'}{'filename'} = "";
$params{'agr'}{'filename'} = "true_" if ($params{'restart'} =~ /results/ );
$params{'agr'}{'filename'} .= sprintf "atoms:%i_%i_%i_", $params{'x'}, $params{'y'}, $params{'z'};
$params{'agr'}{'filename'} .= $params{'method'}{'type'};
$params{'agr'}{'filename'} .= sprintf ":%.2f_%.2f", $params{'target'}, $params{'delta'}
  if( $params{'method'}{'type'} =~ /target/ );
$params{'agr'}{'filename'} .= sprintf ":%.2f", $params{'delta'} if( $params{'method'}{'delta'} =~ /BestOf/ );

$params{'breeding'} =~ s/ /_/g;
$params{'agr'}{'filename'} .= $params{'breeding'};
$params{'breeding'} =~ s/_/ /;

$params{'agr'}{'filename'} .= sprintf "_pop:%i", $params{'GA'}{'popsize'}
    if ( $params{'GA'}{'popsize'} != 100 );
$params{'agr'}{'filename'} .= sprintf "_rep:%.2f", $params{'GA'}{'rate'}
    if ( $params{'GA'}{'rate'} != 0.1 );
$params{'agr'}{'filename'} .= "_partition"
    if ( $params{'GA string'} =~ /partition/ );

$params{'agr'}{'filename'} .= sprintf "_x=%.2f", $params{'GA'}{'x'} 
  if( $params{'GA string'} =~ /x=(\S+)/ );
if( $params{'taboo'} =~ /concentration:(\S+)-(\S+)/ )
{
  $params{'agr'}{'filename'} .= sprintf "_c:%.4f-%.4f", $params{'taboo'}{'c'} 
}

$params{'xml'}{"filename"} = $params{'agr'}{'filename'};
$params{'agr'}{"filename"} .= ".agr";
$params{'xml'}{"filename"} .= ".xml";

print $params{'agr'}{"filename"}, "\n"; 

# begin work
my %structure;

my $cosa = 0;
#system "rm -f $params{'agr'}{'filename'}";
read_structure(); # data passes from PIfile to %structure hash
launch_iaga(); 


exit;


sub launch_iaga()
{
  for( my $count =0; $count < $params{'max GA iters'}; $count++)
  {
    write_lamarck_input();
    system "$params{'iaga call'}";
    system "cat convex_hull.agr >> $params{'agr'}{'filename'}";
    system "cp convex_hull.xml $params{'xml'}{'filename'}  ";
    system "cp $params{'xml'}{'filename'} $params{'agr'}{'filename'}  $params{'directory'}{'result'} ";
  }
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
      printf OUT "  <GA ";
      printf OUT "maxgen=\"%i\" " if( $params{'GA'}{'maxgen'} != 0 );
      printf OUT "populate=\"partition\" " if( $params{'GA string'} =~ /partition/ );
      printf OUT "x=\"%.2f\" ", $params{"GA"}{"x"} if( $params{'GA string'} =~ /x=(\S+)/ );
      printf OUT "popsize=\"%i\" rate=\"%.2f\">\n",
                 $params{'GA'}{'popsize'}, 
                 $params{'GA'}{'rate'};
      printf OUT "    <Method type=\"convexhull\" "
        if ( $params{'method'}{"type"} =~ /CH/ );
      printf OUT "    <Method type=\"%s\" ", $params{'method'}{'type'}
        if ( $params{'method'}{"type"} !~ /CH/ );
      printf OUT "target=\"%.2f\" ", $params{'method'}{'target'}
        if( $params{'method'}{'type'} =~ /target/ );
      printf OUT "delta=\"%.2f\" ", $params{'method'}{'delta'}
        if( $params{'method'}{'type'} =~ /(target|Best Of|BO)/ );
      printf OUT " />\n";
      
      printf OUT "    <Breeding>\n";
      if ( $params{'taboos'} =~ /(pop|history|concentration)/ )
      {
        printf OUT "      <TabooOp>\n";
        printf OUT "        <UtterRandom/>\n" if ( $params{"breeding"} =~ /UtterRandom/ );
        printf OUT "        <Krossover />\n" if ( $params{"breeding"} =~ /Krossover/ );
        printf OUT "        <Crossover/>\n" if ( $params{"breeding"} =~ /Crossover/ );
        printf OUT "      </TabooOp>\n";
      }
      else
      { 
        printf OUT "      <UtterRandom/>\n" if ( $params{"breeding"} =~ /UtterRandom/ );
        printf OUT "      <Krossover/>\n" if ( $params{"breeding"} =~ /Krossover/ );
        printf OUT "      <Crossover/>\n" if ( $params{"breeding"} =~ /Crossover/ );
      }
      printf OUT "      <Minimizer type=\"VA\" />\n" if ( $params{"breeding"} =~ /VA/ );
      printf OUT "      <Minimizer type=\"SA\" />\n" if ( $params{"breeding"} =~ /SA/ );
      printf OUT "      <Minimizer type=\"Beratan\" />\n" if ( $params{"breeding"} =~ /Beratan/ );
      printf OUT "    </Breeding>\n";


      if( $params{'taboos'} =~ /(pop|history|concentration)/ )
      {
        printf OUT "    <Taboos>\n"; 
        printf OUT "      <PopTaboo/>\n"  if ( $params{'taboos'} =~ /pop/ );  
        printf OUT "      <OffspringTaboo/>\n" if ( $params{'taboos'} =~ /pop/ );   
        printf OUT "      <HistoryTaboo/>\n" if ( $params{'taboos'} =~ /history/ );   
        if ( $params{'taboo'} =~ /concentration:(\S+)-(\S+)/ )
        {
          printf OUT "      <Concentration morethan=\"%.4f\" lessthan=\"%.4f\" />\n", $1, $2;
        }
        printf OUT "    </Taboos>\n"; 
      }

      printf OUT "    <History/>\n" if(    $params{'other'} =~ /history/ 
                                        or $params{'taboos'} =~ /history/ );


      printf OUT "    <Terminator ref=\"evaluation\" value=%i/>\n", $params{'terminator'}; 
      printf OUT "    <Save what=\"all\"/>\n";
      printf OUT "    <Restart what=\"results\"/>\n" if ( $params{'restart'} =~ /results/ );
      printf OUT "    <Filenames save=\"convex_hull.xml\" />\n";

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

