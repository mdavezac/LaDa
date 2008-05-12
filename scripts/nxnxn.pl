#!/usr/bin/perl
#
#  Version: $Id$
#

my %params;
read_params();

$params{"home"} = `cd; pwd`; chomp $params{"home"};

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

if ( $params{'GA string'} =~ /optimum/ )
  { $params{'method'}{'type'} = "optimum"; }
elsif ( $params{'GA string'} =~ /BestOf:(\S+)/ )
{
  $params{'method'}{'type'} =  "BestOf";
  $params{'method'}{'delta'} = $1;
}
elsif ( $params{'GA string'} =~ /target:(\S+)/ )
{
  $params{'method'}{'type'} = "Target";
  $params{'method'}{'target'} = $1;
  if( $params{'method'}{'string'} !~ /delta:(\S+)/ )
    { print " Forgot to include delta in method target \n"; exit; }
  $params{'method'}{'delta'} = $1;
}
elsif ( $params{'GA string'} =~ /CH/ )
  { $params{'method'}{'type'} = "CH"; }
else
  { print "Couldn't determine method on input!!\n"; exit; }


$params{'breeding'} .= " utterRandom" if ( $params{'GA string'} =~ /utterRandom/ );
$params{'breeding'} .= " crossover"   if ( $params{'GA string'} =~ /crossover/ );
$params{'breeding'} .= " krossover"   if ( $params{'GA string'} =~ /krossover/ );
$params{'breeding'} .= " SA" if ( $params{'GA string'} =~ /SA/ );
$params{'breeding'} .= " VA" if ( $params{'GA string'} =~ /VA/ );
$params{'breeding'} .= " Beratan" if ( $params{'GA string'} =~ /Beratan/ );
if ( $params{'GA string'} =~ /mut(ation|):(\S+)(-|_)(\S+)/ )
 { $params{'breeding'} .= " mut:$2-$4"; }
  
$params{'taboos'} .= " pop" if ( $params{'GA string'} =~ /poptaboo/ );
$params{'taboos'} .= " history" if ( $params{'GA string'} =~ /history/ );
if ( $params{'GA string'} =~ /c:(\S+)-(\S+)/ )
{
  $params{'taboo'} .= " concentration:$1-$2";
}

if ( $params{'GA string'} =~ /terminator:(\d+)/ )
  { $params{'terminator'} = $1; }

$params{'other'} = "";
$params{'other'} = " history " if ( $params{'GA string'} =~ /history/ );
$params{'other'} = " stats " if ( $params{'GA string'} =~ /stats/ );

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

my $taboos = $params{'taboos'};
$taboos =~ s/^(\s+)//;
$taboos =~ s/\s/_/g;
$params{'agr'}{'filename'} .= "_taboos:" . $taboos;

$params{'xml'}{"filename"} = $params{'agr'}{'filename'};

print $params{'agr'}{"filename"}, "\n"; 

# begin work
my %structure;

my $cosa = 0;
#system "rm -f $params{'agr'}{'filename'}";
read_structure(); # data passes from PIfile to %structure hash

write_lamarck_input();

exit;


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
      my $space = "";
      if ( $params{'taboos'} =~ /(pop|history|concentration)/ )
      {
        $space = " ";
        printf OUT "      <TabooOp>\n";
      }
      printf OUT "  %s      <UtterRandom/>\n", $space if ( $params{"breeding"} =~ /utterRandom/ );
      printf OUT "  %s      <Krossover />\n", $space if ( $params{"breeding"} =~ /krossover/ );
      printf OUT "  %s      <Crossover/>\n", $space if ( $params{"breeding"} =~ /crossover/ );
      if ( $params{'breeding'} =~ /mut(ation|):(\S+)-(\S+)/ )
        { printf OUT "        <Mutation prob=\"%5.4f\" value=\"%5.4f\" />\n", $2, $3; }
      printf OUT "      </TabooOp>\n"
        if ( $params{'taboos'} =~ /(pop|history|concentration)/ );
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

      printf OUT "    <Statistics/>\n" if(    $params{'other'} =~ /stats/ );

      printf OUT "    <Terminator ref=\"evaluation\" value=%i/>\n", $params{'terminator'}; 
      printf OUT "    <Save what=\"all\"/>\n";
      printf OUT "    <Filenames save=\"%s.xml\"/>\n", $params{'xml'}{'filename'};
      printf OUT "    <Filenames xmgrace=\"%s/%s.agr\" />\n", 
                      $params{"directory"}{'working'}, $params{'xml'}{'filename'};

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


sub read_params()
{
  open IN, "COMMANDS";
  $params{'agr'}{'filename'}      =  <IN>;
  $params{'directory'}{'result'}  =  <IN>;
  $params{'directory'}{'working'}  =  <IN>;
  while( ($_=<IN>) )
  {
     my $text = $_; chomp $text;

     $text =~ s/(^\s+|\s+$)//g;
     $params{'GA string'} .= " $text";
  }
  close IN;

  chomp $params{'agr'}{'filename'};
  chomp $params{'directory'}{'result'};
  chomp $params{'directory'}{'working'};
  chomp $params{'GA string'};
  $params{'agr'}{'filename'}     =~ s/(^\s+|\s+$)//g;
  $params{'directory'}{'result'} =~ s/(^\s+|\s+$)//g;
  $params{'directory'}{'working'} =~ s/(^\s+|\s+$)//g;
  $params{'GA string'}           =~ s/(^\s+|\s+$)//g;
}



