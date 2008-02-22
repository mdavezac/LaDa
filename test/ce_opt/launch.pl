#! /usr/bin/perl

exit if( not -e "../../darwin/ce_opt" );
my $argument = shift;

system "rm calc/krossover_VA_gen:20.agr" if ( -e "calc/krossover_VA_gen:20.agr" );
system "cd calc; ln -s ../../../darwin/ce_opt " if ( not -e "calc/ce_opt" );
chdir "calc";
if( $argument =~ /kdbg/ )
{
  system "ln -s ../input.xml" if( not -e "input.xml" );
  system "kdbg ce_opt";
}
elsif( $argument =~ /make/ )
{
  system "cd ../../../; make";
}
elsif( $argument =~ /one/ )
{
    system "./ce_opt -i ../input.xml";
}
else
{ 
  for( my $i = 0; $i < 200; $i++)
  {
    system "./ce_opt -i ../input.xml";
    system "cat ../calc/out.xmg >> calc/krossover_VA_gen:20.agr";
  }
}
system "rm calc/ce_opt" if ( -e "calc/ce_opt" );

