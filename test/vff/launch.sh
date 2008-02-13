#! /bin/bash
#

cd calc
if [ ! -e vff ]; then 
  ln -s ../../../vff/vff .
fi
if [ ! -e layered_vff ]; then 
  ln -s ../../../vff/layered_vff .
fi

./vff -i ../quaternary.xml >quaternary
cp atomic.config quatconfig
grep "Energy" quaternary
./layered_vff -i ../quaternary.xml #>layered_quaternary
# cp atomic.config layered_quatconfig
# grep "Energy" layered_quaternary
# ./vff -i ../ternary.xml> ternary
# cp atomic.config ternconfig
# grep "Energy" ternary
# ./layered_vff -i ../ternary.xml >layered_ternary
# grep "Energy" layered_ternary
# cp atomic.config layered_ternconfig


if [ -e vff ]; then 
  rm vff
fi
if [ -e layered_vff ]; then 
  rm layered_vff
fi
