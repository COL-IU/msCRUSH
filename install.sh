#!/bin/bash

echo ""
echo "Installing SPLASH."

CDIR=`pwd`
TMPDIR="src/app"
cd $TMPDIR

echo 'Step 1:'
bash compile_splash_on_general_charge.sh
echo ''
echo 'Step 2:'
bash compile_generate_consensus_spectrum_for_splash.sh

cd $CDIR

if test -d bin; then
  read -n1 -p 'A local bin dir already exists, overwrite? (Y/N)' booleanYorN
  case $booleanYorN in
   y|Y) echo "" ; rm -f -r bin ; mkdir bin ;;
   n|N) echo "" ; echo "Only replacing binary files in local bin dir" ;;
   *) echo "" ; echo "Invalid Input "; exit 1 ;;
  esac
else 
  mkdir bin
fi




mv $TMPDIR/splash_on_general_charge bin/
mv $TMPDIR/generate_consensus_spectrum_for_splash bin/

echo ""
echo "Executables are now installed under bin/"
echo ""

exit 0

