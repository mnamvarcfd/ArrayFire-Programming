#!/bin/sh

# The first the script is execute, it may be need to execute first
# "sed -i 's/\r//' make.sh" 
# to convert endofline to unix format
# May also need to execute "chmod 755 make.sh" 
# to give right to write to the script make.sh

mkdir build_linux
cd build_linux
cmake ..
make
mv MAIN_MONTECARLO ../workingDirectory

