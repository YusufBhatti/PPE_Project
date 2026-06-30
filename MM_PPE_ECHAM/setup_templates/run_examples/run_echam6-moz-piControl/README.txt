This template provides a configuration to run ECHAM6-MOZ without HAM in a  
piControl setup from 1850 to 1880. The setup has been tested using a reduced 
set of chemical equations that differ from the default full chemical mechanism 
distributed with ECHAM6-HAMMOZ. To use the reduced set of chemical equations 
do the following BEFORE compiling ECHAM6-HAMMOZ

svn checkout https://svn.iac.ethz.ch/external/echam-hammoz/mozart-preproc
cd mozart-preproc/trunk/src/
./make_mozpp_dev # needs ifort
cd ../inputs
../bin/mozpp hammoz_gam001.in
cp -p ../output/gam/*.f90 <your echam root directory>/src/ 
# now recompile ECHAM6-HAMMOZ

template provided by Sebastian Wahl, GEOMAR Kiel (swahl@geomar.de)
