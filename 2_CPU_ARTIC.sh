#!/bin/bash

set -e

echo 'make sure you...'
echo '1 log in twnia3.nchc.org.tw'
echo '2 basecalling and barcoding completed'
echo '3 change the parameters in mini_cov.sh'
echo ''
echo 'press Ctl+c if NOT'

sleep 10

echo 'start ARTIC pipeline'
sbatch sbatch_minicov.sh 
