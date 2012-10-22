#!/bin/sh

cd src
make

cp ppt ../scripts
cp ppt ../S1

cd ../scripts
./do.sh

cd ../S1
./do.sh
