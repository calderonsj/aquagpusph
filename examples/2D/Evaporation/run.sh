#!/bin/bash

rm -f Fluid.dat
rm -f AQUAgpusph.save.*.xml
rm -f log.*.html
rm -f output.*.vtu
rm -f output.pvd
rm -f Timing.dat
rm -f Forces.dat
rm -f Performance.dat
./Create.py
/home/calderon/SPH/aquagpusph.cmake/bin/AQUAgpusph2D -i Main.xml
