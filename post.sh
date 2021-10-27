#!/bin/bash
# poat.sh
e4shared --post --job=heatshield --vtk-xml --tindx-plot=all \
	 --add-vars="mach,pitot,enthalpy, total-p, total-h," \ 
	 
# 
e4shared --post --job=heatshield --tindx-plot=last --add-vars="enthalpy" \
    --output-file= surface-blk0.dat \
    --slice-list="0,$,:,0"

e4shared --post --job=heatshield --tindx-plot=last \
    --output-file= surface-blk1.dat \
    --slice-list="1,$,:,0"

e4shared --post --job=heatshield --tindx-plot=last \
    --output-file= surface-blk2.dat \
    --slice-list="2,$,:,0"

gnuplot surface-pressure.gnuplot

awk -f scale-heat-flux.awk ./loads/t10-loads.dat > stanton.data
gnuplot surface-heat-transfer.gnuplot



	


