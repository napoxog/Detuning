#!/bin/bash

export filename=$1.res
./Detuning.py $1 $2 $3 > $filename

gnuplot -e '
file=system("echo $filename");
print file;
set title file;
plot 
file  using  ($1):($2) w l t "Composite Amplitude",
file  using ($1):($3) w l t "No-tune curve",
file  using ($1):($5) w l t "Transfer Function";
pause mouse button1 "Left Mouse button will terminate"'