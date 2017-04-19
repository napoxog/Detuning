python Detuning.py %1 "250" "0.5" > %1.res
set Detuning_res=%1.res
set Detuning_res=%Detuning_res:\=\\%
gnuplot -e "set title \"%Detuning_res%\"; plot  \"%Detuning_res%\"   using  ($1):($2) w l title 'Composite Amplitude', \"%Detuning_res%\"  using ($1):($3) w l title 'No-tune curve', \"%Detuning_res%\"  using ($1):($5) w l title 'Transfer Function'; pause mouse button1 \"Left Mouse button will terminate\""
