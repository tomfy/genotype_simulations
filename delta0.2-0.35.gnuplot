set terminal png
set out 'delta0.2-0.35.png


set key bottom
set grid
# set style data dots
set pointsize 0.5

set lmargin 8
set xlabel 'hgmr' offset 0,0.5
set ylabel 'agmr' offset 9,0


plot \
     'FS_delta0.2' using 7:6 t'FS0.2' pt 12 lc 1, 'FS_delta0.35' using 7:6 t'FS0.35' pt 4 lc 1, \
 'HS_delta0.2' using 7:6 t'HS0.2' pt 12 lc 2, 'HS_delta0.35' using 7:6 t'HS0.35' pt 4 lc 2, \
 'FC_delta0.2' using 7:6 t'FC0.2' pt 12 lc 5, 'FC_delta0.35' using 7:6 t'FC0.35' pt 4 lc 5, \
'UN_delta0.2' using 7:6 t'UN0.2' pt 12 lc 7, 'UN_delta0.35' using 7:6 t'UN0.35' pt 4 lc 7
# pause -1

#  'GPGC_delta0.2' using 7:6 t'GPGC0.2' pt 7 lc 4, 'GPGC_delta0.35' using 7:6 t'GPGC0.35' pt 4 lc 4, \