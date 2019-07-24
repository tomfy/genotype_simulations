#set terminal png
#set out 'inv0.17-0.5.png'

set lmargin 8

set grid
set key bottom

red = 1
darkred = 2
yellow = 3
green = 4
darkgreen = 5
blue = 6
purple = 7
black = 8

set style data points
set pointsize 0.25

set title '0.17<m<0.5, 1/m'
set xlabel 'HGMR' offset 0,0.1
set ylabel 'AGMR' offset 9,0

shape = 5

# 
plot \
        'PO_inv0.17-0.5_400' using 7:6 t'Parent offspring' pt shape lc 3, \
	'FS_inv0.17-0.5_400' using 7:6 t'Full Sib' pt shape lc 1, \
	'AUNN_inv0.17-0.5_400' using 7:6 t'AU-NN' pt shape lc 2, \
	'GPGC_inv0.17-0.5_400' using 7:6 t'GPGC' pt shape lc 4, \
	'HS_inv0.17-0.5_400' using 7:6 t'Half sib' pt shape lc 5, \
	'FC_inv0.17-0.5_400' using 7:6 t'1st cousin' pt shape lc 6, \
        'HAUNN_inv0.17-0.5_400' using 7:6 t'Half AU-NN' pt shape lc 8, \
        'UN_inv0.17-0.5_400' using 7:6 t'Unrelated' pt shape lc 7
pause -1

plot \
        'PO_inv0.17-0.5' using 7:6 t'Parent offspring' pt shape lc 3, \
	'FS_inv0.17-0.5' using 7:6 t'Full Sib' pt shape lc 1, \
	'AUNN_inv0.17-0.5' using 7:6 t'AU-NN' pt shape lc 2, \
	'GPGC_inv0.17-0.5' using 7:6 t'GPGC' pt shape lc 4, \
	'HS_inv0.17-0.5' using 7:6 t'Half sib' pt shape lc 5, \
	'FC_inv0.17-0.5' using 7:6 t'1st cousin' pt shape lc 6, \
        'HAUNN_inv0.17-0.5' using 7:6 t'Half AU-NN' pt shape lc 8, \
        'UN_inv0.17-0.5' using 7:6 t'Unrelated' pt shape lc 7

# pause -1