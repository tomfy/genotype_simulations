set key bottom
set grid
# set style data dots
set pointsize 0.4


plot \
     'FS_delta0.2' using 7:6, 'FS_delta0.35' using 7:6, \
 'HS_delta0.2' using 7:6, 'HS_delta0.35' using 7:6, \
 'GPGC_delta0.2' using 7:6, 'GPGC_delta0.35' using 7:6, \
 'FC_delta0.2' using 7:6, 'FC_delta0.35' using 7:6
pause -1

splot \
     'FS_delta0.2' using 2:5:4, 'FS_delta0.35' using 2:5:4, \
 'HS_delta0.2' using 2:5:4, 'HS_delta0.35' using 2:5:4, \
 'GPGC_delta0.2' using 2:5:4, 'GPGC_delta0.35' using 2:5:4, \
 'FC_delta0.2' using 2:5:4, 'FC_delta0.35' using 2:5:4
pause -1

set xlabel 'A mismatches'
set zlabel 'H mismatches'
set ylabel 'H denom'
splot [*:*][*:*] \
	'FS_inv0.17-0.5' using 2:5:4, \
	'AUNN_inv0.17-0.5' using 2:5:4, \
	'GPGC_inv0.17-0.5' using 2:5:4, \
	'FC_inv0.17-0.5' using 2:5:4
pause -1

set xlabel 'H denom'
set ylabel 'H mismatches'

plot \
	'FS_inv0.17-0.5' using 5:4, \
	'AUNN_inv0.17-0.5' using 5:4, \
	'GPGC_inv0.17-0.5' using 5:4, \
	'HS_inv0.17-0.5' using 5:4, \
	'FC_inv0.17-0.5' using 5:4, \
	'HAUNN_inv0.17-0.5' using 5:4
pause -1

set xlabel 'A mismatches'
set zlabel 'H mismatches'
set ylabel 'H denom'
splot [*:*][*:*] \
	'FS_inv0.17-0.5' using 2:5:4 t'FS', \
	'AUNN_inv0.17-0.5' using 2:5:4 t'AUNN', \
	'GPGC_inv0.17-0.5' using 2:5:4 t'GPGC', \
	'HS_inv0.17-0.5' using 2:5:4 t'HS', \
	'FC_inv0.17-0.5' using 2:5:4 t'FC', \
	'HAUNN_inv0.17-0.5' using 2:5:4 t'HAUNN'
pause -1


set xlabel 'HGMR'
set ylabel 'AGMR'
plot \
	'FS_inv0.17-0.5' using 7:6 t'Full Sib', \
	'AUNN_inv0.17-0.5' using 7:6 t'AU-NN', \
	'GPGC_inv0.17-0.5' using 7:6 t'GPGC', \
	'HS_inv0.17-0.5' using 7:6 t'Half sib', \
	'FC_inv0.17-0.5' using 7:6 t'1st cousin', \
        'HAUNN_inv0.17-0.5' using 7:6 t'Half AU-NN'

pause -1
