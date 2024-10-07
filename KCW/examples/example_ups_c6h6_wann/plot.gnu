set xlab "Binding Energy [eV]"
set ylab "Intensity [a.u.]"
pl [18:8] 'results/photospectra_ki.dat'    u ($1):($2/0.7770) w l lw 2  lc rgb 'red'   tit 'KI\@MLWFs', \
          'results/photospectra_pkipz.dat' u ($1):($2/0.7123) w l lw 2  lc rgb 'blue'  tit 'pKIPZ\@MLWFs', \
          'reference/exp.dat'                       w p pt 7 ps 2       lc rgb 'black' tit 'Exp.'
pause -1 

pl [18:8] 'reference/photospectra_ki.dat'    u ($1):($2/0.7770) w l lw 2 lc rgb 'red'   tit 'KI\@MLWFs', \
          'reference/photospectra_pkipz.dat' u ($1):($2/0.7123) w l lw 2 lc rgb 'blue'  tit 'pKIPZ\@MLWFs', \
          'reference/exp.dat'                       w p pt 7 ps 2        lc rgb 'black' tit 'Exp.'
pause -1

