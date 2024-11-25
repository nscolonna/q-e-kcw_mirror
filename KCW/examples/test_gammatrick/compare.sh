#!/bin/bash
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

tol_ene=0.000001
tol_eig=0.005
tol_a=0.005
check=1
echo 
echo " Testing the results ..." 
echo "  tolerances: "
printf "     energy %10.7f Ry  \n" $tol_ene
printf "     eigval %10.7f eV  \n" $tol_eig
printf "     alpha  %10.7f     \n" $tol_a

## CHECK PWSCF
a=`grep ! results_gamma/h2o.scf.out | awk '{print $5}'`
a_ref=`grep ! results/h2o.scf.out | awk '{print $5}'`
err=`echo $a $a_ref | awk '{printf "%20.15f \n", sqrt(($1-$2)*($1-$2))}'`
if (( $(echo "$err > $tol_ene" |bc -l) )); then
 echo -e "${RED}  PWSCF  WARNING: ${NC} Total energy does not match: $a $a_ref $err"
 check=0
fi
if (( $check )); then echo -e "  ${GREEN}PWSCF       OK!${NC}"; fi


## CHECK PWSCF Eig
check=1
a=`grep high results_gamma/h2o.scf.out | awk '{print $7}'`
a_ref=`grep high results/h2o.scf.out | awk '{print $7}'`
err=`echo $a $a_ref | awk '{printf "%20.15f \n", sqrt(($1-$2)*($1-$2))}'`
if (( $(echo "$err > $tol_eig" |bc -l) )); then
 echo -e "${RED}  PWSCF  WARNING: ${NC} e_HO does not match: $a $a_ref $err"
 check=0
fi
a=`grep high results_gamma/h2o.scf.out | awk '{print $8}'`
a_ref=`grep high results/h2o.scf.out | awk '{print $8}'`
err=`echo $a $a_ref | awk '{printf "%20.15f \n", sqrt(($1-$2)*($1-$2))}'`
if (( $(echo "$err > $tol_eig" |bc -l) )); then
 echo -e "${RED}PWNSCF WARNING: ${NC} e_LU does not match: $a $a_ref $err"
 check=0
fi
if (( $check )); then echo -e "  ${GREEN}PWNSCF      OK!${NC}"; fi


for which_odd in "qki" "ki" "pkipz"; do
check=1

grep "   KI    " results_gamma/h2o.kcw-ham_${which_odd}.out > pp
grep "   KI    " results/h2o.kcw-ham_${which_odd}.out >> pp
#cat pp
for j in `seq 1 8`; do 
 col=`echo $j+1 | bc`
 eig=`head -1 pp | awk -v col=$col '{print $col}'`
 eig_ref=`tail -1 pp | awk -v col=$col '{print $col}'`
 err=`echo $j $eig $eig_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
 if (( $(echo "$err > $tol_eig" |bc -l) )); then
  echo -e "  ${RED}KC_HAM $which_odd   WARNING: ${NC}eig does not match for orbital $j: $eig $eig_ref $err"
  check=0
 fi
done
if (( $check )); then echo -e "  ${GREEN}KC_HAM  $which_odd    OK!${NC}"; fi
rm pp

done

echo 
