before=$1
after=$2

cd ~/git-repos/fempar/Drivers/par_nsi/modules

for i in $(ls); do j=$(cat $i |  grep $before );  if [ ! "$j" == "" ]; then cat $i | sed s/"$before"/"$after"/g > /tmp/kk; cp /tmp/kk $i; fi; done;
