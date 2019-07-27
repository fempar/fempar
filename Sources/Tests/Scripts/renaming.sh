

# for i in $(find . -name *.[if]90); do grep femap $i; done

before=$1
after=$2

# Sources
#cd ~/git-repos/fempar/Sources
#for i in $(find . -name *.[if]90); do j=$(cat $i | grep $before );  if [ ! "$j" == "" ]; then cat $i | sed s/"$before"/"$after"/g > /tmp/kk; cp /tmp/kk $i; fi; done;

# Test scripts
#cd ~/git-repos/fempar/Sources/Tests/Scripts
for i in $(ls); do j=$(cat $i |  grep $before );  if [ ! "$j" == "" ]; then cat $i | sed s/"$before"/"$after"/g > /tmp/kk; cp /tmp/kk $i; fi; done;

# Drivers
#cd ~/git-repos/fempar/Drivers
#for i in $(find . -name *.[if]90); do j=$(cat $i | grep $before );  if [ ! "$j" == "" ]; then cat $i | sed s/"$before"/"$after"/g > /tmp/kk; cp /tmp/kk $i; fi; done;
