#
cd LAUNCH_DIR
#
############################
# Generate machinefile
srun hostname > hostfile_q
#mpirun -n NNO hostname > hostfile_q
hosts=$(sort -u hostfile_q)
nodes=NNO
odd=ODD
if [ $odd ]
then
    let nodes=nodes-1
fi
echo oddness is $odd and there are $nodes nodes
i=1
while [ $i -le $nodes ]
do
    host=`echo $hosts | sed s/" "" "*/#/g | cut -f $i -d#`
    first=0
    last=$(expr CPT - 1 )
    map="$first-$last"
    #map=$first
    if [ CPT -lt 25  ]
    then
        for j in `seq CPT CPT 23`
        do
	    let k=$j+CPT-1
	    map="$map,$j-$k"
	    #map="$map,$j"
        done
        for j in `seq 24 CPT 47`
        do
	    let k=$j+CPT-1
	    map="$map,$j-$k"
	    #map="$map,$j"
        done
    fi
    echo "$host:TPN binding=map=$map" >> machinefile
    #echo "$host:TPN binding=map=$map;domain=CPT" >> machinefile
    #echo "$host:TPN binding=domain=CPT" >> machinefile
    let i=i+1
done
# Coarse task on an extra node when ODD is defined
LT=-1
if [ $odd ] 
then
    host=`echo $hosts | sed s/" "" "*/#/g | cut -f $i -d#`
    echo "$host:1 binding=map=0-47;domain=48" >> machinefile
    let LT=LT+TNT
fi
###########################

# Run NUM_REP repetitions of the experiment
ID_REP=1
while  [ $ID_REP -le NRP ]
do
    echo $command_to_run
    eval $command_to_run > "stdout$ID_REP" 2> "stderr$ID_REP"
    let ID_REP=ID_REP+1
done
