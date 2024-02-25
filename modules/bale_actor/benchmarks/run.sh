#!/bin/bash

start_node=2
max_node=4
ranks=2
size=s
numtimes=5
agi=0
upc=0
if [[ -z "$SRUN" ]]; then
   SRUN=srun
fi

while getopts s:n:r:t:m:a:u:h option
do
  case "${option}"
  in
    s) size=${OPTARG};;
    n) numtimes=${OPTARG};;
    r) ranks=${OPTARG};;
    t) start_node=${OPTARG};;
    m) max_node=${OPTARG};;
    a) agi=${OPTARG};;
    u) upc=${OPTARG};;
    h) echo "Options
             -s <s|m|l>  default is s i.e. use small input
             -n <1..inf> default is 5 i.e. run each experiment 5 times
             -r <1..inf> default is 2 i.e. two ranks per node
             -t <1..inf> default is 2 i.e. start number of nodes is 2
             -m <1..inf> default is 4 i.e. max number of nodes is 4
             -a <0|1> default is 0 i.e do not run AGI version
             -u <0|1> default is 0 i.e do not run UPC version
             -h show options"; exit 0;;
  esac
done

EXES=(histo ig permute randperm toposort transpose triangle)

if [ "$size" = s ]; then
    echo "Using small input"
    NS=(10000 10000 1000 10000 1000 1000 100)
elif [ "$size" = m ]; then
    echo "Using medium input"
    NS=(1000000 1000000 10000 100000 10000 10000 1000)
elif [ "$size" = l ]; then
    echo "Using large input"
    NS=(10000000 10000000 100000 1000000 100000 100000 10000)
else
    echo "Using small input (default)"
    NS=(10000 10000 1000 10000 1000 1000 100)
fi

for i in $(seq 0 6)
do

  EXE=${EXES[$i]}
  N=${NS[$i]}
  echo $EXE $N $i

  nodes=$start_node
  while [ $nodes -le $max_node ]; do

    count=$((nodes*ranks))

    echo "At $i , nodes $nodes"

    if [ $agi -eq 1 ]; then
    echo "${EXE}_agi"

    for c in $(seq 1 $numtimes)
    do
      echo "in_${EXE}_agi $c"
      $SRUN -n $count ./${EXE}_agi -n $N
      sleep 10
    done
    fi

    if [ $upc -eq 1 ]; then
    echo "${EXE}_upc"

    for c in $(seq 1 $numtimes)
    do
      echo "in_${EXE}_upc $c"
      $SRUN -n $count ./${EXE}_upc -n $N
      sleep 10
    done
    fi

    echo "${EXE}_shmem"

    for c in $(seq 1 $numtimes)
    do
      echo "in_${EXE}_shmem $c"
      $SRUN -n $count ./${EXE}_shmem -n $N
      sleep 10
    done

    echo "${EXE}_conveyor"
    for c in $(seq 1 $numtimes)
    do
      echo "in_${EXE}_conveyor $c"
      $SRUN -n $count ./${EXE}_conveyor -n $N
      sleep 10
    done

    echo "${EXE}_selector"

    for c in $(seq 1 $numtimes)
    do
      echo "in_${EXE}_selector $c"
      $SRUN -n $count ./${EXE}_selector -n $N
      sleep 10
    done

    nodes=$((nodes*2))
  done
done

