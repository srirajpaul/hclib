#!/bin/bash

start_node=2
max_node=4
ranks=2
numtimes=5
agi=0
upc=0

while getopts i:n:r:t:m:a:u:h option
do
  case "${option}"
  in
    i) FILE=${OPTARG};;
    n) numtimes=${OPTARG};;
    r) ranks=${OPTARG};;
    t) start_node=${OPTARG};;
    m) max_node=${OPTARG};;
    a) agi=${OPTARG};;
    u) upc=${OPTARG};;
    h) echo "Options
             -i <input file>
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

count=$(cat $FILE |grep seconds|wc -l)
#echo $count

versions=$((3+agi+upc))
#echo "ver " $versions
expt_node=$((max_node*2/start_node))
#echo "expt " $expt_node
num_node_count=$(echo $expt_node | awk '{print log($1)/log(2)}')
actual_count=$((numtimes*num_node_count*versions*7))
#echo $actual_count

if [ $count -ne $actual_count ]; then
echo
echo "WARNING: Requires " $actual_count " readings but only " $count " found"
echo
fi

avg_times=$(cat $FILE |grep seconds|awk -v n_times="$numtimes" '{sum+=$1} NR%n_times==0 {print sum/n_times; sum=0}')
#echo $avg_times

i=0
bm_count=$((versions*num_node_count))
bm_id=0
node_count=2

echo
echo "            HClib Actor/Selector"
echo
#echo "        SHMEM        Conveyor     Selector"

if [ "$agi" = 1 ]; then
    printf "         AGI"
fi
if [ "$upc" = 1 ]; then
    printf "         UPC"
fi
echo "         SHMEM         Conveyor      Selector"

for str in $avg_times
do
    #print benchmark name
    rem=$((i%bm_count))
    if [ "$rem" = 0 ]; then
    echo
        echo ${EXES[$bm_id]}
        bm_id=$((bm_id+1))
        node_count=2
    fi

    #print time
    printf "%12.4f  " $str

    #print node count after $versions reading
    i=$((i+1))
    rem=$((i%versions))
    if [ "$rem" = 0 ]; then
        num_ranks=$((node_count*ranks))
        printf "%4d PEs\n" $num_ranks
        node_count=$((node_count*2))
    fi
done

echo

