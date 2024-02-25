#!/bin/bash

for app in histo ig permute randperm transpose triangle toposort;
do
    ./trans.sh ${app}
    make ${app}_selector_lambda.trans
done

./trans.sh ig 3
make ig_selector_lambda3.trans

./trans.sh permute 2
make permute_selector_lambda2.trans

./trans.sh randperm 2
make randperm_selector_lambda2.trans
