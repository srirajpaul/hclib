## A top-down BFS implementation of actor

### Step 1: Clone the Graph500 reference implementation repository

```
git clone git@github.com:graph500/graph500.git
```

### Step 2: Replace `Makefile`, `main.c`, and `bfs_custom.c` with the ones in this repo

```
cd graph500/src
cp ../../Makefile .
cp ../../bfs_custom.c .
cp ../../main.c .
```

### Step 3: Build

```
make graph500_custom_bfs
```

### Step 4: Run

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD
$OSHRUN -n 2 ./graph500_custom_bfs 10
```