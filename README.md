# Quicktree Parallelism
We parallelize the Neighbor-Joining (NJ) algorithm from QuickTree on the GPU with CUDA.

## Results
| Input | Speedup vs. QuickTree (CPU) | Accuracy vs. reference tree |
|---|---|---|
| 10,000 taxa | **20×** | Normalized Robinson-Foulds distance **0.0013** (9,250 branches matched, 26 differ) |

Measured on an NVIDIA A30 (MIG 2g.12gb slice) on UCSD's DSMLP cluster. The full run log is in [`outputs/output.log`](outputs/output.log) and the tree comparison in [`outputs/tree_10k_RFdistances.txt`](outputs/tree_10k_RFdistances.txt).

## How it works
NJ repeatedly finds the closest pair of nodes in an N×N distance matrix, merges them, and updates the matrix. On the CPU each step is O(N²), and it runs N times.

- **Min-pair search on the GPU.** Each thread block reduces its slice of the matrix to a local minimum in shared memory (`findBlockMin`), then a second kernel reduces the block results to the global minimum (`findGlobalMin`).
- **Incremental updates.** After a merge, only the merged node's row changes, so we update that row and the per-node distance sums in O(N) (`updateDistanceMatrix`, `updateU`) instead of recomputing the matrix.
- **Active-node mask.** Merged nodes are marked inactive instead of compacting the matrix, which avoids moving memory on every step.

Team project for ECE 213 (Parallel Computing in Bioinformatics), UC San Diego, Winter 2026.

## Table of Contents
- [Quicktree Parallelism](#quicktree-parallelism)
  - [Table of Contents](#table-of-contents)
  - [Dataset](#dataset)
  - [ Requirements](#-requirements)
  - [Build Instructions](#build-instructions)
  - [Run the Program](#run-the-program)
  - [Execution Workflow on DSMLP Server](#execution-workflow-on-dsmlp-server)

## Dataset
Distance matrix with size 1K, 5K, and 10K are available within the data folder in the repository. Distance matrix in 20K size is available on [Google Drive](https://drive.google.com/file/d/1fgA6SpcN3XCDpbHKMmuSvPujOlT3G2lN/view?usp=sharing). Please download and move it within the data folder.

## <a name="req"></a> Requirements
The following dependencies are required:
- C++ compiler (GCC / Clang)
- CMake
- CUDA Toolkit


## Build Instructions
Create a build directory and compile the project:

```bash
mkdir -p build
cd build
cmake ..
make -j4
```

This will generate the executable:

```bash
quicktree
```

## Run the Program

Run the program using a distance matrix in PHILIP format as input (output will be printed to the console):

```bash
./quicktree -m ../data/matrix_10k.phy
```

Run the program using a distance matrix in PHILIP format as input (output will be written to the file):
```bash
./quicktree -m ../data/matrix_10k.phy -o tree.nwk
```

The program arguments are as follows:

| Argument | Description |
|----------|-------------|
| `-m` | input distance matrix (must in PHILIP format)|
| `-o` | output phylogenetic tree in Newick format |


## Execution Workflow on DSMLP Server 
1. Clone the project repository in HOME directory, and decompress the data files
   ```bash
    git clone https://github.com/wangly7/quicktree-parallelism.git
    cd quicktree-parallelism/data
    unzip matrix_10k.phy.zip
    cd ..
   ```

2. Run the code
     ```bash
    /opt/launch-sh/bin/launch.sh -v a30 -c 8 -g 1 -m 8 -i yatisht/ece213-wi26:latest -f ./quicktree-parallelism/run-commands.sh
    ```
