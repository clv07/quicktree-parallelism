#include "buildtree.cuh"
#include <cfloat>
#include <cstdint>
#include <limits>
#include <tuple>
#define BLOCKSIZE 256
#define NUMBLOCKS 8


__global__ void initActive(uint32_t d_size, uint32_t* active) {
    for (uint32_t i=bs*bx+tx; i < d_size; i+=bs*gs) {
        active[i] = 1;
    }
}

__global__ void blockReduction (uint32_t* d_min_ij, double* d_min_dst) {
    // parallel reduction across block
    for (int s=NUMBLOCKS/2; s>0; s>>=1) { 
        if (bx < s) {
            if(d_min_dst[bx].dst < d_min_dst[bx+s].dst) {
                d_min_dst[bx] = d_min_dst[bx+s];
                d_min_ij[bx] = d_min_ij[bx+s];        
            }
        } 
    }
}

/**
 * Building phylogenetic tree on GPU
 */
 __global__ void buildTreeOnCPU(uint32_t d_size, double *mat, double* u_i, uint32_t* active, uint32_t* d_min_ij, double* d_min_dst,  double *d_merged_dst){
    // Kernel Configuration
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int bs = blockDim.x;
    int gs = gridDim.x;
    
    struct rdc_dst {
        uint32_t index_i;
        uint32_t index_j;
        double dst;
    };
    __shared__ rdc_dst reduction[BLOCKSIZE]; 

    // one grid for one distance matrix
    // 1. for each leaf compute u_i = summation from j!= i to N, (d_ij) / (N-2)
    for (uint32_t i = bs*bx+tx; i<d_size; i+=bs*gs) {
        if (active[i] == 0) continue;
        double sum_ui = 0.0;
        for (uint32_t j=0; j<i; j++){
            if (active[j] == 0) continue;
            uint32_t index = i*(i-1)/2+j;
            double d_ij = mat[index]; 
            sum_ui += d_ij;
        }
        u_i[i] = sum_ui / (d_size-2);
    }

    // initialize reduction within blocks 
    reduction[tx].index_i = -1;
    reduction[tx].index_j = -1;
    reduction[tx].dst = DBL_MAX;

    // 2. find leaf pair i and j for which d_ij - u_i - u_j is minimum 
    for (uint32_t i = bs*bx+tx; i<d_size; i+=bs*gs) {
        if (active[i] == 0) continue;
        for (uint32_t j=0; j<i; j++) {
            if (active[j] == 0) continue;
            uint32_t index = i*(i-1)/2+j;
            // if (d_nodes[index] == 0) continue; // skip inactive nodes
            double leaf = mat[index] - u_i[i] - u_i[j];
            if (leaf < reduction[tx].dst) { 
                reduction[tx].dst = leaf;
                reduction[tx].index_i = i;
                reduction[tx].index_j = j;
            }
        }
    }

    // parallel reduction within a block
    for (int s=BLOCKSIZE/2; s>0; s>>=1) { 
        if (tx < s) {
            if(reduction[tx].dst > reduction[tx+s].dst) reduction[tx] = reduction[tx+s];
        }
        __syncthreads();
    }

    if (tx == 0) {
        d_min_ij[2*bx] = reduction[tx].index_i;
        d_min_ij[2*bx+1] = reduction[tx].index_j;
        d_min_dst[bx] = reduction[tx].dst;
    }
    __syncthreads();

    // Problem 1
    // parallel reduction across block
    for (int s=NUMBLOCKS/2; s>0; s>>=1) { 
        if (bx < s) {
            if(d_min_dst[bx].dst < d_min_dst[bx+s].dst) {
                d_min_dst[bx] = d_min_dst[bx+s];
                d_min_ij[bx] = d_min_ij[bx+s];        
            }
        __syncthreads();
    } 

    // 3. Join i and j to form a common node (ij) such that
    // the distance of i from the (ij), d_i(ij) = 0.5 * (d_ij + (u_i-u_j)) and 
    // the distance of j from the (ij), d_j(ij) = d_ij - d_i(ij)
    double min_dst = d_min_dst[0];
    double min_i = d_min_ij[0][0]; // Problem 2
    double min_j = d_min_ij[0][1];
    double d_i_ij = 0.5 * (min_dst + (u_i[min_i] - u_i[min_j]));
    double d_j_ij = mid_dst - d_i_ij;

    active[min_i] = 1; // which indicate new nodes is active
    active[min_j] = 0; // merged 

    // adjustment for negative branch length
    if (d_i_ij < 0.0) {
        d_i_ij = 0.0;
        d_j_ij = (min_dst < 0.0)? 0.0 : d_ij;
    }
    else if (d_j_ij < 0.0) {
        d_j_ij = 0.0;
        d_i_ij = (min_dst < 0.0)? 0.0 : d_ij;
    }

    // copy distances of i and j from ij to variable in global memory
    d_merged_dst[0] = d_i_ij;
    d_merged_dst[1] = d_j_ij;

    // Update distance matrix
    // 4. Treat (ij) as a new tip, ignoring previous tips i and j
    // 5. Distance of (ij) from other tips k, d_k(ij) = (d_ik+d_jk-d_ij) * 0.5 
    uint32_t a = std::min(min_i, min_j);
    uint32_t b = std::max(min_i, min_j);

    // set new merged distance for (a, k) for all tips k 
    for (uint32_t k=bs*bx+tx; k<d_size; k+=gs*bs) {
        if (active[k] == 0) continue;
        uint32_t ak_index = a*(a-1)/2+k;
        uint32_t bk_index = b*(b-1)/2+k;
        if (k != a && k != b){
             double new_d = (mat[ak_index] + mat[bk_index] - min_dst) * 0.5;
             mat[ak_index] = (new_d < 0.0)? 0.0: new_d; // adjust for negative branch length 
        }
    }

    d_size = d_size - 1; // Problem 3
    
 }



/**
 * Main orchestration method.
 * 1. Allocates GPU memory
 * 2. Transfers data from CPU
 * 3. Launches Kernel
 * 4. Retrieves results and reconstructs trees
 */
 void GpuTree::build() {
    int numBlocks = NUMBLOCKS;  // i.e. number of thread blocks on the GPU
    int blockSize = BLOCKSIZE; // i.e. number of GPU threads per thread block

    // 1. Allocate memory on Device
    allocateMem();

    // 2. Transfer sequence to device
    transferMat2Device();

    // 3. Perform building tree on GPU
    buildTreeOnCPU<<<numBlocks, blockSize>>>(d_size, mat, d_nodes, u_i, d_min, d_minj);
    

    // 4. Transfer the final merged nodes from device
    transferTree2Host();

    // 5. Get the phylogenetic tree with merged nodes
    
 }

void GpuTree::getPhylogenneticTree () {
    
}

/**
 * Allocates GPU memory for distance matrix 1D array, nodes status, and merged nodes.
 */
void GpuTree::allocateMem() {
    cudaError_t err;

    // 1. allocate memory to distance matrix
    err = cudaMalloc(&mat, (h_size * (h_size - 1) / 2) * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 2. allocate memory to d_size
    err = cudaMalloc(&d_size, sizeof(uint32_t));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 3. allocate memory to nodes status array 
    err = cudaMalloc(&active, h_size * sizeof(uint32_t));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }
    
    // 5. allocate memory for u_i
    err = cudaMalloc(&u_i, h_size * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 6. allocate struct for min_i, min_j
    err = cudaMalloc(&d_min_ij, 2 * NUMBLOCKS * sizeof(uint32_t));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 7. allocate struct for min_dst
    err = cudaMalloc(&d_min_dst, NUMBLOCKS * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 8. allocate struct for distance of min i and min j from merged ij
    err = cudaMalloc(&d_merged_dst, 2 * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }
}

/*
Transfer distance matrix data and size to device
*/
void GpuTree::transferMat2Device() {
    cudaError_t err;

    // Problem 4
    err = cudaMemcpy(mat, data, (h_size*(h_size-1)/2) * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // Problem 5
    err = cudaMemcpy(d_size, h_size, sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }
}

/**
 * Copies the computed merged nodes information from GPU back to Host.
 */
std::vector<uint32_t> GpuTree::transferTree2Host() {
    std::vector<uint32_t> nodes;

    cudaError_t err = cudaMemcpy(h_nodes, d_nodes, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // transfer merged nodes distance from device to host
    cudaError_t err = cudaMemcpy(h_merged_dst, d_merged_dst, 2 * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }
    return nodes;
}

void GpuTree::clearAndReset() {
    cudaFree(mat);
    cudaFree(u_i);
    cudaFree(d_nodes);
    cudaFree(d_leftovers);
}