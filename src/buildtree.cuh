#include <vector>
#include <tuple>
#include "tree.hpp"

struct GpuTree {
    // CPU memory
    std::vector<double>& data;
    std::vector<TNode*> nodes;
    std::vector<string>& identifies;
    uint32_t N;
    Tree* theTree;

    // GPU memory
    double *mat;  // input
    double *u_i; // input

    uint32_t *d_block_min_ij; // min ij from different blocks
    double* d_block_min_dst; // min distance from different blocks

    uint32_t *d_min_ij; // min i and min j packed into one
    double *d_merged_dst; // [0] distance from i to ij
                         // [1] distance from j to ij

    uint32_t *active;  // active nodes that need to merge

    void allocateMem();
    void transferMat2Device();
    void build();
    void initNodesOnCPU();
    MergeInfo transferNode2Host();
    void buildInternalNode(const MergeInfo& info);
    void clearAndReset();
};

struct RdcDst {
    uint32_t index_i;
    uint32_t index_j;
    double dst;
};


struct MergeInfo{
    uint32_t min_ij[2];
    double branch_len[2];
};