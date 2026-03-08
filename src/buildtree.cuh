#include <vector>
#include <tuple>

struct GpuTree {
    // CPU memory
    std::vector<double>& data; 
    uint32_t h_size;
    uint32_t *h_nodes;
    uint32_t *h_min_ij;
    double *h_merged_dst;

    // GPU memory
    double *mat;  // input
    double *u_i; // input
    uint32_t d_size;  // input

    uint32_t *d_min_ij; // min i and min j packed into one
    double *d_min_dst; // min distance
    double *d_merged_dst; // [0] distance from i to ij
                         // [1] distance from j to ij

    uint32_t *active;  // active nodes that need to merge

    void allocateMem();
    void transferMat2Device();
    std::vector<double> buildTree (std::vector<double>);
    std::vector<double> transferTree2Host(std::vector<double>& data, uint32_t size);
    void clearAndReset();
};