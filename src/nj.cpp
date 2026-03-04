// NJ Algorithm
#include "distancemat.hpp"
#include <climits> 

int compute_tree(DistanceMatrix mat) {

    // 6. Repeat step 1 to 5 till a single node (root) remains
    int cycle = mat.size;
    for (int h=0; h<cycle; h++) { 
        int N = mat.size
        int u_i[N];
        int u_j[N];
    
        int min_leaf = INT_MAX;
        int min_i;
        int min_j;
    
        // 1. for each leaf compute u_i = summation from j!= i to N, (d_ij) / (N-2)
        for (int i=0; i<N; i++) {
            int sum_ui = 0;
            int sum_uj = 0;
    
            // compute u_i
            for (int j=i; j<N; j++){
                int d_ij = mat.get(i, j);
                int temp = d_ij / (N-2);
                sum_ui += temp;
            }
            // compute u_j
            for (int j=0; j<i; j++){
                int d_ij = mat.get(i, j);
                int temp = d_ij / (N-2);
                sum_uj += temp;
            }
    
            u_i[i] = sum_ui;
            u_j[i] = sum_uj;
        }
        
        // 2. find leaf pair i and j for which d_ij - u_i - u_j is minimum
        for (int i=0; i<N; i++) {
            int ui = u_i[i];
            for (int j=0; j<i; j++) {
                int uj = u_j[j];
                int d_ij = mat.get(i, j); 
                int leaf = d_ij - ui - uj;
                
                if (leaf < min_leaf) {
                    min_leaf = leaf;
                    min_i = i;
                    min_j = j;
                }
    
            }
        }
    
        // 3. Join i and j to form a common node (ij) such that
        // the distance of i from the (ij), d_i(ij) = 1/2 (d_ij + (u_i-u_j)) and 
        // the distance of j from the (ij), d_j(ij) = d_ij - d_i(ij)
        int d_ij = mat.get(min_i, min_j); 
        int ui = u_i[min_i];
        int uj = u_i[min_j];
        int d_i_ij =  0.5 * (d_ij + (ui - uj));
        int d_j_ij = d_ij - d_i_ij;
    
        // 4. Treat (ij) as a new tip, ignoring previous tips i and j
        // 5. Distance of (ij) from other tips k, d_k(ij) = (d_ik+d_jk-d_ij)/2
        DistanceMatrix new_mat(N-1);
        for (int k=0; k<N-1; k++) {
            for (int i=0; i<N-1; i++) {
                for (int j=0; j<i; j++) {
                    int d_ik = mat.get(i, k);
                    int d_jk = mat.get(j, k);
                    int d_ij = mat.get(i, j);
                    new_mat.set(i, k, (d_ik + d_jk - d_ij) / 2);
                }
            }
        }
        mat = new_mat;
        // problem: how to store constructed tree? -> save intermediate value?
    }
    
    return 0
}