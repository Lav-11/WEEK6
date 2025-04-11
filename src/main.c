#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../include/cpx_utils.h"

int main() {
    // Create a random TSP instance
    instance inst;
    inst.nnodes = 200; // number of nodes
    inst.xcoord = (double *)malloc(inst.nnodes * sizeof(double));
    inst.ycoord = (double *)malloc(inst.nnodes * sizeof(double));
    inst.integer_costs = 0; // Use continuous costs (Euclidean distances)
    inst.max_coord = 100.0; // Limit for coordinates (for random generation)

    // Initialize random coordinates for the nodes with a fixed seed for reproducibility
    srand(42); // Fixed seed for repeatability
    for (int i = 0; i < inst.nnodes; i++) {
        inst.xcoord[i] = (rand() % (int)inst.max_coord) + 1.0;
        inst.ycoord[i] = (rand() % (int)inst.max_coord) + 1.0;
    }

    // Solve the TSP problem with CPLEX
    TSPopt(&inst);

    // Free memory
    free(inst.xcoord);
    free(inst.ycoord);

    return 0;
}