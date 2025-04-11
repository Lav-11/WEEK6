#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <ilcplex/cplex.h>

#define VERBOSE 1
#define DEBUG 80
#define EPS 1e-5

// Define instance structure
typedef struct {
    int nnodes;
    double *xcoord;
    double *ycoord;
    int integer_costs;
    double zbest;
    double max_coord;
} instance;

// Function prototypes

// Error handling function
void print_error(const char *err);

// Compute distance between two nodes
double dist(int i, int j, instance *inst);

// Compute position index for x(i, j)
int xpos(int i, int j, instance *inst);

// Build the CPLEX model
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);

// Solve the TSP using CPLEX
int TSPopt(instance *inst);

// Compute the successors of each node in the solution
void compute_successors(instance *inst, double *xstar, int *succ);

// Function to compute the successors of each node in the solution
void dfs(int node, int *visited, int **adj_matrix, int n, int comp_id);

// Function to find and count connected components
int find_connected_components(int nnodes, double *xstar, instance *inst);

// Function to ensure the directory exists
void create_directory(const char *path);

// Function to plot the graph and save it as an image
void plot_graph_to_image(int nnodes, double *xcoord, double *ycoord, double *xstar, instance *inst, double max_coord, double padding);

// Function to build the solution from xstar
void add_SEC_constraints(instance *inst, CPXENVptr env, CPXLPptr lp, double *xstar);

void patch_solution(double *xstar, instance *inst);