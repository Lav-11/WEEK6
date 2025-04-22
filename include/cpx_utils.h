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
    int ncols;
    double *xcoord;
    double *ycoord;
    int integer_costs;
    double zbest;
    double max_coord;
    int use_callback;  // 1 to use callback, 0 otherwise
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

void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);

// Solve the TSP using CPLEX
int TSPopt(instance *inst);

// Function to ensure the directory exists
void create_directory(const char *path);

// Function to plot the graph and save it as an image
void plot_graph_to_image(int nnodes, double *xcoord, double *ycoord, double *xstar, instance *inst, double max_coord, double padding);

// Function to use benders
void benders_loop(instance *inst, CPXENVptr env, CPXLPptr lp, double **xstar_ptr, int *succ, int *comp, double start_time, double timelimit, bool *apply_patching);

// Function to build the solution from xstar
void add_SEC_constraints(instance *inst, CPXENVptr env, CPXLPptr lp, double *xstar,
CPXCALLBACKCONTEXTptr context, int contextid);

void patch_solution(double *xstar, instance *inst);

void invert_path(int start, int end, int *succ, double *xstar, instance *inst);

