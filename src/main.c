#include "../include/cpx_utils.h"



int main() {
    // Define problem parameters
    int nnodes = 50; // Number of nodes
    double max_coord = 100.0; // Upper limit for node coordinates
    instance inst;
    inst.nnodes = nnodes;
    inst.integer_costs = 1; // Use integer costs

    // Allocate memory for node coordinates
    inst.xcoord = (double *)malloc(nnodes * sizeof(double));
    inst.ycoord = (double *)malloc(nnodes * sizeof(double));
    if (!inst.xcoord || !inst.ycoord) {
        print_error("Memory allocation failed.");
    }

    // Generate random coordinates for the nodes
    srand(42); // Seed for reproducibility
    for (int i = 0; i < nnodes; i++) {
        inst.xcoord[i] = rand() % (int) max_coord; // Random x-coordinates
        inst.ycoord[i] = rand() % (int) max_coord; // Random y-coordinates
        printf("Node %d: (%.2f, %.2f)\n", i + 1, inst.xcoord[i], inst.ycoord[i]);
    }

    // Initialize CPLEX environment
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    if (error) print_error("Could not open CPLEX environment.");
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP Model");
    if (error) print_error("Could not create CPLEX problem.");

    // Define variables
    int nvars = nnodes * (nnodes - 1) / 2; // Variables for each edge in the graph
    double *obj = (double *)malloc(nvars * sizeof(double));   // Objective coefficients
    double *lb = (double *)malloc(nvars * sizeof(double));    // Lower bounds
    double *ub = (double *)malloc(nvars * sizeof(double));    // Upper bounds
    char *ctype = (char *)malloc(nvars * sizeof(char));       // Variable types
    char **colnames = (char **)malloc(nvars * sizeof(char *));
    int *indices = (int *)malloc(nvars * sizeof(int));

    int var_index = 0;
    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            obj[var_index] = dist(i, j, &inst); // Distance as cost
            lb[var_index] = 0.0;
            ub[var_index] = 1.0;
            ctype[var_index] = 'I'; // Integer variable
            colnames[var_index] = (char *)malloc(50 * sizeof(char));
            sprintf(colnames[var_index], "x(%d,%d)", i + 1, j + 1);
            indices[var_index] = var_index;
            var_index++;
        }
    }

    // Add variables to the model
    if (CPXaddcols(env, lp, nvars, 0, obj, NULL, NULL, NULL, lb, ub, colnames)) {
        print_error("Failed to add variables.");
    }

    // Set variables as integer explicitly
    if (CPXchgctype(env, lp, nvars, indices, ctype)) {
        print_error("Failed to set variable types.");
    }

    // Degree constraints: each node must have degree 2
    for (int i = 0; i < nnodes; i++) {
        int *rmatind = (int *)malloc((nnodes - 1) * sizeof(int));
        double *rmatval = (double *)malloc((nnodes - 1) * sizeof(double));
        int rmatbeg[1] = {0}; // Single constraint starts at position 0
        int count = 0;

        for (int j = 0; j < nnodes; j++) {
            if (i != j) {
                int idx = xpos(i, j, &inst); // Index of variable
                if (idx < 0 || idx >= nvars) {
                    printf("Invalid index for variable (%d, %d): %d\n", i + 1, j + 1, idx);
                    continue;
                }
                rmatind[count] = idx; // Variable index
                rmatval[count] = 1.0; // Coefficient in the constraint
                count++;
            }
        }

        double rhs = 2.0; // Degree of each node must equal 2
        char sense = 'E'; // Equality constraint
        if (CPXaddrows(env, lp, 0, 1, count, &rhs, &sense, rmatbeg, rmatind, rmatval, NULL, NULL)) {
            print_error("Failed to add degree constraints.");
        }

        free(rmatind);
        free(rmatval);
    }

    // Solve the problem
    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Enable screen output
    CPXsetdblparam(env, CPX_PARAM_TILIM, 3600.0);  // Set time limit

    error = CPXmipopt(env, lp);
    if (error) print_error("Failed to optimize the problem.");

    // Retrieve the solution
    int ncols = CPXgetnumcols(env, lp);
    double *xstar = (double *)calloc(ncols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, ncols - 1)) {
        print_error("Failed to retrieve solution.");
    }

    printf("Optimal solution:\n");
    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            if (xstar[xpos(i, j, &inst)] > 0.5) {
                printf("    Edge (%d, %d) selected.\n", i + 1, j + 1);
            }
        }
    }

    // Finds successors
    int *succ = (int *)calloc(nnodes, sizeof(int));
    compute_successors(&inst, xstar, succ);

    // Prints the successors
    printf("\nArray of successors:\n");
    for (int i = 0; i < nnodes; i++) {
        printf("    Node %d -> Node %d\n", i + 1, succ[i] + 1);
    }

    // Find and print connected components
    int num_components = find_connected_components(nnodes, xstar, &inst);

    if (num_components == 1) {
        printf("The graph is fully connected.\n");
    } else {
        printf("The graph has %d connected components.\n", num_components);
    }

    plot_graph_to_image(nnodes, inst.xcoord, inst.ycoord, xstar, &inst, max_coord, 5.0);



    // Clean up
    free(xstar);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    for (int i = 0; i < nvars; i++) free(colnames[i]);
    free(colnames);
    free(obj);
    free(lb);
    free(ub);
    free(ctype);
    free(indices);
    free(inst.xcoord);
    free(inst.ycoord);

    return 0;
}








