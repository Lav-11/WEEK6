#include "../include/cpx_utils.h"


// Error handling function
void print_error(const char *err) {
    printf("\n\n ERROR: %s \n\n", err);
    fflush(NULL);
    exit(1);
}

// Compute distance between two nodes
double dist(int i, int j, instance *inst) {
    double dx = inst->xcoord[i] - inst->xcoord[j];
    double dy = inst->ycoord[i] - inst->ycoord[j];
    if (!inst->integer_costs) return sqrt(dx * dx + dy * dy);
    int dis = sqrt(dx * dx + dy * dy) + 0.499999999;
    return dis + 0.0;
}


// Solve the TSP using CPLEX
int TSPopt(instance *inst) {
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    if (error) print_error("CPXopenCPLEX() error");
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP Model");
    if (error) print_error("CPXcreateprob() error");

    build_model(inst, env, lp);
    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
    CPXsetdblparam(env, CPX_PARAM_TILIM, 3600.0);

    error = CPXmipopt(env, lp);
    if (error) print_error("CPXmipopt() error");

    int ncols = CPXgetnumcols(env, lp);
    double *xstar = (double *)calloc(ncols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, ncols - 1)) print_error("CPXgetx() error");
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = i + 1; j < inst->nnodes; j++) {
            if (xstar[xpos(i, j, inst)] > 0.5)
                printf("  ... x(%3d,%3d) = 1\n", i + 1, j + 1);
        }
    }
    free(xstar);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return 0;
}

// Compute position index for x(i, j)
int xpos(int i, int j, instance *inst) {
    if (i == j) print_error("i == j in xpos");
    if (i > j) return xpos(j, i, inst);
    return i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;
}

// Build the CPLEX model
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp) {
    double zero = 0.0;
    char binary = 'B';
    char **cname = (char **)calloc(1, sizeof(char *));
    cname[0] = (char *)calloc(100, sizeof(char));

    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = i + 1; j < inst->nnodes; j++) {
            sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
            double obj = dist(i, j, inst);
            double lb = 0.0, ub = 1.0;
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
                print_error("CPXnewcols() error");
            if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst))
                print_error("Incorrect x variable position");
        }
    }

    free(cname[0]);
    free(cname);
}

// Compute the successors of each node in the solution
void compute_successors(instance *inst, double *xstar, int *succ) {
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = 0; j < inst->nnodes; j++) {
            if (i != j && xstar[xpos(i, j, inst)] > 0.5) {
                succ[i] = j; // Node j is successor of node i
                break;
            }
        }
    }
}

// Depth-First Search (DFS) to explore nodes
void dfs(int node, int *visited, int **adj_matrix, int n, int comp_id) {
    visited[node] = comp_id; // Mark the node with the connected component ID
    for (int j = 0; j < n; j++) {
        if (adj_matrix[node][j] == 1 && visited[j] == 0) { // Adjacent unvisited node
            dfs(j, visited, adj_matrix, n, comp_id); // Recursive call
        }
    }
}

// Function to find and count connected components
int find_connected_components(int nnodes, double *xstar, instance *inst) {
    // Create an adjacency matrix
    int **adj_matrix = (int **)malloc(nnodes * sizeof(int *));
    for (int i = 0; i < nnodes; i++) {
        adj_matrix[i] = (int *)calloc(nnodes, sizeof(int)); // Initialize to 0
    }

    // Populate the adjacency matrix based on edges selected in xstar
    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            if (xstar[xpos(i, j, inst)] > 0.5) { // If edge (i, j) is selected
                adj_matrix[i][j] = 1;
                adj_matrix[j][i] = 1; // Undirected graph
            }
        }
    }

    // Array to track visited nodes
    int *visited = (int *)calloc(nnodes, sizeof(int));

    // Count connected components
    int comp_id = 0;
    for (int i = 0; i < nnodes; i++) {
        if (visited[i] == 0) { // Unvisited node
            comp_id++;
            dfs(i, visited, adj_matrix, nnodes, comp_id); // Explore the component
        }
    }

    // Print the results
    printf("\nNumber of connected components: %d\n", comp_id);
    for (int i = 1; i <= comp_id; i++) {
        printf("Component %d: ", i);
        for (int j = 0; j < nnodes; j++) {
            if (visited[j] == i) {
                printf("%d ", j + 1); // Node belongs to the component
            }
        }
        printf("\n");
    }

    // Free memory
    for (int i = 0; i < nnodes; i++) {
        free(adj_matrix[i]);
    }
    free(adj_matrix);
    free(visited);

    return comp_id; // Return the number of connected components
}


// Function to ensure the directory exists
void create_directory(const char *path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        #ifdef _WIN32
        mkdir(path);
        #else
        mkdir(path, 0700);
        #endif
    }
}

void plot_graph_to_image(int nnodes, double *xcoord, double *ycoord, double *xstar, instance *inst, double max_coord, double padding) {
    // Ensure the ../data directory exists
    const char *dir_path = "../data";
    create_directory(dir_path);

    // Open a file to write node coordinates and edges in ../data/graph_data.dat
    char data_file_path[256];
    snprintf(data_file_path, sizeof(data_file_path), "%s/graph_data.dat", dir_path);
    FILE *file = fopen(data_file_path, "w");
    if (!file) {
        printf("Error: Unable to open file for writing.\n");
        return;
    }

    // Write the selected edges
    fprintf(file, "# Selected edges (lines)\n");
    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            if (xstar[xpos(i, j, inst)] > 0.5) { // If edge (i, j) is selected
                fprintf(file, "%.2f %.2f\n", xcoord[i], ycoord[i]); // Start of edge
                fprintf(file, "%.2f %.2f\n\n", xcoord[j], ycoord[j]); // End of edge
            }
        }
    }

    // Write the node coordinates
    fprintf(file, "\n\n# Node coordinates (points)\n");
    for (int i = 0; i < nnodes; i++) {
        fprintf(file, "%.2f %.2f\n", xcoord[i], ycoord[i]); // X Y format
    }

    fclose(file);

    // Compute dynamic bounds for the plot
    double xmin = 0 - padding;
    double xmax = max_coord + padding;
    double ymin = 0 - padding;
    double ymax = max_coord + padding;

    // Call GNUplot to save the graph as an image
    FILE *gnuplot = popen("gnuplot", "w");
    if (!gnuplot) {
        printf("Error: Unable to open GNUplot.\n");
        return;
    }

    // Image file path
    char image_file_path[256];
    snprintf(image_file_path, sizeof(image_file_path), "%s/graph_plot.png", dir_path);

    // Setup GNUplot commands for saving the image
    fprintf(gnuplot, "set terminal png size 800,600\n"); // PNG format, resolution 800x600
    fprintf(gnuplot, "set output '%s'\n", image_file_path);
    fprintf(gnuplot, "set title 'Graph Visualization'\n");
    fprintf(gnuplot, "set grid\n");
    fprintf(gnuplot, "set key off\n");
    fprintf(gnuplot, "set size square\n");
    fprintf(gnuplot, "set xrange [%f:%f]\n", xmin, xmax); // Dynamically set X range
    fprintf(gnuplot, "set yrange [%f:%f]\n", ymin, ymax); // Dynamically set Y range

    // Plot the edges (red) and the nodes (blue)
    fprintf(gnuplot, "plot \\\n");
    fprintf(gnuplot, "  '%s' index 0 using 1:2 with lines lc rgb 'red' title 'Edges', \\\n", data_file_path);
    fprintf(gnuplot, "  '%s' index 1 using 1:2 with points pt 7 lc rgb 'blue' title 'Nodes'\n", data_file_path);

    pclose(gnuplot);

    printf("Graph saved as an image: %s\n", image_file_path);
}
