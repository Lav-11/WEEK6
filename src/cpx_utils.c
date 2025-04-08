#include "../include/cpx_utils.h"


void print_error(const char *err) 
{ 
	printf("\n\n ERROR: %s \n\n", err); 
	fflush(NULL); 
	exit(1); 
}   

double dist(int i, int j, instance *inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j]; 
	if ( !inst->integer_costs ) return sqrt(dx*dx+dy*dy);
	int dis = sqrt(dx*dx+dy*dy) + 0.499999999; 					// nearest integer 
	return dis+0.0;
}        


/**************************************************************************************************************************/
int TSPopt(instance *inst)
/**************************************************************************************************************************/
{  

	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if ( error ) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);
	
	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if ( VERBOSE >= 60 ) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);	
	CPXsetdblparam(env, CPX_PARAM_TILIM, 3600.0); 
	// ...

	error = CPXmipopt(env,lp);
	if ( error ) 
	{
		printf("CPX error code %d\n", error);
		print_error("CPXmipopt() error"); 
	}

	// use the optimal solution found by CPLEX
	
	int ncols = CPXgetnumcols(env, lp);
	double *xstar = (double *) calloc(ncols, sizeof(double));
	if ( CPXgetx(env, lp, xstar, 0, ncols-1) ) print_error("CPXgetx() error");	
    for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			if ( xstar[xpos(i,j,inst)] > 0.5 ) 
            if (VERBOSE >= 60) printf("x(%3d,%3d) = 1\n", i+1,j+1);
		}
	}

        // Finds successors
        int *succ = (int *)calloc(inst->nnodes, sizeof(int));
        compute_successors(inst, xstar, succ);
    
        // Prints the successors
        if (VERBOSE >= 60) {
            printf("\nArray of successors:\n");
            for (int i = 0; i < inst->nnodes; i++) {
                printf("    Node %d -> Node %d\n", i + 1, succ[i] + 1);
            }
        }
        // Find and print connected components
        int num_components = find_connected_components(inst->nnodes, xstar, inst);
    
        if (num_components == 1) {
            printf("The graph is fully connected.\n");
        } else {
            printf("The graph has %d connected components.\n", num_components);
        }
    
        plot_graph_to_image(inst->nnodes, inst->xcoord, inst->ycoord, xstar, inst, inst->max_coord, inst->max_coord * 0.1);
    

	free(xstar);
	
	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

	return 0; // or an appropriate nonzero error code

}

/***************************************************************************************************************************/
int xpos(int i, int j, instance *inst)      // to be verified                                           
/***************************************************************************************************************************/
{ 
	if ( i == j ) print_error(" i == j in xpos" );
	if ( i > j ) return xpos(j,i,inst);
	int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;
	return pos;
}
	

/***************************************************************************************************************************/
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{    

	double zero = 0.0;  
	char binary = 'B'; 

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

// add binary var.s x(i,j) for i < j  

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);  		// ... x(1,2), x(1,3) ....
			double obj = dist(i,j,inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, inst) ) print_error(" wrong position for x var.s");
		}
	} 

// add the degree constraints 

	int *index = (int *) calloc(inst->nnodes, sizeof(int));
	double *value = (double *) calloc(inst->nnodes, sizeof(double));

	for ( int h = 0; h < inst->nnodes; h++ )  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h+1);   
		int nnz = 0;
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = xpos(i,h, inst);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");
	} 

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	if ( VERBOSE >= 100 ) CPXwriteprob(env, lp, "model.lp", NULL);   

}

#define DEBUG    // comment out to avoid debugging 
#define EPS 1e-5

/*********************************************************************************************************************************/
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp) // build succ() and comp() wrt xstar()...
/*********************************************************************************************************************************/
{   


	int *degree = (int *) calloc(inst->nnodes, sizeof(int));
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			int k = xpos(i,j,inst);
			if ( fabs(xstar[k]) > EPS && fabs(xstar[k]-1.0) > EPS ) print_error(" wrong xstar in build_sol()");
			if ( xstar[k] > 0.5 ) 
			{
				++degree[i];
				++degree[j];
			}
		}
	}
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		if ( degree[i] != 2 ) print_error("wrong degree in build_sol()");
	}	
	free(degree);


	*ncomp = 0;
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		succ[i] = -1;
		comp[i] = -1;
	}
	
	for ( int start = 0; start < inst->nnodes; start++ )
	{
		if ( comp[start] >= 0 ) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while ( !done )  // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for ( int j = 0; j < inst->nnodes; j++ )
			{
				if ( i != j && xstar[xpos(i,j,inst)] > 0.5 && comp[j] == -1 ) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}	
		succ[i] = start;  // last arc to close the cycle
		
		// go to the next component...
	}
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
    for (int i = 1; i <= comp_id; i++) {
        if (VERBOSE >= 60) printf("Component %d: ", i);
        for (int j = 0; j < nnodes; j++) {
            if (visited[j] == i) {
                if (VERBOSE >= 60) printf("%d ", j + 1); // Node belongs to the component
            }
        }
        if (VERBOSE >= 60) printf("\n");
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
