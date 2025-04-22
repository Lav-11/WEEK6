#include "../include/cpx_utils.h"
#include "../include/callback.h"

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

int TSPopt(instance *inst) {
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    if (error) print_error("CPXopenCPLEX() error");
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
    if (error) print_error("CPXcreateprob() error");

    // Build the TSP model
    build_model(inst, env, lp);
    inst->ncols = CPXgetnumcols(env, lp);

    // Set CPLEX parameters
    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
    if (VERBOSE >= DEBUG) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Enable CPLEX screen output
    CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);
    double timelimit = 600.0;
    CPXsetdblparam(env, CPX_PARAM_TILIM, timelimit);

    double start_time = 0.0;
    CPXgettime(env, &start_time);

    bool apply_patching = false;
    bool use_callback = true;  // Modify this variable to choose whether to use the callback or not

    if (use_callback) {
        // Register the callback
        CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
        if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst))
            print_error("Error while registering the callback");
    }

    // Solve the model
    error = CPXmipopt(env, lp);
    if (error) {
        printf("CPX error code %d\n", error);
        print_error("CPXmipopt() error");
    }

    // Retrieve the solution
    double *xstar = (double *)calloc(inst->ncols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, inst->ncols - 1)) print_error("CPXgetx() error");

    // Build the solution and check for subtours
    int *succ = (int *)calloc(inst->nnodes, sizeof(int));
    int *comp = (int *)calloc(inst->nnodes, sizeof(int));
    int ncomp;
    build_sol(xstar, inst, succ, comp, &ncomp);

    // Perform Benders' loop to eliminate subtours
    benders_loop(inst, env, lp, &xstar, succ, comp, start_time, timelimit, &apply_patching);

    // Apply patching if necessary
    if (apply_patching) {
        printf("Patching solution...\n");
        patch_solution(xstar, inst);
    }

    // Print the solution details if verbose
    if (VERBOSE >= DEBUG) {
        printf("Selected edges:\n");
        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (xstar[xpos(i, j, inst)] > 0.5)
                    printf("x(%3d,%3d) = 1\n", i + 1, j + 1);
            }
        }

        printf("\nArray of successors:\n");
        for (int i = 0; i < inst->nnodes; i++) {
            printf("    Node %d -> Node %d\n", i + 1, succ[i] + 1);
        }
    }

    // Rebuild the solution and print the number of components
    build_sol(xstar, inst, succ, comp, &ncomp);
    printf("Number of components:  %d \n", ncomp);

    // Get and print the final objective value
    double objval;
    if (CPXgetobjval(env, lp, &objval)) {
        print_error("CPXgetobjval() error");
    }
    printf("Final objective value: %f\n", objval);

    // Plot the graph and save it as an image
    plot_graph_to_image(inst->nnodes, inst->xcoord, inst->ycoord, xstar, inst, inst->max_coord, inst->max_coord * 0.1);

    // Free allocated memory and close CPLEX
    free(comp);
    free(succ);
    free(xstar);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return 0;
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


void benders_loop(instance *inst, CPXENVptr env, CPXLPptr lp, double **xstar_ptr,
                  int *succ, int *comp, double start_time, double timelimit,
                  bool *apply_patching) {
    
    int ncols = CPXgetnumcols(env, lp);
    double *xstar = (double *)calloc(ncols, sizeof(double));
    if (xstar == NULL) print_error("calloc() error for xstar");

    // Solve the initial model
    if (CPXmipopt(env, lp)) print_error("CPXmipopt() error (initial)");
    if (CPXgetx(env, lp, xstar, 0, ncols - 1)) print_error("CPXgetx() error (initial)");

    int ncomp = -1;
    build_sol(xstar, inst, succ, comp, &ncomp);

    // Loop until there are no subtours
    while (ncomp > 1) {
        double end_time = 0.0;
        CPXgettime(env, &end_time);
        double elapsed_time = end_time - start_time;

        printf("Subtours detected: %d components. Adding SEC constraints...\n", ncomp);

        // Add Subtour Elimination Constraints (SEC)
        add_SEC_constraints(inst, env, lp, xstar, NULL, -1);

        free(xstar);
        xstar = (double *)calloc(ncols, sizeof(double));
        if (xstar == NULL) print_error("calloc() error for xstar (loop)");

        // Re-solve the model after adding SEC
        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error (after SEC)");
        if (CPXgetx(env, lp, xstar, 0, ncols - 1)) print_error("CPXgetx() error (after SEC)");

        build_sol(xstar, inst, succ, comp, &ncomp);

        // Check if the time limit has been exceeded
        if (elapsed_time >= timelimit) {
            printf("Time limit exceeded. Initiating patching process...\n");
            *apply_patching = true;
            break;
        }
    }

    *xstar_ptr = xstar;  // Pass the result back
}

void add_SEC_constraints(instance *inst, CPXENVptr env, CPXLPptr lp, double *xstar,
                         CPXCALLBACKCONTEXTptr context, int contextid) {
    int *succ = (int *)calloc(inst->nnodes, sizeof(int));
    int *comp = (int *)calloc(inst->nnodes, sizeof(int));
    int ncomp;

    build_sol(xstar, inst, succ, comp, &ncomp);

    // If there is only one component, no SEC is needed
    if (ncomp <= 1) {
        free(succ);
        free(comp);
        return;
    }

    printf("Adding %d SEC constraints%s...\n", ncomp, context ? " via callback" : "");

    // Add SEC constraints for each component
    for (int k = 1; k <= ncomp; k++) {
        int size = 0;
        for (int i = 0; i < inst->nnodes; i++)
            if (comp[i] == k) size++;
        if (size <= 1) continue;

        int max_nz = size * (size - 1) / 2;
        int *index = (int *)malloc(max_nz * sizeof(int));
        double *value = (double *)malloc(max_nz * sizeof(double));
        int nz = 0;

        // Collect edges within the component
        for (int i = 0; i < inst->nnodes; i++) {
            if (comp[i] != k) continue;
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (comp[j] != k) continue;
                index[nz] = xpos(i, j, inst);
                value[nz++] = 1.0;
            }
        }

        double rhs = size - 1.0;
        char sense = 'L';
        int izero = 0;

        if (context) {
            // Add SEC via callback
            if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
                if (CPXcallbackrejectcandidate(context, 1, nz, &rhs, &sense, &izero, index, value)) {
                    print_error("CPXcallbackrejectcandidate() error");
                }
            } else if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) {
                int purgeable = CPX_USECUT_FILTER;
                int local = 0;
                if (CPXcallbackaddusercuts(context, 1, nz, &rhs, &sense, &izero, index, value, &purgeable, &local)) {
                    print_error("CPXcallbackaddusercuts() error");
                }
            }
        } else {
            // Add SEC in the standard phase
            char *rowname = (char *)malloc(100 * sizeof(char));
            sprintf(rowname, "SEC_comp_%d", k);
            if (CPXaddrows(env, lp, 0, 1, nz, &rhs, &sense, &izero, index, value, NULL, &rowname)) {
                print_error("CPXaddrows() failed for SEC");
            }
            free(rowname);
        }

        free(index);
        free(value);
    }

    free(succ);
    free(comp);
}

void invert_path(int start, int end, int *succ, double *xstar, instance *inst) {
    int current = start;
    int prev = -1;

    // Follow the path and invert the successors
    while (current != end) {
        int next = succ[current]; // Save the current successor
        succ[current] = prev;    // Invert the successor
        prev = current;          // Update the previous node
        current = next;          // Move to the next node
    }

    // Update the last node
    succ[current] = prev;

    // Update xstar to reflect the inversion
    current = start;
    while (current != end) {
        int next = succ[current];
        xstar[xpos(current, next, inst)] = 1.0;  // Add the inverted edge
        xstar[xpos(next, current, inst)] = 0.0; // Remove the original edge
        current = next;
    }
}

void patch_solution(double *xstar, instance *inst) {
    int *succ = (int *) malloc(inst->nnodes * sizeof(int));
    int *comp = (int *) malloc(inst->nnodes * sizeof(int));
    int ncomp;

    build_sol(xstar, inst, succ, comp, &ncomp);

    // If there is only one component, no patching is needed
    if (ncomp <= 1) {
        printf("No patching needed. Solution is a single tour.\n");
        free(succ);
        free(comp);
        return;
    }

    printf("Patching needed. Components found: %d\n", ncomp);

    while (ncomp > 1) {

        double best_cost = DBL_MAX;
        int best_i = -1, best_j = -1;
        int succ_i = -1, succ_j = -1;
        bool swap = false;

        // Find the best pair of edges to connect components
        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (comp[i] != comp[j]) {
                    swap = false; // Reset swap flag
                    // Find successors
                    succ_i = succ[best_i];
                    succ_j = succ[best_j];
                    double cij = dist(i, j, inst) + dist(succ_j, succ_i, inst); 
                    if (cij < best_cost) {
                        best_cost = cij;
                        best_i = i;
                        best_j = j;
                    }
                    double cij_swap = dist(i, succ_j, inst) + dist(j, succ_i, inst);
                    if (cij_swap < best_cost) {
                        best_cost = cij_swap;
                        best_i = i;
                        best_j = j;
                        swap = true; // Indicate that we are swapping
                    }
                }
            }
        }

        if (best_i == -1 || best_j == -1) {
            printf("No valid edge pairs found. Aborting patching.\n");
            break;
        }

        // Remove two edges: (best_i, succ_i) and (best_j, succ_j)
        xstar[xpos(best_i, succ_i, inst)] = 0.0;
        xstar[xpos(best_j, succ_j, inst)] = 0.0;

        if  (swap == false){

            // Add two new edges: (best_i, best_j) and (succ_i, succ_j)
            xstar[xpos(best_i, best_j, inst)] = 1.0;
            xstar[xpos(succ_j, succ_i, inst)] = 1.0;
        }
        else if (swap == true) { // Invert the path
            printf("Inverting path from %d to %d\n", best_j, succ_j);
        
            // Invert path from best_j to succ_j
            invert_path(best_j, succ_j, succ, xstar, inst);
        
            // Add new edges
            xstar[xpos(best_i, succ_j, inst)] = 1.0;
            xstar[xpos(best_j, succ_i, inst)] = 1.0;
        }

        // Rebuild the solution
        build_sol(xstar, inst, succ, comp, &ncomp);
        printf("Components after patching: %d\n", ncomp);
    }

    free(succ);
    free(comp);
}