#include "../include/callback.h"
#include "../include/cpx_utils.h"

// Callback function
int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle) {
    instance *inst = (instance *)userhandle;

    // Get the number of variables (edges) from the model
    int ncols = inst->ncols;

    // Allocate the xstar vector for the variables
    double *xstar = (double *)malloc(ncols * sizeof(double));
    if (!xstar) print_error("Memory allocation error");

    double objval = CPX_INFBOUND;
    // Extract the current solution from the callback (candidate or relaxation)
    if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
        if (CPXcallbackgetcandidatepoint(context, xstar, 0, ncols - 1, &objval)) {
            print_error("CPXcallbackgetcandidatepoint error");
        }
    } else if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) {
        if (CPXcallbackgetrelaxationpoint(context, xstar, 0, ncols - 1, &objval)) {
            print_error("CPXcallbackgetrelaxationpoint error");
        }
    } else {
        free(xstar);
        return 0; // We are not interested in handling other types of callbacks
    }

    // Optional debug information
    int mythread = -1;
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
    int mynode = -1;
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
    double incumbent = CPX_INFBOUND;
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);

    // Add SEC constraints if there are subtours
    add_SEC_constraints(inst, NULL, NULL, xstar, context, contextid);

    free(xstar);
    return 0;
}