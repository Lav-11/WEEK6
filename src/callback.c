#include "../include/callback.h"
#include "../include/cpx_utils.h"



// Callback function
int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle) {
    instance *inst = (instance *)userhandle;
    double *xstar = (double *)malloc(inst->nnodes * sizeof(double));
    if (!xstar) print_error("Memory allocation error");

    double objval = CPX_INFBOUND;
    if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE &&
        CPXcallbackgetcandidatepoint(context, xstar, 0, inst->nnodes - 1, &objval)) {
        print_error("CPXcallbackgetcandidatepoint error");
    }

    if (contextid == CPX_CALLBACKCONTEXT_RELAXATION &&
        CPXcallbackgetrelaxationpoint(context, xstar, 0, inst->nnodes - 1, &objval)) {
        print_error("CPXcallbackgetrelaxationpoint error");
    }

    // Example of using node information
    int mythread = -1;
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
    int mynode = -1;
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
    double incumbent = CPX_INFBOUND;
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);

    // Logic to check infeasibility and add cuts
    int nnz = 0; // Number of non-zero elements in the cut
    double rhs = 0.0; // Right-hand side of the cut
    char sense = 'L'; // Sense of the cut ('L' for <=, 'G' for >=, 'E' for =)
    int *index = NULL; // Indices of the variables in the cut
    double *value = NULL; // Coefficients of the variables in the cut

    // Example: if you find a violated cut
    if (nnz > 0) {
        int izero = 0;
        if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE &&
            CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value)) {
            print_error("CPXcallbackrejectcandidate() error");
        }

        if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) {
            int purgeable = CPX_USECUT_FILTER;
            int local = 0;
            if (CPXcallbackaddusercuts(context, 1, nnz, &rhs, &sense, &izero, index, value, &purgeable, &local)) {
                print_error("CPXcallbackaddusercuts() error");
            }
        }
    }

    free(xstar);
    return 0;
}
