#include "../include/callback.h"
#include "../include/cpx_utils.h"

// Callback function
int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle) {
    instance *inst = (instance *)userhandle;

    // Ottieni numero di variabili (archi) dal modello
    int ncols = inst->ncols;

    // Alloca il vettore xstar delle variabili
    double *xstar = (double *)malloc(ncols * sizeof(double));
    if (!xstar) print_error("Memory allocation error");

    double objval = CPX_INFBOUND;
    // Estrai la soluzione corrente dal callback (candidate o relaxation)
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
        return 0; // Non ci interessa gestire altri tipi di callback
    }

    // Informazioni di debug opzionali
    int mythread = -1;
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
    int mynode = -1;
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
    double incumbent = CPX_INFBOUND;
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);

    // Aggiunta vincoli SEC se ci sono sottogiri
    add_SEC_constraints(inst, NULL, NULL, xstar, context, contextid);

    free(xstar);
    return 0;
}