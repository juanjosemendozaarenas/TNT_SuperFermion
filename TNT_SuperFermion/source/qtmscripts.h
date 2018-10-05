/*
Authors: Juan Jose Mendoza-Arenas and Stephen Clark
Date:    August 2018
(c) Universidad de los Andes 2018
*/

/*! This header file contains prototypes for all the functions for the time evolution
    of an open quantum system analyzed with the super-fermion approach
 */

#define tntMpsExOp tntExOp

void tntMpsExpecOutputSuperFermion(tntNetwork mps,      /* The network representing the MPS. Unchanged by the function. */
                                   tntNetwork left_vacuum, /* Left vacuum of the super-fermion method for open systems */
                                   tntMpsExOp *Op,      /* The operators for calculating expectation values */
                                   unsigned printOutput,/* Set to 1 to print expectation values to screen, 0 otherwise */
                                   unsigned saveOutput, /* Set to 1 to save values to output file, 0 otherwise */
                                   char *savename,      /* Path to the output file. If the extension is missing ".mat" is added. Only used if saveOutput is 1. */
                                   unsigned counter);   /* can be used if saving multiple batches of expectation values e.g. for different timesteps.
                                                         It appends the values to a pre-exisiting array in the position given by counter.
                                                         Pass zero if a counter is not required. */

tntNetwork tntMpsCreatePropST2scPartialStar(unsigned L,                /* Length of the system. */
                                            unsigned N,                /* Number of lead modes. */
                                            tntComplex dtc,            /* Size of the time step. See the main description for information on how real and imaginary parts are applied */
                                            tntNodeArray *nnL,         /* Array of nearest-neighbour operators for left site. Send NULL if there are no nearest neighbour terms. */
                                            tntNodeArray *nnR,         /* Array of nearest-neighbour operators for right site. Send NULL if there are no nearest neighbour terms. */
                                            tntComplexArray *nnparam,  /* Array of parameters for nearest-neighbour operators. Send NULL if there are no nearest neighbour terms. */
                                            tntNodeArray *os,          /* Array of onsite operators. Send NULL if there are no on-site operators. */
                                            tntComplexArray *osparam,  /* Parameters for the on-site operators. Send NULL if there are no on-site operators. */
                                            tntNodeArray swap_gate);   /* Swap gate, two swap last site of the system across the lead */

tntNetwork tntMpsPropArrayToST2scPartialStar(unsigned L, /* Length of the system. */
                                             unsigned N, /* Number of lead modes. */
                                             tntNodeArray Proparr, /* Array of propagators. Uncahnged by the function - copies of all nodes are used. */
                                             tntNodeArray swap_gate); /* Swap gate, two swap last site of the system across the lead */

tntNodeArray tntMpsCreatePropArrayPartialStar(unsigned L, /* Length of system. */
                                              unsigned N, /* Number of lead modes. */
                                              tntComplex h, /* Uniform scale factor to apply to all terms. See the main description for information on how real and imaginary parts are applied */
                                              tntNodeArray *nnl, /* Array of nearest-neighbour operators for left site. Send NULL if there are no nearest neighbour terms. Unchanged by function - copies are used. */
                                              tntNodeArray *nnr, /* Array of nearest-neighbour operators for right site. Send NULL if there are no nearest neighbour terms. Unchanged by function - copies are used. */
                                              tntComplexArray *nnparam, /* Array of parameters for nearest-neighbour operators. Send NULL if there are no nearest neighbour terms. */
                                              tntNodeArray *os,  /* Array of onsite operators. Send NULL if there are no on-site operators. Unchanged by function - copies are used. */
                                              tntComplexArray *osparam);  /* Parameters for the on-site operators. Send NULL if there are no on-site operators. */
