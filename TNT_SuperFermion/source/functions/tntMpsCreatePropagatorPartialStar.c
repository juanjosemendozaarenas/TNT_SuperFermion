/*
Authors: Juan Jose Mendoza Aenas and Stephen Clark

(c) Universidad de los Andes 2018
*/

/* Include the headers for the TNT library */
#include "tnt.h"
#include "tntMps.h"

/* Include the header for all network functions included in djscripts */
#include "../qtmscripts.h"

/*! Based on the TNT function tntMpsCreatePropagator.c, but reduced to the situation requred here, namely
 * that there are on-site terms and that the parameters are inhomogeneous
 *
 * Creates the array of two site gates required for evolution under the effective superfermion Hamiltonian
 * for a system + lead + environment setup.
 */
tntNodeArray tntMpsCreatePropArrayPartialStar(unsigned L, /*!< Length of system. */
                                   unsigned N, /*!< Number of lead modes. */
                                   tntComplex h, /*!< Uniform scale factor to apply to all terms. See the main description for information on how real and imaginary parts are applied */
                                   tntNodeArray *nnl, /*!< Array of nearest-neighbour operators for left site. Send NULL if there are no nearest neighbour terms. Unchanged by function - copies are used. */
                                   tntNodeArray *nnr, /*!< Array of nearest-neighbour operators for right site. Send NULL if there are no nearest neighbour terms. Unchanged by function - copies are used. */
                                   tntComplexArray *nnparam, /*!< Array of parameters for nearest-neighbour operators. Send NULL if there are no nearest neighbour terms. */
                                   tntNodeArray *os,  /*!< Array of onsite operators. Send NULL if there are no on-site operators. Unchanged by function - copies are used. */
                                   tntComplexArray *osparam)  /*!< Parameters for the on-site operators. Send NULL if there are no on-site operators. */
{
    unsigned M;                         /* Total number of sites */
    tntNodeArray Parr;                  /* Array of propagators to be returned by function */
    tntComplex hs;                      /* Uniform scale factor multiplied by -i for real part, i for imaginary part */
    unsigned i, j;                      /* Used for looping over terms and sites respectively */
    tntComplexArray two_site_params;    /* Array for parameters on each site */

    /* Check that there are the same number of left and right nearest neighbour operators */
    if ((NULL == nnl && NULL != nnr) || (NULL != nnl && NULL == nnr) || ((NULL != nnl && NULL != nnr) && nnl->sz != nnr->sz)) 
        tntErrorPrint("Cannot create a propagator!|The number of left nearest-neighbour terms is not equal to the number of right nearest-neighbour terms!");  /* NO_COVERAGE */

    /* Check that there is a parameter for each of the left and right operators, or for each of the operators and each of the sites */
    if (((nnl != NULL)&&(nnl->sz != 0))&&(NULL == nnparam || nnl->sz != nnparam->numrows)) {
        tntErrorPrint("Cannot create a propagator!|The number of nearest neighbour parameters (%d) is not equal to the number of nearest neighbour terms (%d)!",nnparam->numrows,nnl->sz);
    } /* NO_COVERAGE */

    /* Check that there is a parameter for each of the onsite operators */
    if (((os != NULL)&&(os->sz != 0)) && (NULL == osparam || os->sz != osparam->numrows)) {
        tntErrorPrint("Cannot create a propagator!|The number of on-site parameters per site (%d) is not equal to the number of on-site terms (%d)!",osparam->numrows,os->sz);
    } /* NO_COVERAGE */

    /* Check that the time step is not non-negligible */
    if (TNT_DEFAULT_TOL > fabs(h.re) && TNT_DEFAULT_TOL > fabs(h.im)) 
        tntErrorPrint("Error! Cannot create spin half two-site gates as the time step size is less than the tolerance.");  /* NO_COVERAGE */

    /* Check that there is a non-zero number of operators */
    if ((NULL == nnl || 0 == nnl->sz) && (NULL == os || 0 == os->sz)) 
        tntErrorPrint("Cannot create a propogator, as there are no terms!");  /* NO_COVERAGE */

    /* Turn off warning since we know we are going to strip QN */
    tntSysQNClearWarnOff();
    
    /* Define total number of sites */
    M = L + N;
    
    /* Allocate the array to be returned */
    Parr = tntNodeArrayAlloc(M-1);
    
    /* First take real part of time step and multiply by -i and take the imaginary part of the timestep and multiply by i */
    hs.re = -h.im;
    hs.im = -h.re;

    /* Create the identity node */
    tntNode eye = tntNodeCreateEyeOp(os->vals[0]); /* Identity node created using basis operator*/
    tntNodeArray Lterms, Rterms; /* Arrays that contain both onsite and nearest neighbour left and right terms */
    
    /* Create a node array that contains the nearest neighbour terms and the onsite terms */
    Lterms = tntNodeArrayAlloc(nnl->sz + os->sz*2); /* twice the number of onsite terms, as require one term for .. */
    Rterms = tntNodeArrayAlloc(nnl->sz + os->sz*2); /* ... operator on left site, and one term for operator on right site */
    
    /* Create a complex array for the parameters */
    two_site_params = tntComplexArrayAlloc(Lterms.sz); /* Array for parameters on each site */

    /* First put the nearest neighbour terms in the array */
    for (i = 0; i < nnl->sz; i++) {
        Lterms.vals[i] = tntNodeCopy(nnl->vals[i]);
        Rterms.vals[i] = tntNodeCopy(nnr->vals[i]);
    }
    
    /* Then put the on-site terms in the array */
    for (i = 0; i < os->sz; i++) {
        /* For the first term put the on-site operator on the left */
        Lterms.vals[i*2 + nnl->sz] = tntNodeCopy(os->vals[i]);
        Rterms.vals[i*2 + nnl->sz] = tntNodeCopy(eye);

        /* For the first term put the on-site operator on the right */
        Rterms.vals[i*2 + 1 + nnl->sz] = tntNodeCopy(os->vals[i]);
        Lterms.vals[i*2 + 1 + nnl->sz] = tntNodeCopy(eye);
    }
    
    /* Every node in the array will be different */
    double scL, scR; /* Amount to scale onsite parameter by - depends on site */
    unsigned ind, indL, indR; /* index in array */

    /* Set the terms for all sites */
    for (j = 0; j<M-1; j++) {
        /* Set the nearest neighbour terms for the central sites */
        for (i = 0; i < nnl->sz; i++) {
            ind = (1 == nnparam->numcols) ? i : i + j*nnl->sz;
            two_site_params.vals[i] = nnparam->vals[ind];
        }

        /* Define scaling factors */
        scL = (0==j) ? 1.0 : 2.0;
        scR = 2.0; /* For the usual Trotter decomposition this is "scR = (M-2 == j) ? 1.0 : 2.0". However here there is no difference between the final lead site and the rest, as all of them are involved in two gates */
            
        /* Set the onsite terms for the central sites */
        if(j<L-1) { /* Do the usual for system degrees of freedom */
            for (i = 0; i < os->sz; i++) {
                indL = (1 == osparam->numcols) ? i : i + j*os->sz;
                indR = (1 == osparam->numcols) ? i : i + (j+1)*os->sz;
                two_site_params.vals[i*2 + nnl->sz].re = osparam->vals[indL].re/scL;
                two_site_params.vals[i*2 + nnl->sz].im = osparam->vals[indL].im/scL;
                two_site_params.vals[i*2 + 1 + nnl->sz].re = osparam->vals[indR].re/scR;
                two_site_params.vals[i*2 + 1 + nnl->sz].im = osparam->vals[indR].im/scR;
            }
        } else { /* Sites on the right appear on the sweep only twice, so their total value is taken */
            for (i = 0; i < os->sz; i++) {
                indL = (1 == osparam->numcols) ? i : i + j*os->sz;
                indR = (1 == osparam->numcols) ? i : i + (j+1)*os->sz;
                two_site_params.vals[i*2 + nnl->sz].re = 0;
                two_site_params.vals[i*2 + nnl->sz].im = 0;
                two_site_params.vals[i*2 + 1 + nnl->sz].re = 2*osparam->vals[indR].re/scR;
                two_site_params.vals[i*2 + 1 + nnl->sz].im = 2*osparam->vals[indR].im/scR;
            }
        }

        /* Create the node and put into the first term in the array */
        Parr.vals[j] = tntMpsCreateTwoSiteOp(&Lterms,&Rterms,&two_site_params);
        /* Take the matrix exponential to form the two-site propagator */
        tntNodeExp(Parr.vals[j], hs, "DE", "UV");
    }
    
    /* Free the arrays and nodes */
    tntComplexArrayFree(&two_site_params);
    tntNodeArrayFree(&Lterms);
    tntNodeArrayFree(&Rterms);
    tntNodeFree(&eye);
    
    /* Turn the warning back on */
    tntSysQNClearWarnOn();

    /* Return the array of nodes */
    return Parr;
}
