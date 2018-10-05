/*
Authors: Juan Jose Mendoza Arenas and Stephen Clark
 
(c) Universidad de los Andes 2018
*/

/*! \file tntMpsExpecOutputSuperFermion.c
 *  \brief This contains the routine for calculating, displaying then saving the expectation value for open systems in the super-fermion technique
 */

/* Include the header for the TNT library */
#include "tnt.h"
#include "tntMps.h"

/* Include the header for all network functions included in djscripts */
#include "../qtmscripts.h"

/* Internal function to avoid repeated code */
void _tnt_MpsExpecOutputOneOp(tntComplexArray *Exval, int printOutput, int saveOutput, char *savename, char *label, unsigned counter);


/*!
 * \ingroup mps
 * 
 * Calculates various types of one-site and two-site expectation values for the super-fermion approach of open quantum systems.
 * These are obtained from the formula $\langle A(t)\rangle = \langle I|A|rho(t) \rangle$, with A the operator of interest,
 * $|rho(t) \rangle$ the evolved state (of systems + lead + ancilla) under the effective non-Hermitian Hamiltonian, and $|I\rangle$
 * the "left" vacuum. This function calculates $A|rho(t) \rangle$, and then calculates the overlap with $|I\rangle$.
 
 * The operators for which are passed in the ::tntMpsExOp structure.
 *
 * The function can output the expectation values to the screen as well as saving them to the named output file. 
 * Use the appropriate flags to specify whether results should be output to the screen or to an output file. 
 * The function will have no effect (and nothing will be calculated) if both the screen output flag is \c 0 and the printOutput flag is \c 0 - 
 * in this case a warning will be printed to screen but execution will continue.
 *
 * The function does not normalize the state, but the expectation values are normalized correctly. Note that the norm of the state $|rho(t) \rangle$ is
 * given by $\langle I|rho(t) \rangle$, i.e. the expectation value when the operator A is the identity
 *
 * \return No return value - the expectation values are output to the screen and/or to an output file.
 *
 * \see tntMpsProcessExpecOptions() for a function which converts command line arguments to the arrays of operators that are required by this function.
 */
void tntMpsExpecOutputSuperFermion(tntNetwork mps,      /*!< The network representing the MPS. Unchanged by the function. */
                                   tntNetwork left_vacuum, /*!< Left vacuum of the super-fermion method for open systems */
                                   tntMpsExOp *Op,      /*!< The operators for calculating expectation values */
                                   unsigned printOutput,/*!< Set to 1 to print expectation values to screen, 0 otherwise */
                                   unsigned saveOutput, /*!< Set to 1 to save values to output file, 0 otherwise */
                                   char *savename,      /*!< Path to the output file. If the extension is missing ".mat" is added. Only used if saveOutput is 1. */
                                   unsigned counter)    /*!< can be used if saving multiple batches of expectation values e.g. for different timesteps.
                                                 It appends the values to a pre-exisiting array in the position given by counter.
                                                 Pass zero if a counter is not required. */
{

    tntNetwork mpsc;        /* Copy to the input network mps */
    tntComplexArray Exval;  /* array for holding expectation values */
    tntComplex normsq;      /* Its real part is the norm squared of the wave function */
    double norm;            /* Real component of normsq */
    unsigned loop;          /* Used to loop through the different expectation value types */
    tntNodeArray op;        /* Node array for holding operators for expectation value currently being calculated */
    tntIntArray sitenum;    /* Array for holding site number for expectation value currently being calculated */
    int L;                  /* Length of the MPS */
    
    if (!printOutput && !saveOutput) {
        tntWarningPrint("No expectation values are being calculated as both output flags are set to off"); /* NO_COVERAGE */
        return; /* NO_COVERAGE */
    }  /* NO_COVERAGE */
    
    tntSysQNClearWarnOff();
    
    /* Find the length of the MPS */
    L = (int) tntMpsLength(mps);
    
    /* Calculate the norm squared */
    normsq = tntMpsMpsProduct(mps,left_vacuum);
    norm = normsq.re;
    
    if (printOutput) {
        tntPrintf("\n~~~~~~~             Norm squared from ExpVal function is % #6.6g             ~~~~~~~~\n\n",norm);
    }
    if (saveOutput) {
        if (counter) {
            tntDoubleParamsUpdate(savename, counter, norm);
        } else {
            tntDoubleParamsSave(savename, norm);
        }
    }
    
    /* ----- Onsite expectation values ----- */
    
    if (Op->ExOpOs.sz > 0) {
    
        /* Allocate arrays of size one for onsite operators */
        op = tntNodeArrayAlloc(1);
        sitenum = tntIntArrayAlloc(1);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~          Onsite expectation values          ~~~~~~~~\n");
        }
    }
    
    /* Find the single site expectation values of the calculated ground state */
    for (loop = 0; loop < Op->ExOpOs.sz; loop++) {
        
        /* Allocate memory for storing single-site expectation values */
        Exval = tntComplexArrayAlloc(L);
        
        /* Copy over the required operator */
        op.vals[0] = tntNodeCopy(Op->ExOpOs.vals[loop],0);
        
        for (sitenum.vals[0] = 0; sitenum.vals[0] < L; (sitenum.vals[0])++) {
            /* Find the expectation value on each site using a product MPO. */
            mpsc = tntNetworkCopy(mps);
            tntMpsPmpoProduct(mpsc, &op, &sitenum);
            Exval.vals[sitenum.vals[0]] = tntMpsMpsProduct(mpsc, left_vacuum);
            tntNetworkFree(&mpsc);
        }
        printf("\n");
        printf("\n");
        
        /* output and free values */
        _tnt_MpsExpecOutputOneOp(&Exval, printOutput, saveOutput, savename, Op->LbOpOs.vals[loop], counter);
        
        /* Free the copy of the operator */
        tntNodeFree(&(op.vals[0]));
    }
    
    if (Op->ExOpOs.sz > 0) {
        
        /* free the arrays */
        tntIntArrayFree(&sitenum);
        tntNodeArrayFree(&op);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        }
    }
    
    /* ---------- Nearest-neighbour expectation values ----------- */
    
    if (Op->ExOp2nn.sz > 0) {
        
        /* Allocate arrays of size two for two-site operator site numbers */
        sitenum = tntIntArrayAlloc(2);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~    Nearest-neighbour expectation values     ~~~~~~~~\n");
        }
    }
    
    for (loop = 0; loop < Op->ExOp2nn.sz; loop+=2) {
        
        /* Allocate memory for storing two-site expectation values */
        Exval = tntComplexArrayAlloc(L-1);
        
        /* Allocate arrays of size two for two-site operators */
        op = tntNodeArrayAlloc(2);
        
        /* Copy over the required operator */
        op.vals[0] = tntNodeCopy(Op->ExOp2nn.vals[loop],0);
        op.vals[1] = tntNodeCopy(Op->ExOp2nn.vals[loop+1],0);
        
        for (sitenum.vals[0] = 0; sitenum.vals[0] < L-1; (sitenum.vals[0])++) {
            /* assign value for second site */
            sitenum.vals[1] = sitenum.vals[0]+1;
            /* Find the expectation value on each site using a product MPO. */
            mpsc = tntNetworkCopy(mps);
            tntMpsPmpoProduct(mpsc, &op, &sitenum);
            Exval.vals[sitenum.vals[0]] = tntMpsMpsProduct(mpsc, left_vacuum);
            Exval.vals[sitenum.vals[0]].re /= norm;
            Exval.vals[sitenum.vals[0]].im /= norm;
            tntNetworkFree(&mpsc);
        }
        
        /* Free the operators */
        tntNodeArrayFree(&op);
        
        /* output and free values */
        _tnt_MpsExpecOutputOneOp(&Exval, printOutput, saveOutput, savename, Op->LbOp2nn.vals[loop/2], counter);
    }
    
    if (Op->ExOp2nn.sz > 0) {
        
        /* free the arrays */
        tntIntArrayFree(&sitenum);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        }
    }
    
    /* ---------- Centre-site expectation values ----------- */
    
    if (Op->ExOp2cs.sz > 0) {
        
        /* Allocate arrays of size two for two-site operator site numbers */
        sitenum = tntIntArrayAlloc(2);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~       Centre-site expectation values        ~~~~~~~~\n");
        }
    }
    
    for (loop = 0; loop < Op->ExOp2cs.sz; loop+=2) {
        
        /* Allocate memory for storing two-site expectation values */
        Exval = tntComplexArrayAlloc(L);

        /* Allocate first site that expectation value will be found on */
        sitenum.vals[0] = L/2;
        
        /* Allocate arrays of size two for two-site operators */
        op = tntNodeArrayAlloc(2);
        
        /* Copy over the required operator */
        op.vals[0] = tntNodeCopy(Op->ExOp2cs.vals[loop]);
        op.vals[1] = tntNodeCopy(Op->ExOp2cs.vals[loop+1]);
        
        for (sitenum.vals[1] = 0; sitenum.vals[1] < L; (sitenum.vals[1])++) {
            /* Find the expectation value on each site using a product MPO. */
            mpsc = tntNetworkCopy(mps);
            tntMpsPmpoProduct(mpsc, &op, &sitenum);
            Exval.vals[sitenum.vals[1]] = tntMpsMpsProduct(mpsc, left_vacuum);
            Exval.vals[sitenum.vals[1]].re /= norm;
            Exval.vals[sitenum.vals[1]].im /= norm;
            tntNetworkFree(&mpsc);
        }
        
        /* Free the operators */
        tntNodeArrayFree(&op);
        
        /* output and free values */
        _tnt_MpsExpecOutputOneOp(&Exval, printOutput, saveOutput, savename, Op->LbOp2cs.vals[loop/2], counter);
    }
    
    if (Op->ExOp2cs.sz > 0) {
        
        /* free the arrays */
        tntIntArrayFree(&sitenum);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        }
    }
    
    
    /* ---------- All pairs of two-site expectation values ----------- */
    
    /* Between the sites of interest there are JW strings, passed in the Maltab intialization file as an entry in ExOp.
       So the expectation values are of the form  < A_{i} O_{i+1} O_{i+2} ... O_{j-1} B_{j} > for operators of interest
       A and B, with O the JW string operators */
    
    int count, ind, x; /* Counter, total matrix index, and total length of string (including sites of interest) */
    int site_l, site_r; /* Sites of left and right operators */
    
    if (Op->ExOp2ap.sz > 0) {
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~~~ All pairs of two-site expectation values ~~~~~~~~~\n");
        }
    }
    
    for (loop = 0; loop < Op->ExOp2ap.sz; loop+=3) {
        
        /* Allocate memory for storing two-site expectation values */
        Exval = tntComplexArrayAlloc(L,L);
        
        /* Find the expectation value on each different pair of sites using a product MPO, excluding operators on the same site. */
        for (site_l = 0; site_l < L-1; site_l++) {
            for (site_r = site_l+1; site_r < L; site_r++) {
                
                /* Total lenght of string */
                x = site_r - site_l + 1;
                
                /* Allocate array for site numbers of operators */
                sitenum = tntIntArrayAlloc(x);
                
                /* Allocate array of operators */
                op = tntNodeArrayAlloc(x);
                
                /* Copy over the required operators */
                op.vals[0] = tntNodeCopy(Op->ExOp2ap.vals[loop]);
                sitenum.vals[0] = site_l;
                for (count = 1; count < x-1; count ++){
                    op.vals[count] = tntNodeCopy(Op->ExOp2ap.vals[loop+1]);
                    sitenum.vals[count] = site_l + count;
                }
                op.vals[x-1] = tntNodeCopy(Op->ExOp2ap.vals[loop+2]);
                sitenum.vals[x-1] = site_r;
                
                ind = site_l + L*site_r;
                mpsc = tntNetworkCopy(mps);
                tntMpsPmpoProduct(mpsc, &op, &sitenum);
                Exval.vals[ind] = tntMpsMpsProduct(mpsc, left_vacuum);
                Exval.vals[ind].re /= norm;
                Exval.vals[ind].im /= norm;
                tntNetworkFree(&mpsc);
                
                /* Free the operators and list of sites */
                tntNodeArrayFree(&op);
                tntIntArrayFree(&sitenum);
            }
        }
        
        /* output and free values */
        _tnt_MpsExpecOutputOneOp(&Exval, printOutput, saveOutput, savename, Op->LbOp2ap.vals[loop/2], counter);
    }
    
    if (Op->ExOp2ap.sz > 0) {
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        }
    }

    tntSysQNClearWarnOn();
}

/* prints, saves and frees the complex array, converting to double if necessary, and taking account of flags and counters */
void _tnt_MpsExpecOutputOneOp(tntComplexArray *Exval, int printOutput, int saveOutput, char *savename, char *label, unsigned counter) {
    
    tntDoubleArray Exvalr; /* array for holding real part of expectation values */
    tntComplexArray Exvalc; /* array for holding real part of expectation values */
    
    /* Now take different action depending on whether the expectation values are real or complex */
    if (tntComplexArrayIsReal(Exval)) {
        /* change complex array to real array */
        Exvalr = tntComplexArrayToReal(Exval);
        /* Print real output */
        if (printOutput) {
            tntNamedDoubleArrayPrint(label,&Exvalr);
        }
        /* Save this expectation value */
        if (saveOutput) {
            if (counter) tntDoubleArraysNamedUpdate(savename, counter,Exvalr, label);
            else         tntDoubleArraysNamedSave(savename, Exvalr, label);

        }
        
        /* free real array */
        tntDoubleArrayFree(&Exvalr);
    } else {
        Exvalc = *Exval;
        /* Print complex output */
        if (printOutput) {
            tntNamedComplexArrayPrint(label,Exval);
        }
        /* Save complex output */
        if (saveOutput) {
            if (counter) tntComplexArraysNamedUpdate(savename, counter, Exvalc, label);
            else         tntComplexArraysNamedSave(savename, Exvalc, label);
        }
        /* Free complex array */
        tntComplexArrayFree(Exval);
    }
}

