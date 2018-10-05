/*
Authors: Juan Jose Mendoza-Arenas and Stephen Clark
Date:    October 2018

(c) Universidad de Los Andes 2018
*/

/* Include the header for the TNT library */
#include "tntMps.h"
 /* Include the functions created by us */
#include "qtmscripts.h"

/*!

Main fuction to calculate time evolution of a fermionic lattice (in Jordan Wigner picture) coupled to a lead, which is itself coupled
to a Lindblad thermal reservoir. The dynamics of this quantum thermal machine (qtm) is simulated within the super-fermion approach

*/

int main(int argc, char **argv)
{
    char loadname[TNT_STRLEN]; /* The name of the initialisation file */
    char saveprefix[TNT_STRLEN]; /* Path to output file */
    char saveprefix_copy[TNT_STRLEN]; /* Copy of path to output file */
    char savefile[TNT_STRLEN]; /* Name of saving file */
    char loadstate[TNT_STRLEN]; /* Name of file from which initial state is read for intermediate simulations */
    char savestate[TNT_STRLEN]; /* Name of file in which evolved state is saved for intermediate simulations */

    tntNetwork wf, left_vacuum, prop; /* Networks representing the time-evolved MPS wave function, the "left" vacuum and propagator */
    int chi, L, N;                      /* Maximum internal dimension for the MPS, length of the system and of the lead */
    tntNodeArray nnl, nnr, os;          /* operators for building the hamiltonian */
    tntComplexArray nnparam, osparam;   /* parameters for building the hamiltonian */
    double err_step, err = 0.0;         /* The truncation error at each time step and total error */
    
    int numsteps, tbigstep;             /* Number of steps for run and between expectation values */
    int tstep = 0, bigstep = 0;         /* Counters for the above */
    tntComplex dt, norm, norm_inv;      /* Time step size, norm of evolved state and its inverse */

    int intermediate; /* Parameter determining if running the first or an intermediate part of the simulation */
    int current_chi, previous_chi, save_state; /* When running parts of simulation, the final current chi and the final previous chi, and flag to save state */
    int printstep, rescalestep; /* Number of time steps between printing advances and rescaling of the state */
    
    tntIntArray chi_list, t_list;  /* Length of each era, and the internal matrix dimension used for each one*/
    tntComplexArray dt_list; /* Time step */

    unsigned loop, loop_total=0, num_chi, era; /* For looping */

    tntMpsExOp ExOp; /* Defines all the operators for calculating expectation values */
    tntNodeArray swap_gate; /* Swap gate */

    /* Initialize the TNT library */
    tntInitialize();

    /* ------------------ Load basic information ------------------ */
    /* Process the command line options */
    tntProcessCLOptions(argc, argv, NULL, loadname, saveprefix, NULL, NULL);

    /* Load the simulation parameters */
    printf("Reading information from file %s.\n", loadname);
    tntIntParamsLoad(loadname, L, N, tbigstep, printstep, rescalestep);
    tntIntParamsLoad(loadname, intermediate, current_chi, previous_chi, numsteps, num_chi);
    tntIntParamsLoad(loadname, save_state);
    tntIntArraysLoad(loadname, chi_list, t_list);
    tntComplexArraysLoad(loadname, dt_list);
    tntNodeArraysNamedLoad(loadname,nnl,"nnlt",nnr,"nnrt",os,"ost",swap_gate,"swap_gate");
    tntComplexArraysNamedLoad(loadname,nnparam,"nnparamt",osparam,"osparamt");
    tntStringsLoad(loadname, savefile);

    /* Load "left" vacuum state and information for expectation values */
    tntNetworksLoad(loadname, left_vacuum);
    tntExOpsLoad(loadname, ExOp); /* Up to two sites */

    /* Print simulation parameters to the screen before starting the TEBD algorithm */
    printf("\n\nSystem size = %d, lead size = %d, final chi = %d.\n", L, N, current_chi);
    printf("numsteps = %d, final dt = %f, tbigstep = %d.\n", numsteps, dt_list.vals[0].re, tbigstep);

    /* Copy of path to saving file */
    sprintf(saveprefix_copy,"%s",saveprefix);

    /* Define complete name of output file */
    strcat(saveprefix, savefile);

    /* Save information of simulation */
    tntIntParamsSave(saveprefix, L, N, numsteps, tbigstep);
    tntIntParamsSave(saveprefix, current_chi, previous_chi);
    tntIntArraysSave(saveprefix, chi_list, t_list);
    tntComplexArraysSave(saveprefix, dt_list);

    /* Define initial state */
    if(intermediate == 0){ /* Make initial state equal to "left" vacuum (but normalized) if starting new simulation */
        tntNetworksLoad(loadname, wf);
    } else if(intermediate == 1){ /* Load wavefunction if intermediate calculation */
        /* Define name of files where state (S) and disorder in magnetic field (MF) are going to be loaded from */
        sprintf(loadstate,"%sSMF_chi%d",saveprefix_copy,previous_chi);
        strcat(loadstate, savefile);
        printf("Reading initial state from %s.\n", loadstate);
        /* Load function */
        tntNetworksLoad(loadstate, wf);
    }

    /* Find the norm of initial state, then rescale */
    norm = tntMpsMpsProduct(wf,left_vacuum);
    norm_inv.re = 1/norm.re; norm_inv.im = 0;
    tntNodeScaleComplex(tntNodeFindFirst(wf), norm_inv);
    printf("\nNorm of loaded initial state is %g.\n", norm.re);
    norm = tntMpsMpsProduct(wf,left_vacuum);
    printf("\nInitial state has been renormalised, with norm = %g.\n", norm.re);

    /* Calculate initial expectation values */
    tntMpsExpecOutputSuperFermion(wf, left_vacuum, &ExOp, 1, 1, saveprefix, bigstep);

    printf("\n\n");

    /* ------------------ Time evolution ------------------ */

    for(era = 0; era < num_chi; era++){
        
        /* Define dt and chi */
        dt = dt_list.vals[era];
        chi = chi_list.vals[era];

        printf("\nStarting simulation for super-fermion description of quantum thermal machine for chi = %d and dt = %g.\n\n", chi, dt.re);

        /* Create the propagator from the system parameters, using a Suzuki-Trotter 2nd-order staircase expansion.
           Up to site L (i.e. for the system) it performs the standard routine, but for sites > L (i.e. for the lead)
           it performs swap operations, moving site L around so when each local gate is applied, site L is nearest-
           neighbor to one of the lead modes. */
        prop = tntMpsCreatePropST2scPartialStar(L,N,dt,&nnl,&nnr,&nnparam,&os,&osparam,swap_gate);

        /* -------- Run the simulation for the current era -------- */
        for (loop = 1; loop <= t_list.vals[era]; loop++){
            
            /* Total time step */
            loop_total = loop_total + 1;

            /* Perform TEBD step */
            err += tntMpsPropST2scProduct(wf, prop, chi);

            /* Rescale state */
            if (loop%rescalestep == 0){
                norm = tntMpsMpsProduct(wf,left_vacuum);
                printf("Norm from main file = %g + i %g\n", norm.re, norm.im);
                norm_inv.re = norm.re/(pow(norm.re,2)+pow(norm.im,2));
                norm_inv.im = -1*norm.im/(pow(norm.re,2)+pow(norm.im,2));
                tntNodeScaleComplex(tntNodeFindFirst(wf), norm_inv);
                norm = tntMpsMpsProduct(wf,left_vacuum);
                printf("Norm from main file after rescaling = %g + i %g\n", norm.re, norm.im);
            }

            /* Calculate expectation values at set time intervals */
            if(loop%tbigstep == 0) {
                bigstep++;
                if(loop_total==numsteps){ /* Print and save in final step */
                    printf("---------- Final expectation values ----------\n");
                    tntMpsExpecOutputSuperFermion(wf, left_vacuum, &ExOp, 1, 1, saveprefix, bigstep);
                } else{ /* Only save for the rest of time steps */
                    tntMpsExpecOutputSuperFermion(wf, left_vacuum, &ExOp, 0, 1, saveprefix, bigstep);
                }

                /* Save information */
                tntDoubleParamsUpdate(saveprefix, bigstep, err);
            }

            /* Print advances */
            if (loop%printstep == 0) {
                printf("\nCompleted time step %d of %d\n", loop_total ,numsteps); /* Only print each certain number of time steps */
            }
        }

        /* Save NESS, if required, at the end of each era */
        if(save_state == 1){
            sprintf(savestate,"%sSMF_",saveprefix_copy);
            strcat(savestate, savefile);
            tntNetworksSave(savestate, wf); /* Save network */
        }

        /* Free propagator, as a new one will be defined for next era */
        tntNetworkFree(&prop);
    }

    /* ------------------ Free memory and finalise ------------------ */

    /*  Free all the dynamically allocated nodes and associated dynamically allocated arrays. */
    tntNetworkFree(&wf);
    tntNetworkFree(&left_vacuum);
    tntNodeArrayFree(&nnl);
    tntNodeArrayFree(&nnr);
    tntNodeArrayFree(&os);
    tntComplexArrayFree(&nnparam);
    tntComplexArrayFree(&osparam);
    tntNodeArrayFree(&swap_gate);
    tntMpsExOpFree(&ExOp);

    /* Finish with the TNT library */
    tntFinalize();

    return 0;
}
