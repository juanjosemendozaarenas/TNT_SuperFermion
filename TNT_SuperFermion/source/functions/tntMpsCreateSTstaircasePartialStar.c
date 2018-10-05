/*
Authors: Juan Jose Mendoza Arenas, Sarah Al-Assam, Stephen Clark and Dieter Jaksch
 
(c) Universidad de los Andes 2018
*/

/* Include the headers for the TNT library */
#include "tnt.h"
#include "tntMps.h"

/* Include the header for all network functions included in djscripts */
#include "../qtmscripts.h"

/*!
 * Creates a network of two-site terms representing a site-wide propagator decomposed using a Suzuki-Trotter second-order staircase expansion,
 * for a system coupled to a non-interacting lead with dissipation, treated with the super-fermion approach.
 */
tntNetwork tntMpsCreatePropST2scPartialStar(unsigned L,                /*!< Length of the system. */
                                            unsigned N,                /*!< Number of lead modes. */
                                            tntComplex dtc,            /*!< Size of the time step. See the main description for information on how real and imaginary parts are applied */
                                            tntNodeArray *nnL,         /*!< Array of nearest-neighbour operators for left site. Send NULL if there are no nearest neighbour terms. */
                                            tntNodeArray *nnR,         /*!< Array of nearest-neighbour operators for right site. Send NULL if there are no nearest neighbour terms. */
                                            tntComplexArray *nnparam,  /*!< Array of parameters for nearest-neighbour operators. Send NULL if there are no nearest neighbour terms. */
                                            tntNodeArray *os,          /*!< Array of onsite operators. Send NULL if there are no on-site operators. */
                                            tntComplexArray *osparam,  /*!< Parameters for the on-site operators. Send NULL if there are no on-site operators. */
                                            tntNodeArray swap_gate)    /*!< Swap gate, two swap last site of the system across the lead */
{
    tntComplex h; /* Uniform scale factor for exponentials */
    tntNodeArray Proparr; /* Array of propagators */
    tntNetwork Prop; /* The propagator network to return */
    
    /* Create the scale factor from the time step */
    h.re = dtc.re/2.0;
    h.im = dtc.im/2.0;
    
    /* Generate an array of propagators required to build the network */
    Proparr = tntMpsCreatePropArrayPartialStar(L,N,h,nnL,nnR,nnparam,os,osparam);
    
    Prop = tntMpsPropArrayToST2scPartialStar(L,N,Proparr,swap_gate);
    
    /* Free the array of propagators (copies of all nodes taken for the network so we can free it) */
    tntNodeArrayFree(&Proparr);
    
    /* Return the network */
    return Prop;

}

/*! \ingroup tebd
 * 
 * Creates a network of two-site terms representing a site-wide propagator from an array of propagators decomposed using a Suzuki-Trotter second-order staircase expansion, for a system coupled to a non-interacting lead with dissipation, treated with the super-fermion approach.
 *
 * If j < L (number of system sites):
 *
 * Array entry `j` is assumed to contain the two-site propagator for sites \f$j,j+1\f$, and two copies of the propagator are placed in the network
 * i.e. one for a right to left sweep, and one for a left-to-right sweep. 
 * Such an array can be created using tntMpsCreatePropArray(), although it does not have to be (e.g. the nodes could be loaded from an initialisation file instead).
 *
 * If L <= j < L+N-1 (coupling between last system site and lead modes)
 *
 * A similar process is followed, but inserting the corresponding swap gates. The exception will be on the last gate, where no swap operation will take place so site L remains on the left of each two-site gate. To keep this order, from the left-to-right sweep the swap gate is applied after the corresponding two-site gate, while from the right-to-left the swap gate is applied before.
 *
 * \return A network representing the site-wide propagator
 */
tntNetwork tntMpsPropArrayToST2scPartialStar(unsigned L, /*!< Length of the system. */
                                             unsigned N, /*!< Number of lead modes. */
                                             tntNodeArray Proparr, /*!< Array of propagators. Uncahnged by the function - copies of all nodes are used. */
                                             tntNodeArray swap_gate) /*!< Swap gate, two swap last site of the system across the lead */
{
    tntNetwork Prop; /* The propagator network to return */
    tntNode Ptop; /* Last inserted node in the top of left-to-right staircase */
    tntNode Pbot; /* Last inserted node in the bottom of right-to-left staircase */
    tntNode local_gate; /* Local gate to be inserted in the propagator, which contains the swap gate when dealing with the lead sites (coupled to site L) */
    tntNode swap_copy; /* Copy of the swap gate */
    unsigned j; /* site counter */
    
    /* Create an empty network */
    Prop = tntNetworkCreate();
    
    /* Copy over first node */
    Ptop = tntNodeCopy(Proparr.vals[0]);
    
    /* Multiply by swap gate if there is only one site in the system; swap gate goes below */
    if(L==1) {
        swap_copy = tntNodeCopy(swap_gate.vals[0],0);
        tntNodeJoin(Ptop, "D", swap_copy, "U");
        tntNodeJoin(Ptop, "E", swap_copy, "V");
        Ptop = tntNodeContract(Ptop,swap_copy,NULL,NULL);
    }
    
    /* Add a leg to the first node in the network for joining to start of network */
    tntNodeAddLeg(Ptop,"L");
    
    /* Place a copy of the first node at the start of the network 
     * This represents the top of the staircase 
     * Left bottom leg points to the end of the network */
    tntNodeInsertAtStart(Ptop,"L","D",Prop);
    
    /* Copy over first node again */
    Pbot = tntNodeCopy(Proparr.vals[0]);
    
    /* Multiply by swap gate if there is only one site in the system; swap gate goes above */
    if(L==1) {
        swap_copy = tntNodeCopy(swap_gate.vals[0],0);
        tntNodeJoin(swap_copy, "D", Pbot, "U");
        tntNodeJoin(swap_copy, "E", Pbot, "V");
        Pbot = tntNodeContract(swap_copy,Pbot,NULL,NULL);
    }
    
    /* Add a leg to the first node in the network for joining to end of network */
    tntNodeAddLeg(Pbot,"L");
    
    /* Place a copy of the first node at the end of the network
     * This represents the end of the staircase or zip back and forth
     * The top left leg of this node connects to the bottom left leg of the topmost propagator that previously pointed to the end of the network */
    tntNodeInsertAtEnd(Pbot, "U", "L", Prop);
    
    /* Join the nodes on their right legs too */
    tntNodeJoin(Ptop, "E", Pbot, "V");
    
    /* Now move from left to right, inserting nodes below and right to the nodes forming the left to right staircase and above and right the nodes forming the right to left staircase. If L-1 <= j < N+L-2 (the first node is j = 0), insert as well the swap gate (for the last gate no swap operation is done). From left to right the swap is after the two-site gate, and from right to left it is before. */
    for (j = 1; j < Proparr.sz; j++) {
        
        /* Copy the node from the propagator array to be inserted */
        local_gate = tntNodeCopy(Proparr.vals[j]);
        
        /* Multiply by swap gate if required; the swap goes below the original two-site gate */
        if((j>=(L-1)) && (j<(N+L-2))) {
            swap_copy = tntNodeCopy(swap_gate.vals[0],0);
            tntNodeJoin(local_gate, "D", swap_copy, "U");
            tntNodeJoin(local_gate, "E", swap_copy, "V");
            local_gate = tntNodeContract(local_gate,swap_copy,NULL,NULL);
        }
        
        /* Insert new top node to the right of and between the current top node and bottom node */
        tntNodeInsert(local_gate, "U", "D", Ptop, "E", Pbot, "V");
        
        /* Reassign the current top node */
        Ptop = tntNodeFindConn(Ptop, "E");
        
        /* Copy again the node from the propagator array to be inserted */
        local_gate = tntNodeCopy(Proparr.vals[j]);        
        
        /* Multiply by swap gate if required; the swap goes above the original two-site gate */
        if((j>=(L-1)) && (j<(N+L-2))) {
            swap_copy = tntNodeCopy(swap_gate.vals[0],0);
            tntNodeJoin(swap_copy, "D", local_gate, "U");
            tntNodeJoin(swap_copy, "E", local_gate, "V");
            local_gate = tntNodeContract(swap_copy,local_gate,NULL,NULL);
        }
        
        /* Insert new bottom node directly under the new current top node and above and to the right of the current bottom node */
        tntNodeInsert(local_gate, "U", "D", Ptop, "D", Pbot, "V");
        
        /* Reassign the current bottom node */
        Pbot = tntNodeFindConn(Pbot, "V");
        
        /* Join the nodes on their right legs too */
        tntNodeJoin(Ptop, "E", Pbot, "V");
        
    }
    
    /* Return the network */
    return Prop;
}
