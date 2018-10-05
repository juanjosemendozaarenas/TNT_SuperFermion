%======================================================================
%> @ingroup matscripts
%> Creates matrices for N-species boson operators. These matrices can then be used to build Hamiltonian terms
%> If you wish to turn these matrices into nodes that can be loaded into the library use tntMatCreateOpArray().
%>
%> Returns the result in a cell having a number of elements equal to the number of species
%>
%> @param nmax Maximum number of bosons allowed on each site
%> @param N The number of species
%>
%> @retval bdag A cell, each entry containing the matrix representing \f$\hat{b}^{\dagger}\f$ for a given species. 
%> @retval b A cell, each entry containing the matrix representing \f$\hat{b}\f$ for a given species. 
%> @retval n A cell, each entry containing the matrix representing \f$\hat{n}\f$ for a given species. 
%======================================================================

function [bdag,b,n] = tntMatBosonOps(nmax,N)
%% tntMatBosonOps Creates N-species boson operators up to a maximum number per site of nmax

% Build sparse matrix version of basic spin operators :
bdags = diag(sqrt(1:nmax),-1);
bs = diag(sqrt(1:nmax),1);
ns = diag((0:nmax)',0);

% Construct spin operators for each spin in the full Hilbert space :
b = tntMatExpandBasis(bs,N);
bdag = tntMatExpandBasis(bdags,N);
n = tntMatExpandBasis(ns,N);


