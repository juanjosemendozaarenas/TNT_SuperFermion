
%======================================================================
%> @ingroup matscripts
%> Creates matrices for N-species spin operators. These matrices can then be used to build Hamiltonian terms
%> If you wish to turn these matrices into nodes that can be loaded into the library use tntMatCreateOpArray().
%>
%> Returns the result in a cell having a number of elements equal to the number of species
%>
%> @param s Required spin
%> @param N The number of species
%>
%> @retval sx A cell, each entry containing the matrix representing \f$\hat{S}^x\f$ for a given species. 
%> @retval sy A cell, each entry containing the matrix representing \f$\hat{S}^y\f$ for a given species. 
%> @retval sz A cell, each entry containing the matrix representing \f$\hat{S}^z\f$ for a given species. 
%> @retval sp A cell, each entry containing the matrix representing \f$\hat{S}^p\f$ for a given species. 
%> @retval sm A cell, each entry containing the matrix representing \f$\hat{S}^m\f$ for a given species. 
%======================================================================

function [sx,sy,sz,sp,sm] = tntMatSpinOps(s,N)
%% tntMatBosonOps Creates N-species spin operators with a spin given by TwoS

% z component of spin
m = s:-1:-s;

% Build sparse matrix version of basic spin operators :
sz = diag(m);
sp = diag(sqrt((s-m(2:end)).*(s+m(2:end)+1)),1);
sm = diag(sqrt((s+m(1:(end-1))).*(s-m(1:(end-1))+1)),-1);

sx = 0.5*(sp+sm);
sy = -0.5*1i*(sp-sm);

% Construct spin operators for each spin in the full Hilbert space :
sx = tntMatExpandBasis(sx,N);
sy = tntMatExpandBasis(sy,N);
sz = tntMatExpandBasis(sz,N);
sp = tntMatExpandBasis(sp,N);
sm = tntMatExpandBasis(sm,N);
