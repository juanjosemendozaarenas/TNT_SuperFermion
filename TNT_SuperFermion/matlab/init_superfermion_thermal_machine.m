%% Creating initialisation file for calculating NESS of fermionnic lattice coupled to a dissipative lead, using super-fermion approach

clear; clc;
path(path,'./tnt_matfiles/'); % Add path for common functions

%% Define simulation parameters

count_file = 3;

% System and lead parameters
L = 80; % Number of system sites.
N = 80; % Number of lead sites (eigenmodes).
M = L + N; % Total number of sites.
t_S = 1; % Hopping in the system.
t_L = 1.2; % Hopping in the lead.
t_c = 1; % System-lead coupling.
mu_s = 0; % The chemical potential of the system.
mu_l = 0; % The chemical potential of the lead.
gamma = 1; % The thermalising rate of the lead.
beta = 10; % Inverse temperature of the bath.

% Time evolution
chi_list = [30 60 90]; % Internal matrix dimension for each era. Expectation values are saved to output file at the end of each era.
num_chi = length(chi_list); % Number of different chi values
t_list = [10 10 8]*1000; % Time steps for each era 
numsteps = sum(t_list); % Total temporal extent of the simulation 
%dt_list = [0.05 0.05 0.04 0.025]; % Size of a timestep
dt_list = [0.03 0.02 0.01];
tbigstep = 2000; % Time between each evaluation of expectation values
rescalestep = 2000; % Number of time steps between rescaling
printstep = 2000; % Number of time steps between printing where we are in simulation
save_state = 2000; % 1 to save state in specidied file, 0 otherwise

% Parameters used when reading the state from a previous simulation 
intermediate = 0; % 0 if this file gives info for the intial part of the simulation, 1 if it gives info for an intermediate part of the simulation
current_chi = chi_list(end); % Final chi of the current simulation
previous_chi = 120; % Final chi of the previous simulation

%% Diagonalize noninteracting lead

% Construct lead single-particle Hamiltonian:
H_lead = -t_L*diag(ones(1,N-1),-1) - t_L*diag(ones(1,N-1),1) - mu_l*eye(N);
H_lead(1,N) = -t_L; % Periodic boundary conditions on the lead.
H_lead(N,1) = -t_L;

% Numerically diagonalise the lead Hamiltonian:
[U_lead,en_lead] = eig(H_lead);
en_lead = diag(en_lead).'; % Eigenenergies in ascending order.
v = U_lead(1,:); % Eigenstate amplitudes on the first lead site give couplings.
f = 1./(1 + exp(beta*en_lead)); % Fermi occupation factors for each bath site.

%% Define important operators

% Create Pauli matrices (single-site operators)
[Sp,Sm,Sx,Sy,Sz,I] = spin_hamiltonian(2);
d = size(Sz,1); % Dimension of a single physical, lead or ancilla site
d_dimer = d*d; % Dimension of a dimer (physical or lead site + corresponding ancilla)
    
FSWAP = fermionic_swap_gate; % Build the fermionic swap gate
%FSWAP = eye(16);
FSWAP = reshape(FSWAP,[d_dimer,d_dimer,d_dimer,d_dimer]);
swap_gate = tntMatCreate2SitesOpArray({FSWAP});

%% Defining the physical basis and symmetry information
% The Hamiltonian does not conserve particle number, so no quantum number
% information is given to the simulation

qnums = [];

tntSystem = tntMatCreateBasisOp(Sz,qnums);
tntSystem.sysnum = 1; % Type of system (bosonic of spin). Not used by the code.
tntSystem.symm_type = 0; % No symmetry is used
tntSystem.symm_num_qn = 0; % Number of conserved quantities. Equal to ns e.g. if the number of each species is conserved

%% Global parameters used in linear algebra routines
% These parameters are used while taking SVDs. They are (in order)
% * The tolerance for zeroing matrix values in automatic blocking
% * The absolute value below which all singular values will be discarded
% * The ratio between smallest singular value that will be kept and the largest singular value
tntSystem.zero_tol = 1e-10;
tntSystem.abs_trunc_tol = -1;
tntSystem.rel_trunc_tol = 1e-10;
tntSystem.trunc_err_tol = -1;

% Define the function that will be used for calculating the truncation
% error. Choose from 'sumsquares', '1norm', '1normscaled', '2norm',
% '2normscaled'.
tntSystem.trunc_err_func = 'sumsquares';

% Define the type of SVD to use.
tntSystem.svdtype = 1;

% Define the maximum number of iterations for the iterative eigenvalue
% solver. You may want to change this if you get non-convergance errors.
tntSystem.maxeigiter = 300;

% Determine whether reshape re-use is on or off. It is best to have it on
% (|1|) if you have many tensors of a similar size (i.e. $\chi$ is uniform
% throughout the lattice) and off (|0|) if this is not true
tntSystem.reshape_reuse = 1;

%% Define the non-Hermitian Hamiltonian of the super-fermion approach
% The whole lattice will be described by dimers (or supersites), consisting 
% of a physical site (from the system or the lead) and the corresponding  
% ancilla site. In the dimer, operators AB come from the direct product of  
% operators from the physical site (A) and its ancilla (B), with order A x B. 

% On-site terms
ost = tntMatCreateOpArray({kron(Sm*Sp,I)-kron(I,Sm*Sp), kron(Sm*Sp,I)+kron(I,Sm*Sp), kron(Sp,Sp), kron(Sm,Sm), kron(I,I)}); % On-site operators
osparamt = [mu_s*ones(1,L-1) mu_s/(N+1) en_lead; % On-site energies. That of site L is divided by N+1, as it will be repeated N+1 times (1 with site L-1, and N with lead modes)
            zeros(1,L) -1i*0.5*gamma*(1-2*f); % Dissipative particle-conserving term, only acting on the lead
            zeros(1,L) 1i*gamma*(1-f); % Dissipative gain process, only acting on the lead
            zeros(1,L) 1i*gamma*f; % Dissipative decay process, only acting on the lead
            zeros(1,L) -1i*gamma*f]; % Constant term       
        
% Nearest-neighbor terms
nnlt = tntMatCreateOpArray({kron(Sm,Sz), kron(Sp,Sz), kron(I,Sm), kron(I,Sp)}); % Left-side nearest-neighbour operators
nnrt = tntMatCreateOpArray({kron(Sp,I), kron(Sm,I), kron(Sz,Sp), kron(Sz,Sm)}); % Right-side nearest-neighbour operators

hopping_list = [t_S*ones(1,L-1) t_c*v];
nnparamt = [-1*hopping_list; % Hopping between system sites, or between last system site and lead modes
            -1*conj(hopping_list); % Hermitian conjugate of the latter
            hopping_list; % Hopping between ancillas corresponding to first line
            conj(hopping_list)]; % Hermitian conjugate of the latter
        
%% Expectation values to take
% For the current, the calculation could be modified so that for system
% sizes it is done as usual, but for the lead (where the current is between 
% each lead mode and site L of the system) a string is calculated!!!!!!!!!!
ExOp.os_operators = tntMatCreateOpArray({kron(Sm*Sp,I)}); % On-site expectation values
ExOp.os_labels = {'density'};
% ExOp.nn_operators = tntMatCreateOpArray({kron(Sm,Sz),kron(Sp,I),kron(Sp,Sz),kron(Sm,I)}); % Two-site nearest neighbours
% ExOp.nn_labels = {'SmSp','SpSm'};
ExOp.nn_operators = tntMatCreateOpArray({}); % Two-site nearest neighbours
ExOp.nn_labels = {};
ExOp.cs_operators = tntMatCreateOpArray({}); % Centered-site
ExOp.cs_labels = {};
ExOp.ap_operators = tntMatCreateOpArray({kron(Sm,Sz),kron(Sz,Sz),kron(Sp,I)}); % All-pairs; first operator for left site, central operator for JW string, third operator for right site 
ExOp.ap_labels = {'G'};

%% Create "left" vacuum state and initial state
% The left vacuum is not normalized. The initial state will be its
% normalized version, created here (THINK OF BETTER CHOICES FOR THE INITIAL STATE LATER ON).

% Create space for state
left_vacuum = cell(1,M);

% Define vacuum for a single (physical or ancilla) site
vac = zeros(d,1); vac(1) = 1; % This corresponds to spin up

% Define state of each dimer, |up,up> + |dn,dn>
upup = kron(vac,vac); % Dimer state in |up,up> (physical and ancilla sites are up)
dndn = kron(Sm*vac,Sm*vac); % Dimer state in |dn,dn> (physical and ancilla sites are down)

% Store state in each dimer
for site = 1:M
     left_vacuum{site} = upup+dndn;
end

% Make the initial state equal to the unnormalized left vacuum
wf = left_vacuum;

% Put in TNT language, with wf normalized but left_vacuum not
left_vacuum = tntMatCreateProdMpsNotNormalized(left_vacuum,qnums);
wf = tntMatCreateProdMps(wf,qnums);
        
%% Save information

% Save output of simulation
savefile = ['Thermal_L' num2str(L) '_N' num2str(N) '_tS' num2str(t_S) '_tL' num2str(t_L) '_tc' num2str(t_c) '_gamma' num2str(gamma) '_beta' num2str(beta) '.mat'];
% Save the parameters
fname = ['../initfiles/init_' num2str(count_file) '_supfermi.mat'];
save(fname);