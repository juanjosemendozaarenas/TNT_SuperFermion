%% Creating initialisation file for calculating exact NESS of fermionnic lattice coupled to a dissipative lead, using super-fermion approach

clear; clc;
path(path,'./tnt_matfiles/'); % Add path for common functions

%% Define simulation parameters

% System and lead parameters
L = 2; % Number of system sites.
N = 3; % Number of lead sites (eigenmodes).
M = L + N; % Total number of sites.
t_S = 1; % Hopping in the system.
t_L = 1.2; % Hopping in the lead.
t_c = 1; % System-lead coupling.
mu_s = 0; % The chemical potential of the system.
mu_l = 0; % The chemical potential of the lead.
gamma = 1; % The thermalising rate of the lead.
beta = 5; % Inverse temperature of the bath.

% Time evolution
T = 100; % Evolution time

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

% Create identities from 0 (number 1) to M-1 supersites (i.e. dimers), as
% well as Jordan-Wigner strings for supersites (sequence of sz x sz)
Identities{1} = 1;
JW{1} = 1;
for j = 2:M
    Identities{j} = kron(Identities{j-1},speye(d_dimer));
    JW{j} = kron(JW{j-1},kron(Sz,Sz));
end

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
        
% Hopping terms
hoplt = tntMatCreateOpArray({kron(Sm,Sz), kron(Sp,Sz), kron(I,Sm), kron(I,Sp)}); % Left-side nearest-neighbour operators
hoprt = tntMatCreateOpArray({kron(Sp,I), kron(Sm,I), kron(Sz,Sp), kron(Sz,Sm)}); % Right-side nearest-neighbour operators

hopping_list_system = t_S*ones(1,L-1);
hopping_list_lead = t_c*v;

hopparamt_system = [-1*hopping_list_system; % Hopping between system sites
            -1*conj(hopping_list_system); % Hermitian conjugate of the latter
            hopping_list_system; % Hopping between ancillas corresponding to first line
            conj(hopping_list_system)]; % Hermitian conjugate of the latter
        
hopparamt_lead = [-1*hopping_list_lead; % Hopping between last system site and lead modes
            -1*conj(hopping_list_lead); % Hermitian conjugate of the latter
            hopping_list_lead; % Hopping between ancillas corresponding to first line
            conj(hopping_list_lead)]; % Hermitian conjugate of the latter            
        
% Create Hamiltonian
H = sparse(zeros(d_dimer^M,d_dimer^M));

for k = 1:size(ost,2) % On-site terms  
    for j = 1:M
        H = H + osparamt(k,j)*kron(Identities{j},kron(ost(k).tensor.elems,Identities{M-j+1}));
    end
end

for k = 1:size(hoplt,2) % Hopping terms in system
    for j = 1:L-1
        H = H + hopparamt_system(k,j)*kron(Identities{j},kron(kron(hoplt(k).tensor.elems,hoprt(k).tensor.elems),Identities{M-j}));
    end
end

for k = 1:size(hoplt,2) % Hopping terms from last system site to lead modes. Between operators we musl put a JW string
    for j = 1:N
        H = H + hopparamt_lead(k,j)*kron(kron(Identities{L},hoplt(k).tensor.elems),kron(JW{j},kron(hoprt(k).tensor.elems,Identities{N-j+1})));
    end
end

%% Create "left" vacuum state and initial state
% The left vacuum is not normalized. The initial state will be its
% normalized version, created here

% Define vacuum for a single (physical or ancilla) site
vac = zeros(d,1); vac(1) = 1; % This corresponds to spin up

% Define state of each dimer, |up,up> + |dn,dn>
upup = kron(vac,vac); % Dimer state in |up,up> (physical and ancilla sites are up)
dndn = kron(Sm*vac,Sm*vac); % Dimer state in |dn,dn> (physical and ancilla sites are down)
superposition = upup+dndn;

% Define left vacuum
left_vacuum = 1;
for j = 1:M
    left_vacuum = kron(left_vacuum,superposition);
end

% Initialize starting state as normalized left vacuum
rho = left_vacuum;
norm = left_vacuum'*rho;
rho = rho/norm;
norm = left_vacuum'*rho;

disp(['Norm of initial state = ' num2str(norm)])

%% Perform time evolution and obtain expectation values

% Define operators for expectation values
density_local = kron(Sm*Sp,I);
for j = 1:M
    density_oper{j} = kron(kron(Identities{j},density_local),Identities{M-j+1}); 
end

% Initial density
for j = 1:M
    density0(1,j) = real(left_vacuum'*density_oper{j}*rho);
end

% Evolve initial time with non-Hermitian Hamiltonian
rho_t = expm(-1i*H*T)*rho;
norm = left_vacuum'*rho_t;
rho_t = rho_t/norm;

% for count = 1:50
%     rho_t = expm(-1i*H*T)*rho_t;
%     norm = left_vacuum'*rho_t
%     rho_t = rho_t/norm;
%     
%     for j = 1:M
%        density(count,j) = real(left_vacuum'*density_oper{j}*rho_t); 
%     end
%     
% end
% 
% density = [density0; density];
% 
% for j =1:M
%    plot(density(:,j),'-o')
%    hold on
% end

opts.disp = 0; opts.maxit = 1000; opts.isreal = 0; opts.issym = 0; % Eigensolver settings.
[rho_t,val] = eigs(H,1,'sm',opts); % Determine the right-eigenvector with zero eigenvalue.
norm = left_vacuum'*rho_t;
rho_t = rho_t/norm;

% Calculate populations
density_local = kron(Sm*Sp,I);
for j = 1:M
    density_oper{j} = kron(kron(Identities{j},density_local),Identities{M-j+1});
    density(j) = real(left_vacuum'*density_oper{j}*rho_t); 
end

% disp(['Densities = ' num2str(density(end,:))])
fprintf('%.10f ', density(end,:));