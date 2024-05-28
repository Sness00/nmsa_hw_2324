function [A,M]=Matrix1D(Data,femregion)
%% [A,M] = Matrix1D(Data,Femregion)
%==========================================================================
% Assembly of the stiffness matrices A and rhs f
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          Data        : (struct)  see DataTest.m
%          femregion   : (struct)  see CreateFemregion.m
%
%    OUTPUT:
%          A           : (sparse(ndof,ndof) real) stiffnes matrix
%          M           : (sparse(ndof,ndof) real) mass matrix


fprintf('Assembling the matrices M and A... \n');


% connectivity infos
ndof         = femregion.ndof;         % degrees of freedom
nln          = femregion.nln;          % local degrees of freedom
ne           = femregion.ne;           % number of elements
connectivity = femregion.connectivity; % connectivity matrix


% shape functions
[basis] = ShapeBasis;

% quadrature nodes and weights for integrals
[nodes_1D, w_1D] = Quadrature(3);

% evaluation of shape bases on quadrature nodes
[Phi,GradPhi] = EvalShapeBasis(basis,nodes_1D);


% Assembly begin ...
A = sparse(ndof,ndof);  % Global Stiffness matrix
M = sparse(ndof,ndof);  % Global mass matrix

for ie = 1 : ne
     
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie);
  
    [BJ, ~] = GetJacobian(femregion.coord(iglo,:), nodes_1D);
    % BJ        = Jacobian of the elemental map 
    % pphys_1D  = vertex coordinates in the physical domain 
   
    %=============================================================%
    % STIFFNESS MATRIX
    %=============================================================%    
    for ii=1:length(nodes_1D)
        mu(ii) = Data.mu(femregion.coord(ie) + femregion.h*nodes_1D(ii));
    end
    [A_loc] = Stiffness(GradPhi, w_1D, nln, BJ, mu);
    A(iglo,iglo) = A(iglo,iglo) + A_loc; 
    
    %=============================================================%
    % MASS MATRIX
    %=============================================================%
    for ii=1:length(nodes_1D)
        ro(ii) = Data.ro(femregion.coord(ie) + femregion.h*nodes_1D(ii));
    end
    [M_loc] = Mass(Phi, w_1D, nln, BJ, ro);
    M(iglo,iglo) = M(iglo,iglo) + M_loc;
end