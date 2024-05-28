function [f]=Rhs1D(Data,femregion)
%% [f] = Rhs1D(Data,Femregion)
%==========================================================================
% Assembly of rhs f for any time t
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          Data        : (struct)  see DataTest.m
%          femregion   : (struct)  see CreateFemregion.m
%
%    OUTPUT:
%          f           : (sparse(ndof,nT) real) load term


fprintf('Assembling the right hand side f... \n');


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
[Phi,~] = EvalShapeBasis(basis,nodes_1D);

time = 0 : Data.dt : Data.T;
f = sparse(ndof,length(time));     % Global Load vector

i = 1;
for t = time
    
    for ie = 1 : ne
        
        % Local to global map --> To be used in the assembly phase
        iglo = connectivity(1:nln,ie);
        
        [BJ, nodes_1D_phys] = GetJacobian(femregion.coord(iglo,:), nodes_1D);
        % BJ        = Jacobian of the elemental map
        % pphys_1D  = vertex coordinates in the physical domain
        
        %==============================================
        % FORCING TERM --RHS
        %==============================================
        
        % Local load vector
        [load] = Load(Data.force, Phi, BJ, w_1D, nodes_1D_phys, nln, t);
        
        % Assembly phase for the load vector
        f(iglo,i) = f(iglo,i) + load;
    end
    
    if strcmp(Data.boundary(1),'N')
        f(1,i) = f(1,i) - Data.gN1(t);
    end
    if strcmp(Data.boundary(2),'N')
        f(end,i) = f(end,i) + Data.gN2(t);
    end
    
    i = i + 1;
    
    
end




