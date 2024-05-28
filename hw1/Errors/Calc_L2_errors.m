function [E_L2] = Calc_L2_errors(femregion, uh, Data)
%% [E_L2, E_SEMI_H1]=C_error_L2_H1(femregion, uh, Data)
%==========================================================================
% Compute L2 and semi-H1 errors
%==========================================================================
%    called in ComputeErrors.m
%
%    INPUT:
%          Data        : (struct)  see Data.m
%          femregion   : (struct)  see CreateFemregion.m
%          uh          : (sparse(nfod,1))  solution vector
%
%    OUTPUT:
%          E_L2        : (real) L2 error 
%          E_SEMI_H1   : (real) H1 error - seminorm


nln = femregion.nln;
ne = femregion.ne;

% shape functions
[basis] = ShapeBasis;

% quadrature nodes and weights for integrals
[nodes_1D, w_1D] = Quadrature(3);

% evaluation of shape bases on quadrature nodes
[Phi, ~] = EvalShapeBasis(basis, nodes_1D);

E_L2_LOC = zeros(ne, 1);


% loop over the mesh elements
 for ie = 1 : ne
         
    iglo = femregion.connectivity(1:nln, ie);   
    
    [BJ, x] = GetJacobian(femregion.coord(iglo,:), nodes_1D);
    % BJ        = Jacobian of the elemental map 
    % pphys_1D  = vertex coordinates in the physical domain 
   
    local_uh = uh(iglo);
    
    % Solution and gradient at physical nodes x
	
    uex_loc = Data.uex(x,Data.T);
       
    % Approximated solution and gradient at quadrature nodes
    uh_loc = zeros(1,length(w_1D));
    
    for k = 1:length(w_1D)
        for s=1:nln
             uh_loc(k) = uh_loc(k) + Phi(1,k,s).*local_uh(s);
        end
    end
    
    % Evaluation integrals for errors
    for k = 1 : length(w_1D)
        Jdet = BJ;               % determinant
        dx = abs(Jdet).*w_1D(k); % wheight

        E_L2_LOC(ie) = E_L2_LOC(ie) + ((uh_loc(k)-uex_loc(k)).^2).*dx;
    end
end

% OUTPUT - h1 and l2 error
E_L2      = sqrt(sum(E_L2_LOC));
