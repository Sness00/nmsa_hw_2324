function [dphiq,Grad] = EvalShapeBasis(basis, nodes_1D)
%% [dphiq,Grad] = EvalShapeBasis(basis,nodes_1D);

%==========================================================================
% Evaluation of basis functions and gradients on quadrature nodes 
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          basis    : (struct)  see C_shape_basis.m
%          node1D   : (array real) [nqn x 1]  coord of quadrature nodes   
%
%    OUTPUT:
%          dphi     : dphiq(:,:,i) vecotr containing the evaluation of the
%                     i-th basis function on the quadrature nodes 
%          Grad     : dphiq(:,1:2,i) vecotr containing the evaluation of the
%                     gradinet of dphi(:,:,i) on the quadrature nodes 


nln = length(basis);

%==========================================================================
% EVALUATION OF BASIS FUNCTIONS AND GRADIENTS
%==========================================================================

for s = 1:nln
    % coordinates
    csi = nodes_1D(:,1);
    % evaluation of basis functions
    dphiq(1,:,s) = eval(basis(s).fbases);
    % evaluation of gradients
    Grad(:,1,s) = eval(basis(s).Gbases);
end
