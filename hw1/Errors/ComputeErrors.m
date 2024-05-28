function [errors] = ComputeErrors(Data, femregion, solutions)
%% [errors] = C_compute_errors(femregion, solutions)
%==========================================================================
% Compute L2 error
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          Data        : (struct)  see DataTest.m
%          femregion   : (struct)  see CreateFemregion.m
%          solutions   : (struct)  see PostProcessing.m
%
%    OUTPUT:
%          errors      : (struct)

[E_L2] = Calc_L2_errors(femregion, solutions.uh, Data);

errors = struct('L2',   E_L2);