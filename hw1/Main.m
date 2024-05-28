function [Errors, Solutions, Femregion, Data] = Main(Data, nEl)
%%
%    INPUT:
%          Data    : (struct) Data struct
%          nEl     : (int)    Number of mesh elements
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) finite element space
%
%          Data        : (struct)  Data struct
%

fprintf('============================================================\n')
fprintf(['Solving test ', Data.name, ' with ', num2str(nEl), ' elements \n']);

%==========================================================================
% MESH GENERATION
%==========================================================================

[Region] = CreateMesh(Data, nEl);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[Femregion] = CreateFemregion(Data, Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES and RIGHT-HAND SIDE
%==========================================================================

[A_no_bc, M_no_bc] = Matrix1D(Data, Femregion);

%==========================================================================
% BUILD FINITE ELEMENTS RHS a time 0
%==========================================================================

[b_nbc] = Rhs1D(Data, Femregion);

%==========================================================================
% BUILD INITIAL CONDITIONS
%==========================================================================

x = Femregion.coord;

% Time scheme parameters
gamma = 1/2;
beta = 1/4;

% I.C.
u0 = Data.u0(x);
v0 = Data.v0(x);

u(:, 1) = u0;
y(:, 1) = [u0; v0];

% Time loop
prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);
for t = Data.dt : Data.dt : Data.T
    M = M_no_bc;
    A = A_no_bc;
    b_bc = b_nbc(:, round(t/Data.dt)); % i = 1 => t = 0[s]
    b_bcpu = b_nbc(:, round(t/Data.dt)+1); % i = 2 => t = dt[s]
    u_g = 0;

    rhs = [Data.dt^2*beta.*b_bcpu + Data.dt^2*(1/2-beta).*b_bc;...
           Data.dt*gamma*b_bcpu + Data.dt*(1-gamma)*b_bc];

    M0 = [M-Data.dt^2*(1/2-beta)*A, Data.dt*M;...
          -Data.dt*(1-gamma)*A, M];

    M1 = [M+Data.dt^2*beta*A, zeros(size(M));...
          Data.dt*gamma*A, M];

    b = (M0*y(:, round(t/Data.dt))+rhs);

    if (strcmp(Data.boundary,'DD') || strcmp(Data.boundary(1),'D') ...
            || strcmp(Data.boundary(2),'D'))
        [M1_dbc, b_dbc, u_g] = BoundaryConditions(M1, b, Femregion, Data, t-Data.dt); 
        M1 = M1_dbc;
        b = b_dbc;
    end
    
    y_step = M1\b;

    y(:, round(t/Data.dt)+1) = y_step + u_g;
    
    u(:, round(t/Data.dt)+1) = y(1:length(u0), round(t/Data.dt)+1);
    
    if Data.snapshot 
        Snapshot(Femregion, u(:, round(t/Data.dt)+1), t);
    end

    prog = ( 100*(t/Data.T) );
    fprintf(1, '\b\b\b\b%3.0f%%', prog);
end

if Data.surf
    [xx,tt] = meshgrid(0:Data.dt:Data.T, Femregion.coord);
    figure
    surf(xx, tt, u, 'EdgeColor', 'None');
    xlabel('time'); ylabel('x'); zlabel('u(x,t)');
    title('u_h in space and time')
    view(2)
end

Solutions = [];
Errors = [];

if Data.visual_graph
    [Solutions] = PostProcessing(Data, Femregion, u(:, end));
end
    
if (Data.calc_errors)
        [Errors] = ComputeErrors(Data,Femregion,Solutions);
end


