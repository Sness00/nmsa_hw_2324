clc
clear
close all

%% Path
addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing

%% Data for Test
Data = DataTest('Q3a');

% Options
Data.surf = 1;
Data.snapshot = 1;
Data.visual_graph = 0;
Data.calc_errors = 0;

[err, sol, fem, D] = Main(Data, 200);

%% Data for Test
close all
Data = DataTest('test');

% Options
Data.surf = 0;
Data.snapshot = 0;
Data.visual_graph = 1;
Data.calc_errors = 1;

% Main routine
[err1, sol1, fem1, D1] = Main(Data, 16);
[err2, sol2, fem2, D2] = Main(Data, 32);
[err3, sol3, fem3, D3] = Main(Data, 64);
[err4, sol4, fem4, D4] = Main(Data, 128);
[err5, sol5, fem5, D5] = Main(Data, 256);
[err6, sol6, fem6, D6] = Main(Data, 512);

%% Plot Errors
close all
hVec   = [fem1.h, fem2.h, fem3.h, fem4.h, fem5.h, fem6.h]; 
eVecL2 = [err1.L2, err2.L2, err3.L2, err4.L2, err5.L2, err6.L2];
figure
loglog(hVec,hVec.^2,'-+b','Linewidth',2); hold on;
loglog(hVec,eVecL2,'-or','Linewidth',2);
legend(sprintf('h^%i',2),'||u-u_h||_{L^2}');
ylabel('L^2-error');
xlabel('h');