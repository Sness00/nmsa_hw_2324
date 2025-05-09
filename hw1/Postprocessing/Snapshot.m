function [uh] = Snapshot(femregion, uh, t)
%% Snapshot(femregion, uh, Dati,t)
%==========================================================================
% PLOT THE EXACT SOLUTION ON THE DOFS
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          femregion   : (struct)  see C_create_femregion.m
%          uh          : (sparse(ndof,1) real) solution vector
%          Dati        : (struct) see C_Dati.m
%


x1 = femregion.domain(1,1);
x2 = femregion.domain(1,2);


M =  2; % max(uh)
m = -2; % min(uh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PLOT OF SOLUTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(100)
plot(femregion.coord(:,1),full(uh));
title(['u_h(x, t) at time t =  ', num2str(t)]); xlabel('x'); ylabel('y', 'Rotation', 0);
axis([x1,x2,m,M]); 
grid on
end
