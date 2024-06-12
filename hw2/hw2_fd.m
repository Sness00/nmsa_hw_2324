clc
clear
close all

%% Data and exact solution

L = 1;
c = 1;
gamma = c/L;

S1 = @(x) 1 + 0.*x;
S1_x = @(x) 0.*x;

S2 = @(x) (1 + 2*x).^2;
S2_x = @(x) 4*(1 + 2*x);

cross = input('Which surface? ');
while(cross ~= 1 && cross ~= 2)
    cross = input('Which surface? ');
end
if cross == 1
    S = S1;
    S_x = S1_x;
else
    S = S2;
    S_x = S2_x;
end

phi_ex = @(x, t) cos(pi(x/2 + 1)).*sin(3*pi*t);

phi_0 = @(x) 0.*x;
phi_1 = @(x) 3*pi*cos(pi*(x*0.5 + 1));

f = @(x, t) gamma^2*S_x(x)*pi/2.*sin(pi*(x*0.5+1)).*sin(3*pi*t) -...
    (9-gamma^2/4)*pi^2*S(x).*cos(pi*(x*0.5+1)).*sin(3*pi*t);

D = [0, L];
T = 2/3;

Nt = 4000;
Nx = 1000;

dt = T/Nt;
dx = (D(2) - D(1))/Nx;

t = 0:dt:T;
x = D(1):dx:D(2);

%% Solution

% Initial Conditions
phi = zeros(length(t), length(x));

phi(1, :) = phi_0(x);

phi(2, :) = phi_1(x)*2*dt;

for k = 2:Nt
    phi(k+1, 1) = 2*phi(k, 1) - phi(k-1, 1) + ...
                  dt^2/dx^2*gamma^2*(2*phi(k, 2)-2*phi(k, 1)) + ...
                  dt^2*f(x(1), t(k))/S(x(1));
    for j = 2:Nx
        phi(k+1, j) = 2*phi(k, j) - phi(k-1, j) + ...
                      dt^2*gamma^2/(2*dx)^2*(S(x(j+1))-S(x(j-1)))/S(x(j))*(phi(k, j+1)-phi(k, j-1)) + ...
                      dt^2/dx^2*gamma^2*(phi(k, j+1)-2*phi(k, j)+phi(k, j-1)) + ...
                      dt^2*f(x(j), t(k))/S(x(j));
    end
    phi(k+1, Nx+1) = phi(k-1, Nx+1);
end