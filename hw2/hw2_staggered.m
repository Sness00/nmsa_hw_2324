clc
clear
close all

%% Data
L = 1;
c = 1;
gamma = c/L;
S1 = @(x) 1 + 0.*x;
S1_x = @(x) 0.*x;

S2 = @(x) (1 + 2.*x).^2;
S2_x = @(x) 4*(1 + 2.*x);

cross = input('Which section? ');
while(cross ~= 1 && cross ~= 2)
    clc
    cross = input('Which section? ');
end
if cross == 1
    S = S1;
    S_x = S1_x;
else
    S = S2;
    S_x = S2_x;
end
clc

p_ex = @(x, t) 1/gamma*3*pi*cos(pi*(x/2+1)).*cos(3*pi*t);
u_ex = @(x, t) S(x).*pi/2.*sin(pi*(x/2+1)).*sin(3*pi*t);

phi_0 = @(x) 0.*x;
phi_1 = @(x) 3*pi*cos(pi*(x/2+1));

f = @(x, t) gamma^2*S_x(x)*pi/2.*sin(pi*(x/2+1)).*sin(3*pi*t) - ...
            (9-gamma^2/4)*pi^2*S(x).*cos(pi*(x/2+1)).*sin(3*pi*t);

D = [0, L];
T = 2;

Nt = 6000;
Nx = 2000;

dt = T/Nt;
dx = (D(2) - D(1))/Nx;

t = 0:dt:T;
t_u = t;
t_p = t(1:end-1) + dt/2;
x_p = D(1):dx:D(2);
x_u = D(1)+dx/2:dx:D(2);

%% Solution

u = zeros(length(t_u), length(x_u));
p = zeros(length(t_p), length(x_p));

% I.C.
p_m12 = 1/gamma*phi_1(x_p);
u_0 = -S(x_u).*phi_0(x_u);

u(1, :) = u_0;

% B.C. (x = 0)
p(1, 1) = -gamma*dt/dx*u(1, 1)/S(x_p(1)) + dt/gamma*f(x_p(1), t_u(1))/S(x_p(1)) + p_m12(1);
for j = 2:Nx
    p(1, j) = -gamma*dt/dx*(u(1, j) - u(1, j-1))/S(x_p(j)) + dt/gamma*f(x_p(j), t_u(1))/S(x_p(j)) + p_m12(j);
end
%B.C. (x = 1)
p(1, end) = 0;

for k = 2:Nt
    for j = 1:Nx
        u(k, j) = -gamma*S(x_u(j))*dt/dx*(p(k-1, j+1)-p(k-1, j)) + u(k-1, j);
    end
    % B.C. (x = 0)
    p(k, 1) = -gamma*dt/dx*u(k, 1)/S(x_p(1)) + dt/gamma*f(x_p(1), t_u(k))/S(x_p(1)) + p(k-1, 1);
    for j = 2:Nx
        p(k, j) = -gamma*dt/dx*(u(k, j) - u(k, j-1))/S(x_p(j)) + dt/gamma*f(x_p(j), t_u(k))/S(x_p(j)) + p(k-1, j);
    end
    %B.C. (x = 1)
    p(k, end) = 0;
end

for j = 1:Nx
    u(Nt+1, j) = -gamma*S(x_u(j))*dt/dx*(p(Nt, j+1)-p(Nt, j)) + u(Nt, j);
end

[xxp, ttp] = meshgrid(x_p, t_p);
[xxu, ttu] = meshgrid(x_u, t_u);

%% Errors

L2_err_p = norm(p_ex(x_p, t_p(end)) - p(end, :), 2)*dx^0.5;
L1_err_p = norm(p_ex(x_p, t_p(end)) - p(end, :), 1)*dx;
Linf_err_p = norm(p_ex(x_p, t_p(end)) - p(end, :), Inf);

disp('Discretization errors for pressure: ')
disp([L2_err_p L1_err_p Linf_err_p])

L2_err_u = norm(u_ex(x_u, t_u(end)) - u(end, :), 2)*dx^0.5;
L1_err_u = norm(u_ex(x_u, t_u(end)) - u(end, :), 1)*dx;
Linf_err_u = norm(u_ex(x_u, t_u(end)) - u(end, :), Inf);

disp('Discretization errors for velocity: ')
disp([L2_err_u L1_err_u Linf_err_u])

%% Graphs

figure
subplot(1, 2, 1)
surf(xxp, ttp, p, 'EdgeAlpha', 0)
view(2)
title('Pressure p(x, t)')
xlabel('x')
ylabel('t', 'Rotation', 0)
colorbar
subplot(1, 2, 2)
surf(xxu, ttu, u, 'EdgeAlpha', 0)
view(2)
title('Velocity u(x, t)')
xlabel('x')
ylabel('t', 'Rotation', 0)
colorbar

figure
sgtitle('Pressure and Velocity at t = T')
subplot(2, 1, 1)
plot(x_p, p(end, :), 'LineWidth', 1.2)
hold on
plot(x_p, p_ex(x_p, t_p(end)), 'LineWidth', 1.2)
title('Pressure')
xlabel('x')
ylabel('p(x)')
legend('Computed', 'Exact', 'Location', 'southeast')
grid on
subplot(2, 1, 2)
plot(x_u, u(end, :), 'LineWidth', 1.2)
hold on
plot(x_u, u_ex(x_u, t_u(end)), 'LineWidth', 1.2)
title('Velocity')
xlabel('x')
ylabel('u(x)')
ylim([1.1*min(u(end, :)) 1.1*max(u(end, :))])
legend('Computed', 'Exact', 'Location', 'northeast')
grid on