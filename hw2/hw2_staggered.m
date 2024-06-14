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

cross = input('Which surface? ');
while(cross ~= 1 && cross ~= 2)
    clc
    cross = input('Which surface? ');
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

Nt = 2000;
Nx = 1000;

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
    u(Nt+1, j) = -gamma*S(x_p(j))*dt/dx*(p(Nt, j+1)-p(Nt, j)) + u(Nt, j);
end

u_exc = zeros(length(t_u), length(x_u));
p_exc = zeros(length(t_p), length(x_p));
for k = 1:Nt+1
    for j = 1:Nx
        u_exc(k, j) = u_ex(x_u(j), t_u(k));
    end
end
for k = 1:Nt
    for j = 1:Nx+1
        p_exc(k, j) = p_ex(x_p(j), t_p(k));
    end
end
[xxp, ttp] = meshgrid(x_p, t_p);
[xxu, ttu] = meshgrid(x_u, t_u);

%% Errors

L2_err_p = norm(p_exc(end, :) - p(end, :), 2)*dx^0.5;
L1_err_p = norm(p_exc(end, :) - p(end, :), 1)*dx;
Linf_err_p = norm(p_exc(end, :) - p(end, :), Inf);

disp('Discretization errors for pressure: ')
disp([L2_err_p L1_err_p Linf_err_p])

L2_err_u = norm(u_exc(end, :) - u(end, :), 2)*dx^0.5;
L1_err_u = norm(u_exc(end, :) - u(end, :), 1)*dx;
Linf_err_u = norm(u_exc(end, :) - u(end, :), Inf);

disp('Discretization errors for velocity: ')
disp([L2_err_u L1_err_u Linf_err_u])

%% Graphs

figure
sgtitle('Velocity')
subplot(1, 2, 1)
surf(xxu, ttu, u, 'EdgeAlpha', 0)
view(2)
title('Computed')
subplot(1, 2, 2)
surf(xxu, ttu, u_exc, 'EdgeAlpha', 0)
view(2)
title('Exact')

figure
sgtitle('Pressure')
subplot(1, 2, 1)
surf(xxp, ttp, p, 'EdgeAlpha', 0)
view(2)
title('Computed')
subplot(1, 2, 2)
surf(xxp, ttp, p_exc, 'EdgeAlpha', 0)
view(2)
title('Exact')

figure
sgtitle('Pressure and Velocity at t = T')
subplot(2, 1, 1)
plot(x_p, p(end, :))
hold on
plot(x_p, p_exc(end, :))
title('Pressure')
legend('Computed', 'Exact')
subplot(2, 1, 2)
plot(x_u, u(end, :))
hold on
plot(x_u, u_exc(end, :))
title('Velocity')
legend('Computed', 'Exact')