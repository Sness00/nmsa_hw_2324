% FINITE VOLUME - BURGER's EQUATION
% CENTRAL FLUX 
clc
clear
close all

%%

N = 100;
x = linspace(0, 1, N+1);
dx = x(2:end) - x(1:end-1);
u_int = -1/(2*pi) * cos(2*pi*x);
u0 = (u_int(2:end)-u_int(1:end-1))./dx;
figure(1);
plot([x(1:end-1); x(2:end)], [u0; u0]);
u = u0;

dt = 0.01;
T = 1;
t_end = 0;

uplot = NaN(N, T/dt);
uplot(:, 1) = u';
j = 2;
fluxes = {'Constant', 'Upwind', 'Lax-Wendroff', 'Fromm'};
% fluxes = {'Constant', 'Upwind'};
[ind, ~] = listdlg('ListString', fluxes, 'SelectionMode', 'single', 'PromptString', 'Flux Reconstruction Scheme');

if ind == 1
    handle = @ddtFiniteVolumeConstant;
elseif ind == 2
    handle = @ddtFiniteVolumeUpwind;
elseif ind == 3
    handle = @ddtFiniteVolumeLW;
else
    handle = @ddtFiniteVolumeFromm;
end

for t = dt : dt : T
    [~, u] = ode45(handle, [0, dt], u);
    u = u(end,:);
    uplot(:,j) = u;
    plot([x(1:end-1); x(2:end)], [u;u], 'k');
    ylim([-1.5, 1.5]);
    xlabel('x'); ylabel('u(x,t)'); 
    t_end = t_end + dt;
    title(['Burgers equation at t = ', num2str(t_end)]);
    drawnow;
    j = j + 1;    
end

figure(2)
xplot = dx(1)/2:dx(1):1-dx(1)/2;
tplot = 0:dt:T;
surf(tplot, xplot, uplot,'EdgeColor','none')
xlabel('t'); ylabel('x', 'Rotation', 0); title('u(x,t)'); colorbar;
view(2)

function [dudt] = ddtFiniteVolumeConstant(~, u)
    N = length(u); % cell number
    x = linspace(0, 1, N+1); % grid
    dx = x(2:end) - x(1:end-1); % cell length
    % Constant Flux Reconstruction
    f = u.^2/2; % Flux    
    fL = f(1:end-1);
    fR = f(2:end);
    ul_gt_ur = u(1:end-1) - u(2:end) > 0;
    f_int = max(fL, fR).* ul_gt_ur + min(fL, fR) .* (~ul_gt_ur);        
    f_int = [0; f_int; 0];
    dudt = (f_int(1:end-1) - f_int(2:end))./dx';
end

function [dudt] = ddtFiniteVolumeUpwind(~, u)
    N = length(u);
    x = linspace(0, 1, N+1);
    dx = x(2:end) - x(1:end-1);
    u_plus = [u(1) - 0.5*u(1); u(2:end-1) - 0.5*(u(2:end-1) - u(1:end-2))];
    u_minus = u(2:end) + 0.5*(u(2:end) - u(1:end-1));  
    fL = u_plus.^2/2;
    fR = u_minus.^2/2;
    ul_gt_ur = u_plus - u_minus > 0;
    f_int = max(fL, fR).* ul_gt_ur + min(fL, fR) .* (~ul_gt_ur);        
    f_int = [0; f_int; 0];
    dudt = (f_int(1:end-1) - f_int(2:end))./(dx)';
end

function [dudt] = ddtFiniteVolumeLW(~, u)
    N = length(u); % cell number
    x = linspace(0, 1, N+1); % grid
    dx = x(2:end) - x(1:end-1); % cell length
    u_plus = u(1:end-2) - 0.5*(u(2:end-1) - u(1:end-2));
    u_minus = u(2:end) + 0.5*(u(3:end-1) - u(2:end));
    fL = u_plus.^2/2;
    fR = u_minus.^2/2;
    ul_gt_ur = u_plus - u_minus > 0;
    f_int = max(fL, fR).* ul_gt_ur + min(fL, fR) .* (~ul_gt_ur);        
    f_int = [0; f_int; 0];
    dudt = (f_int(1:end-1) - f_int(2:end))./dx';
end

function [dudt] = ddtFiniteVolumeFromm(~, u)
    N = length(u); % cell number
    x = linspace(0, 1, N+1); % grid
    dx = x(2:end) - x(1:end-1); % cell length
    u_plus = [u(1); u(2:end-1) - 0.25*(u(3:end) - u(1:end-2))];
    u_minus = [u(2:end-1) + 0.25*(u(3:end) - u(1:end-2)); u(end)];
    fL = u_plus.^2/2;
    fR = u_minus.^2/2;
    ul_gt_ur = u_plus - u_minus > 0;
    f_int = max(fL, fR).* ul_gt_ur + min(fL, fR) .* (~ul_gt_ur);        
    f_int = [0; f_int; 0];
    dudt = (f_int(1:end-1) - f_int(2:end))./dx';
end