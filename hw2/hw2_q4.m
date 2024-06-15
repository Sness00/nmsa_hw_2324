clc
clear
close all

%% Data

L = 1;
c = 2;
gamma = c/L;

phi_0 = @(x) 0.*x;
phi_1 = @(x) 0.*x;
f = @(x, t) 0.*x.*t;

u_b = @(t) 0.5*(sin(pi*(t-0.1)/0.05) + abs(sin(pi*(t-0.1)/0.05)));

D = [0, L];
T = 1;

Nt = 4000;
Nx = 1000;

dt = T/Nt;
dx = (D(2) - D(1))/Nx;

t = 0:dt:T;
t_u = t;
t_p = t(1:end-1) + dt/2;
x_p = D(1):dx:D(2);
x_u = D(1)+dx/2:dx:D(2);

sec = input('Which section? ');
if sec == 1
    [S_p, S_u] = sections(x_p, sec, dx);
elseif sec == 2
    [S_p, S_u] = sections(x_p, sec, dx);
else
    S_p = ones(size(x_p));
    S_u = ones(size(x_u));
end
%% Solution

u = zeros(length(t_u), length(x_u));
p = zeros(length(t_p), length(x_p));

% I.C.
p_m12 = 1/gamma*phi_1(x_p);
u_0 = -S_u.*phi_0(x_u);

u(1, :) = u_0;

% B.C. (x = 0)
p(1, 1) = -gamma*dt/dx*(u(1, 1) - u_b(t_u(1)))/S_p(1) + dt/gamma*f(x_p(1), t_u(1))/S_p(1) + p_m12(1);
for j = 2:Nx
    p(1, j) = -gamma*dt/dx*(u(1, j) - u(1, j-1))/S_p(j) + dt/gamma*f(x_p(j), t_u(1))/S_p(j) + p_m12(j);
end
%B.C. (x = 1)
p(1, end) = 0;

for k = 2:Nt
    for j = 1:Nx
        u(k, j) = -gamma*S_u(j)*dt/dx*(p(k-1, j+1)-p(k-1, j)) + u(k-1, j);
    end
    % B.C. (x = 0)
    p(k, 1) = -gamma*dt/dx*(u(k, 1) - u_b(t_u(k)))/S_p(1) + dt/gamma*f(x_p(1), t_u(k))/S_p(1) + p(k-1, 1);
    for j = 2:Nx
        p(k, j) = -gamma*dt/dx*(u(k, j) - u(k, j-1))/S_p(j) + dt/gamma*f(x_p(j), t_u(k))/S_p(j) + p(k-1, j);
    end
    %B.C. (x = 1)
    p(k, end) = 0;
end

for j = 1:Nx
    u(Nt+1, j) = -gamma*S_u(j)*dt/dx*(p(Nt, j+1)-p(Nt, j)) + u(Nt, j);
end

N = 2^ceil(log2(Nt+1));
fs = 1/dt;
U_in = fft(u_b(t_u), N).';
U_out = fft(u(:, end), N);

frf = abs(U_out./U_in);
frf_norm_dB = 20*log10(frf(1:N/2)/max(frf));
freq = 0:fs/N:fs*(1-1/N);

[xxp, ttp] = meshgrid(x_p, t_p);
[xxu, ttu] = meshgrid(x_u, t_u);

P_out = 20*log10(abs(fft(p(:, 1), N)));

%% Graphs

figure
sgtitle('Pressure and Velocity')
subplot(1, 2, 1)
surf(xxp, ttp, p, 'EdgeAlpha', 0)
view(2)
title('p(x, t)')
xlabel('x')
ylabel('t', 'Rotation', 0)
colorbar
subplot(1, 2, 2)
surf(xxu, ttu, u, 'EdgeAlpha', 0)
view(2)
title('u(x, t)')
xlabel('x')
ylabel('t', 'Rotation', 0)
colorbar

figure
plot(freq(1:N/2), frf_norm_dB, 'LineWidth', 1.2)
ylim([1.1*min(frf_norm_dB) 10])
title('|FFT\{u_{out}(t)\} / FFT\{u_{in}(t)\}|')
xlabel('frequency')
ylabel('Magnitude [dB]')
grid on

figure
plot(freq(1:N/2), P_out(1:N/2), 'LineWidth', 1.2)
title('|FFT\{p_{out}(t)\}|')
xlabel('frequency')
ylabel('Magnitude [dB]')
grid on