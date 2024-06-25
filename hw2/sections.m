function [S_int_p, S_int_u] = sections(x, n, dx)

% vocal tract profile, non-dimensional [pos S] pairs
% /E/
S1 = [0 1;0.09 0.4;0.11 2.4;0.24 2.4;0.26 3.2;0.29 3.2;0.32 4.2;...
0.41 4.2;0.47 3.2;0.59 1.8;0.65 1.6;0.71 1.6;0.74 1;0.76 0.8;...
0.82 0.8;0.88 2;0.91 2;0.94 3.2;1 3.2];
% /A/
S2 = [0 1;0.03 0.60;0.09 0.4;0.12 1.6;0.18 0.6;0.29 0.2;0.35 0.4;...
    0.41 0.8;0.47 1;0.50 0.6;0.59 2;0.65 3.2;0.85 3.2;0.94 2;1 2];

if n == 1
    S = S1;
else
    S = S2;
end

S_int_p = zeros(size(x));
S_int_u = zeros(size(x));

for ii = 1:(length(S(:, 1))-1)
    m = (S(ii+1, 2) - S(ii, 2))/(S(ii+1, 1) - S(ii, 1));
    q = (S(ii+1, 1)*S(ii, 2) - S(ii, 1)*S(ii+1, 2))/(S(ii+1, 1) - S(ii, 1));
    ind1 = find(abs(x - S(ii, 1)) < 1e-6);
    ind2 = find(abs(x - S(ii+1, 1)) < 1e-6);
    x_cut_p = x(ind1:ind2);
    x_cut_u = x_cut_p + dx/2;
    S_int_p(ind1:ind2) = m.*x_cut_p + q;
    S_int_u(ind1:ind2) = m.*x_cut_u + q;
end

S_int_u = S_int_u(1:end-1);

% figure
% subplot(2, 1, 1)
% plot(S1(:, 1), -sqrt(S1(:, 2)), 'k', S1(:, 1), sqrt(S1(:, 2)), 'k', 'LineWidth', 1.2)
% xlabel('x')
% ylabel('$\sqrt{S(x)}$', Interpreter='latex')
% ylim([-2.2 2.2])
% grid on
% title('Vocal Tract Profile, "E" Sound')
% subplot(2, 1, 2)
% plot(S2(:, 1), -sqrt(S2(:, 2)), 'k', S2(:, 1), sqrt(S2(:, 2)), 'k', 'LineWidth', 1.2)
% xlabel('x')
% ylabel('$\sqrt{S(x)}$', Interpreter='latex')
% grid on
% ylim([-2.2 2.2])
% title('Vocal Tract Profile, "A" Sound')
% 
% figure
% subplot(2, 1, 1)
% stem(x(101:131), S_int_p(101:131), 'Marker', '*', 'Color', '#EDB120', 'LineWidth', 1.2, 'MarkerSize', 8)
% hold on
% stem(x(101:131), -S_int_p(101:131), 'Marker', '*', 'Color', '#EDB120', 'LineWidth', 1.2, 'MarkerSize', 8)
% hold on
% plot(S(:, 1), -S(:, 2), S(:, 1), S(:, 2), 'LineWidth', 1.2, 'Color', 'black')
% grid on
% xlim([0.1 0.13])
% xlabel('x')
% ylabel('S(x)')
% subplot(2, 1, 2)
% stem(x(101:130)+dx/2, S_int_u(101:130), 'Marker', '*', 'Color', '#EDB120', 'LineWidth', 1.2, 'MarkerSize', 8)
% hold on
% stem(x(101:130)+dx/2, -S_int_u(101:130), 'Marker', '*', 'Color', '#EDB120', 'LineWidth', 1.2, 'MarkerSize', 8)
% hold on
% plot(S(:, 1), -S(:, 2), S(:, 1), S(:, 2), 'LineWidth', 1.2, 'Color', 'black')
% grid on
% xlim([0.1 0.13])
% xlabel('x')
% ylabel('S(x)')