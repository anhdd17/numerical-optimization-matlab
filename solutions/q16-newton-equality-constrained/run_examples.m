% RUN_EXAMPLES  Cac vi du minh hoa thuat toan Newton cho quy hoach loi voi rang
%               buoc tuyen tinh A x = b: diem xuat phat chap nhan duoc va bat ky (cau 16).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q06-bisection-search'));   % bisection_min.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: ham bac hai - mot buoc Newton
%   min 1/2 ||x||^2   s.t.  x1 + x2 + x3 = 3,  x1 - x3 = 1
fprintf('==================== VI DU 1 ====================\n');
A1 = [1 1 1; 1 0 -1];  b1 = [3; 1];
f1 = @(x) 0.5 * (x' * x);  g1 = @(x) x;  h1 = @(x) eye(3);
fprintf('-- Xuat phat chap nhan duoc x0 = (2, 0, 1) --\n');
[x1a, nu1a, o1a] = newton_eq(f1, g1, h1, A1, b1, [2; 0; 1], struct('display', true));
fprintf('x* = %s, nu = %s\n', mat2str(x1a', 6), mat2str(nu1a', 6));
fprintf('-- Xuat phat bat ky x0 = (0, 0, 0) --\n');
[x1b, nu1b, o1b] = newton_eq(f1, g1, h1, A1, b1, [0; 0; 0], struct('display', true));
fprintf('x* = %s, nu = %s (giai tay: x* = (1.5, 1, 0.5), nu = (-1, -0.5))\n\n', mat2str(x1b', 6), mat2str(nu1b', 6));

%% Vi du 2: tam giai tich (analytic centering), n = 100, p = 30
%   min -sum log x_i   s.t.  A x = b
fprintf('==================== VI DU 2 ====================\n');
n = 100;  p = 30;
[I, J] = ndgrid(1:p, 1:n);
A = sin(I .* J + 0.5 * I) + 0.3 * cos(2 * J);   % du lieu tat dinh
A(1, :) = 1 + 0.5 * sin(1:n);                     % dong duong: mien chap nhan duoc bi chan
xh = 1 + 0.5 * sin((1:n)').^2;                    % diem chap nhan duoc, xh > 0
b = A * xh;
f2 = @(x) -sum(log(max(x, 0)));                   % = +Inf neu co x_i <= 0
g2 = @(x) -1 ./ x;
h2 = @(x) diag(1 ./ x.^2);
fprintf('rank A = %d\n', rank(A));
fprintf('-- (a) xuat phat chap nhan duoc x0 = xh --\n');
[x2a, nu2a, o2a] = newton_eq(f2, g2, h2, A, b, xh, struct('display', true));
fprintf('-- (b) xuat phat bat ky x0 = (1, ..., 1), ||A x0 - b|| = %.4g --\n', norm(A * ones(n, 1) - b));
[x2b, nu2b, o2b] = newton_eq(f2, g2, h2, A, b, ones(n, 1), struct('display', true));
fprintf('-- (c) xuat phat bat ky x0 = (0.1, ..., 0.1) --\n');
[x2c, nu2c, o2c] = newton_eq(f2, g2, h2, A, b, 0.1 * ones(n, 1), struct('display', true));
fprintf('Sai khac nghiem: ||x_b - x_a|| = %.1e, ||x_c - x_a|| = %.1e\n', norm(x2b - x2a), norm(x2c - x2a));
fprintf('KKT: ||grad f + A''nu|| = %.1e, ||A x - b|| = %.1e, min x = %.4f\n', ...
  norm(g2(x2a) + A' * nu2a), norm(A * x2a - b), min(x2a));
fprintf('(c): ||A x_k - b||, k = 0..4: %s; ty so lien tiep %s, 1 - t: %s\n', mat2str(o2c.feas(1:5), 5), ...
  mat2str(o2c.feas(2:4) ./ o2c.feas(1:3), 5), mat2str(1 - o2c.ts(1:3), 5));
% Doi chieu: phuong phap khu bien x = xh + F z (F: co so cua null(A)) cho cung buoc Newton
F = null(A);
x = xh;  H = h2(x);  g = g2(x);
dz = -(F' * H * F) \ (F' * g);
sol = [H, A'; A, zeros(p)] \ [-g; zeros(p, 1)];
fprintf('Buoc Newton tai xh: he KKT va khu bien (null(A), %d x %d) sai khac %.1e\n\n', size(F, 1), size(F, 2), ...
  norm(sol(1:n) - F * dz));

%% Vi du 3: cuc dai entropy tren {1, ..., 10} voi ky vong cho truoc
%   min sum x_i log x_i   s.t.  sum x_i = 1,  sum i x_i = 3.5
fprintf('==================== VI DU 3 ====================\n');
m = 10;  v = (1:m)';
A3 = [ones(1, m); v'];  b3 = [1; 3.5];
f3 = @(x) sum(x .* log(max(x, realmin))) - log(double(all(x > 0)));   % = +Inf neu co x_i <= 0
g3 = @(x) log(x) + 1;
h3 = @(x) diag(1 ./ x);
% Nghiem dang dong: x_i ~ exp(-theta i); theta cuc tieu ham doi ngau loi
%   q(theta) = log sum exp(-theta i) + 3.5 theta,   q'(theta) = 3.5 - E_theta[i]
phi  = @(th) log(sum(exp(-th * v))) + 3.5 * th;
dphi = @(th) 3.5 - sum(v .* exp(-th * v)) / sum(exp(-th * v));
[th, ~, ob] = bisection_min(phi, -5, 5, struct('method', 'derivative', 'df', dphi, 'tol', 1e-14));
xg = exp(-th * v) / sum(exp(-th * v));
fprintf('Doi ngau + chia doi (cau 6): theta = %.10f (%d buoc), x = %s\n', th, ob.iterations, mat2str(xg', 5));
u = ones(m, 1) / m;  q = [0.5; 0.5; zeros(m - 2, 1)];
x03 = 0.5 * u + 0.5 * q;                        % tron hai phan phoi: ky vong 0.5*5.5 + 0.5*1.5 = 3.5
fprintf('-- (a) xuat phat chap nhan duoc x0 = (u + q)/2 = %s --\n', mat2str(x03', 3));
[x3a, nu3a, o3a] = newton_eq(f3, g3, h3, A3, b3, x03, struct('display', true));
fprintf('-- (b) xuat phat bat ky: phan phoi deu u (ky vong 5.5) --\n');
[x3b, nu3b, o3b] = newton_eq(f3, g3, h3, A3, b3, u, struct('display', true));
fprintf('||x_a - x_doingau|| = %.1e, ||x_b - x_doingau|| = %.1e, nu = %s (nu_2 = theta: %.1e)\n', ...
  norm(x3a - xg), norm(x3b - xg), mat2str(nu3a', 8), abs(nu3a(2) - th));
fprintf('Entropy: phan phoi deu %.6f (khong thoa ky vong), nghiem %.6f, log 10 = %.6f\n\n', ...
  -sum(u .* log(u)), -sum(x3a .* log(x3a)), log(10));

%% Hinh 1: vi du 2 - toc do hoi tu
fig = figure('visible', 'off', 'Position', [100 100 1150 440]);
subplot(1, 2, 1);
semilogy(0:numel(o2a.dec) - 1, max(o2a.dec, 1e-20), 'k.-', 'LineWidth', 1.3, 'MarkerSize', 14); hold on; grid on;
semilogy(0:numel(o2b.rnorm) - 1, o2b.rnorm, 'b.-', 'LineWidth', 1.3, 'MarkerSize', 14);
semilogy(0:numel(o2c.rnorm) - 1, o2c.rnorm, 'r.-', 'LineWidth', 1.3, 'MarkerSize', 14);
xlabel('buoc lap k'); ylabel('\lambda^2/2  hoac  ||r||');
legend('(a) chap nhan duoc: \lambda(x_k)^2/2', '(b) x_0 = 1: ||r(x_k, \nu_k)||', '(c) x_0 = 0.1: ||r(x_k, \nu_k)||', ...
  'Location', 'northeast');
ylim([1e-19 1e18]); set(gca, 'YTick', 10.^(-18:6:18));
title('Vi du 2: tam giai tich, n = 100, p = 30');
subplot(1, 2, 2);
semilogy(0:o2b.iterations, max(o2b.feas, 1e-16), 'b.-', 'LineWidth', 1.3, 'MarkerSize', 14); hold on; grid on;
semilogy(0:o2c.iterations, max(o2c.feas, 1e-16), 'r.-', 'LineWidth', 1.3, 'MarkerSize', 14);
text(4.5, 1e-6, sprintf('(c): t_0..t_3 = 1/%g, 1/%g, 1/%g, %g', 1 ./ o2c.ts(1:3), o2c.ts(4)), 'FontSize', 9, 'Color', 'r');
xlabel('buoc lap k'); ylabel('||A x_k - b||');
legend('(b) x_0 = 1', '(c) x_0 = 0.1', 'Location', 'northeast');
title('Do vi pham rang buoc: ||A x_{k+1} - b|| = (1 - t_k) ||A x_k - b||');
ylim([1e-16 1e4]); set(gca, 'YTick', 10.^(-16:4:4));
print(fig, '-dpng', '-r150', fullfile('figures', 'example2_convergence.png'));
close(fig);

%% Hinh 2: vi du 3 - phan phoi entropy cuc dai
fig = figure('visible', 'off', 'Position', [100 100 760 400]);
hb = bar(v, [x03, x3a]);
set(hb(1), 'FaceColor', [0.7 0.7 0.7]);  set(hb(2), 'FaceColor', [0.3 0.6 1]);
hold on; grid on;
hl = plot(v, xg, 'r.--', 'MarkerSize', 16);
xlabel('i'); ylabel('x_i');
legend([hb(1), hb(2), hl], {'x_0 = (u + q)/2', 'Newton', 'C e^{-\theta i} (doi ngau)'}, 'Location', 'northeast');
title(sprintf('Vi du 3: entropy cuc dai voi E[i] = 3.5, \\theta = %.4f', th));
print(fig, '-dpng', '-r150', fullfile('figures', 'example3_max_entropy.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
