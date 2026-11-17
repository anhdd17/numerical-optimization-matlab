% RUN_EXAMPLES  Cac vi du minh hoa thuat toan ham chan Newton (cau 17).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q16-newton-equality-constrained'));   % newton_eq.m
addpath(fullfile('..', 'q05-lagrange-duality'));              % solve_qp.m
addpath(fullfile('..', 'q04-kkt-conditions'));                % solve_nlp.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: quy hoach toan phuong (vi du 1 cua cau 15) - duong trung tam
%   min 2x1^2 + 2x2^2 - 2x1x2 - 4x1 - 6x2  s.t.  x1 + x2 <= 2, x1 + 5x2 <= 5, x >= 0
fprintf('==================== VI DU 1 ====================\n');
P1.f0    = @(x) 2*x(1)^2 + 2*x(2)^2 - 2*x(1)*x(2) - 4*x(1) - 6*x(2);
P1.grad0 = @(x) [4*x(1) - 2*x(2) - 4; 4*x(2) - 2*x(1) - 6];
P1.hess0 = @(x) [4 -2; -2 4];
P1.fi    = @(x) [x(1) + x(2) - 2; x(1) + 5*x(2) - 5; -x(1); -x(2)];
P1.jac   = @(x) [1 1; 1 5; -1 0; 0 -1];
ps1 = -222 / 31;  xs1 = [35; 24] / 31;
[x1, f1, o1] = barrier_method(P1, [0.5; 0.5], struct('display', true));
fprintf('x = %s, sai so so voi (35/31, 24/31): %.1e, f = %.10g (p* = -222/31 = %.10g)\n', ...
  mat2str(x1', 10), norm(x1 - xs1), f1, ps1);
fprintf('lambda = %s (nghiem doi ngau: (0, 32/31 = %.6f, 0, 0))\n', mat2str(o1.lambda', 6), 32/31);
fprintf('Khe thuc te f0(x*(t)) - p* so voi m/t:\n');
for k = 1:4
  fprintf('  t = %6g: f0 - p* = %.4e <= m/t = %.1e\n', o1.t(k), o1.F(k) - ps1, o1.gap(k));
end
% Diem doi ngau tu duong trung tam: q(lambda(t)) = f0(x*(t)) - m/t (kiem tra truc tiep tai t = 1)
xt = o1.X(:, 1);  lam = -1 ./ (1 * P1.fi(xt));
Lag = @(z) P1.f0(z) + lam' * P1.fi(z);
zL = [4 -2; -2 4] \ ([4; 6] - P1.jac(0)' * lam);          % cuc tieu cua ham Lagrange (bac hai)
fprintf('t = 1: lambda(t) = %s >= 0, q(lambda) = min_x L = %.8f, f0(x*(1)) - m/1 = %.8f\n\n', ...
  mat2str(lam', 5), Lag(zL), o1.F(1) - 4);
% Duong trung tam x*(t), t tu 0.05 den 1e4 (khoi dong tu diem truoc)
ts = logspace(log10(0.05), 4, 80);
CP = zeros(2, numel(ts));  z = [0.5; 0.5];
for k = 1:numel(ts)
  t = ts(k);
  Ft = @(x) t * P1.f0(x) - sum(log(max(-P1.fi(x), 0)));
  Gt = @(x) t * P1.grad0(x) + P1.jac(x)' * (1 ./ -P1.fi(x));
  Ht = @(x) t * P1.hess0(x) + P1.jac(x)' * diag(1 ./ P1.fi(x).^2) * P1.jac(x);
  z = newton_eq(Ft, Gt, Ht, zeros(0, 2), zeros(0, 1), z, struct('start', 'feasible'));
  CP(:, k) = z;
end
fprintf('Tam giai tich (t -> 0): x*(0.05) = %s\n\n', mat2str(CP(:, 1)', 5));

%% Vi du 2: rang buoc phi tuyen (vi du 1 cau 13), xuat phat khong chap nhan duoc -> pha 1
%   min (x1-3)^2 + (x2-2)^2  s.t.  x1^2 + x2^2 <= 5,  x1 + 2x2 <= 4,   x0 = (3, 3)
fprintf('==================== VI DU 2 ====================\n');
P2.f0    = @(x) (x(1) - 3)^2 + (x(2) - 2)^2;
P2.grad0 = @(x) [2 * (x(1) - 3); 2 * (x(2) - 2)];
P2.hess0 = @(x) 2 * eye(2);
P2.fi    = @(x) [x(1)^2 + x(2)^2 - 5; x(1) + 2 * x(2) - 4];
P2.jac   = @(x) [2 * x(1), 2 * x(2); 1, 2];
P2.hessc = @(x, w) 2 * w(1) * eye(2);
fprintf('x0 = (3, 3): f_i(x0) = %s (vi pham)\n', mat2str(P2.fi([3; 3])'));
[x2, f2, o2] = barrier_method(P2, [3; 3], struct('display', true));
fprintf('Pha 1: diem chap nhan chat %s, f_i = %s\n', mat2str(o2.phase1.x', 5), mat2str(P2.fi(o2.phase1.x)', 4));
fprintf('x = %s, f = %.10g, lambda = %s (cau 13: (2, 1), 2, (1/3, 2/3))\n\n', mat2str(x2', 10), f2, mat2str(o2.lambda', 6));

%% Vi du 3: pha 1 phat hien bai toan khong co diem chap nhan chat
fprintf('==================== VI DU 3 ====================\n');
% (a) vo nghiem: x1^2 + x2^2 <= 1 va x1 + x2 >= 2
P3 = P2;
P3.fi  = @(x) [x(1)^2 + x(2)^2 - 1; 2 - x(1) - x(2)];
P3.jac = @(x) [2 * x(1), 2 * x(2); -1, -1];
[~, ~, o3] = barrier_method(P3, [0; 0]);
a3 = (-2 + sqrt(28)) / 4;
fprintf('(a) exitflag = %d, pha 1: s* = %.6f tai %s; giai tay s* = 2 - 2a = %.6f, a = %.6f\n', ...
  o3.exitflag, o3.phase1.s, mat2str(o3.phase1.x', 6), 2 - 2 * a3, a3);
% (b) hai hinh tron tiep xuc (vi du 2 cau 13): co nghiem nhung khong co diem trong
P3b = P2;
P3b.fi  = @(x) [(x(1) - 1)^2 + x(2)^2 - 1; (x(1) + 1)^2 + x(2)^2 - 1];
P3b.jac = @(x) [2 * (x(1) - 1), 2 * x(2); 2 * (x(1) + 1), 2 * x(2)];
P3b.hessc = @(x, w) 2 * sum(w) * eye(2);
[~, ~, o3b] = barrier_method(P3b, [0.3; 0.4]);
fprintf('(b) exitflag = %d, pha 1: s* = %.2e (giai tay s* = 0: khong co diem Slater)\n\n', o3b.exitflag, o3b.phase1.s);

%% Vi du 4: danh muc Markowitz n = 50 (cau 14, 15), anh huong cua mu
fprintf('==================== VI DU 4 ====================\n');
n = 50;  i = (1:n)';
mu  = 0.05 + 0.1 * sin(i).^2;
Fm  = cos(i * (1:3)) / sqrt(n);
Sig = Fm * Fm' + diag(0.01 + 0.01 * (1 + cos(3 * i)));
P4.f0 = @(x) x' * Sig * x - mu' * x;  P4.grad0 = @(x) 2 * Sig * x - mu;  P4.hess0 = @(x) 2 * Sig;
P4.fi = @(x) -x;  P4.jac = @(x) -eye(n);
P4.A = ones(1, n);  P4.b = 1;
[xq4, fq4] = solve_qp(2 * Sig, -mu, [], [], ones(1, n), 1, zeros(n, 1));
fprintf('Solver qp: f* = %.12g\n', fq4);
mus = [2 10 50 200];
o4 = cell(1, numel(mus));
fprintf('%6s %10s %14s %14s %12s %12s\n', 'mu', 'buoc ngoai', 'tong Newton', 'f - f_qp', 'm/t cuoi', 'tai san');
for j = 1:numel(mus)
  [x4, f4, o] = barrier_method(P4, ones(n, 1) / n, struct('mu', mus(j), 'tol', 1e-8));
  fprintf('%6g %10d %14d %14.2e %12.1e %12d\n', mus(j), o.outer, sum(o.newton), f4 - fq4, o.gap(end), nnz(x4 > 1e-6));
  o4{j} = o;
end
fprintf('mu = 10: so buoc Newton moi buoc ngoai: %s\n', mat2str(o4{2}.newton));
fprintf('mu = 200: nu = %.8f; KKT: ||grad f + A''nu - lambda|| = %.1e, min lambda = %.1e, sum lambda_i x_i = %.1e\n\n', ...
  o.nu, norm(P4.grad0(x4) + P4.A' * o.nu - o.lambda), min(o.lambda), o.lambda' * x4);

%% Hinh 1: duong trung tam (vi du 1)
fig = figure('visible', 'off', 'Position', [100 100 760 560]);
hold on; grid on; box on;
fill([0 2 1.25 0], [0 0 0.75 1], [0.85 0.92 1], 'EdgeColor', [0.2 0.4 0.8], 'LineWidth', 1.2);
[X1, X2] = meshgrid(linspace(-0.3, 2.3, 300), linspace(-0.3, 1.6, 300));
Z = 2 * X1.^2 + 2 * X2.^2 - 2 * X1 .* X2 - 4 * X1 - 6 * X2;
contour(X1, X2, Z, [-7.1 -6.8 -6 -5 -4 -3 -2 -1], 'LineColor', [0.65 0.65 0.65]);
hc = plot(CP(1, :), CP(2, :), 'm-', 'LineWidth', 1.8);
hb = plot(o1.X(1, :), o1.X(2, :), 'ko', 'MarkerFaceColor', [0 0.6 0], 'MarkerSize', 7);
text(o1.X(1, 1) - 0.05, o1.X(2, 1) - 0.07, 't = 1', 'FontSize', 9, 'HorizontalAlignment', 'right', ...
  'BackgroundColor', 'w');
text(o1.X(1, 2) - 0.08, o1.X(2, 2) + 0.12, 't = 10, 100, ...', 'FontSize', 9, 'HorizontalAlignment', 'right', ...
  'BackgroundColor', 'w');
plot(CP(1, 1), CP(2, 1), 'm^', 'MarkerFaceColor', 'm', 'MarkerSize', 8);
plot(xs1(1), xs1(2), 'kp', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0]);
text(CP(1, 1) - 0.05, CP(2, 1) - 0.1, 'x^*(0.05)', 'FontSize', 9, 'HorizontalAlignment', 'right');
text(1.2, 0.88, 'x^* = (35/31, 24/31)', 'FontSize', 10, 'BackgroundColor', 'w');
axis equal; axis([-0.3 2.3 -0.3 1.5]);
xlabel('x_1'); ylabel('x_2');
legend([hc, hb], {'duong trung tam x^*(t), 0.05 \leq t \leq 10^4', 'cac diem trung tam cua thuat toan (\mu = 10)'}, ...
  'Location', 'northeast');
title('Vi du 1: duong trung tam');
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_central_path.png'));
close(fig);

%% Hinh 2: khe doi ngau theo so buoc Newton (vi du 4)
fig = figure('visible', 'off', 'Position', [100 100 800 460]);
cols = {'b', [0 0.6 0], 'r', 'm'};
hh = zeros(1, numel(mus));
for j = 1:numel(mus)
  o = o4{j};
  cum = cumsum(o.newton);
  % bac thang: trong buoc ngoai k (cum(k-1) -> cum(k) buoc Newton) khe con la m/t_{k-1}
  xs = [cum(1), kron(cum(2:end), [1 1])];
  ys = [kron(o.gap(1:end-1), [1 1]), o.gap(end)];
  hh(j) = semilogy(xs, ys, '-', 'Color', cols{j}, 'LineWidth', 1.4); hold on;
end
grid on;
xlabel('tong so buoc Newton'); ylabel('khe doi ngau m/t');
legend(hh, arrayfun(@(v) sprintf('\\mu = %g', v), mus, 'UniformOutput', false), 'Location', 'northeast');
title('Vi du 4: Markowitz n = 50, m = 50 rang buoc');
print(fig, '-dpng', '-r150', fullfile('figures', 'example4_gap_vs_newton.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
