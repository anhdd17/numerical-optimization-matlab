% RUN_EXAMPLES  Cac vi du minh hoa giai thuat lat cat vang (cau 7).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q06-bisection-search'));   % bisection_min.m (so sanh)
if ~exist('figures', 'dir'), mkdir('figures'); end
tau = (sqrt(5) - 1) / 2;

%% Vi du 1: ham bac 4 tren [0, 2] (vi du 1 cua cau 6)
fprintf('==================== VI DU 1 ====================\n');
f1 = @(x) x.^4 - 14*x.^3 + 60*x.^2 - 70*x;
[x1, fx1, o1] = golden_section(f1, 0, 2, struct('tol', 1e-5, 'display', true));
kTh = ceil(log(1e-5 / 2) / log(tau));
fprintf('x* = %.8f, f* = %.8f, so buoc = %d (ly thuyet %d), so lan tinh f = %d\n', ...
  x1, fx1, o1.iterations, kTh, o1.nfev);
[~, ~, ob] = bisection_min(f1, 0, 2, struct('tol', 1e-5));
fprintf('Chia doi (cau 6) cung tol: %d buoc, %d lan tinh f\n\n', ob.iterations, ob.nfev);

%% Vi du 2: so sanh voi chia doi tren f(x) = e^{-x} + x^2, [0, 1]
fprintf('==================== VI DU 2 ====================\n');
f2  = @(x) exp(-x) + x.^2;
df2 = @(x) -exp(-x) + 2*x;
tol = 1e-6;
xr2 = fzero(df2, [0 1]);
[xg, ~, og] = golden_section(f2, 0, 1, struct('tol', tol));
[xa, ~, oa] = bisection_min(f2, 0, 1, struct('tol', tol));
[xd, ~, od] = bisection_min(f2, 0, 1, struct('tol', tol, 'method', 'derivative', 'df', df2));
fprintf('%-26s %12s %6s %10s %12s\n', 'Phuong phap', 'x*', 'buoc', 'tinh f/f''', '|x - x_ref|');
fprintf('%-26s %12.8f %6d %10d %12.1e\n', 'Lat cat vang', xg, og.iterations, og.nfev, abs(xg - xr2));
fprintf('%-26s %12.8f %6d %10d %12.1e\n', 'Chia doi (gia tri ham)', xa, oa.iterations, oa.nfev, abs(xa - xr2));
fprintf('%-26s %12.8f %6d %10d %12.1e\n', 'Chia doi (dao ham)', xd, od.iterations, od.ndfev, abs(xd - xr2));
fprintf('Nghiem tham chieu: x* = %.8f\n\n', xr2);

%% Vi du 3: tim kiem theo tia (line search) cho ham hai bien
%   F(x) = (x1 - 2)^4 + (x1 - 2x2)^2,  x0 = (0, 3),  d = -grad F(x0)
%   phi(t) = F(x0 + t d),  t in [0, 1]
fprintf('==================== VI DU 3 ====================\n');
F    = @(x) (x(1) - 2)^4 + (x(1) - 2*x(2))^2;
grad = @(x) [4*(x(1) - 2)^3 + 2*(x(1) - 2*x(2)); -4*(x(1) - 2*x(2))];
x0 = [0; 3];
d  = -grad(x0);
phi = @(t) F(x0 + t * d);
[t3, p3, o3] = golden_section(phi, 0, 1, struct('tol', 1e-8));
[tr, pr] = fminbnd(phi, 0, 1, optimset('TolX', 1e-12));
xnew = x0 + t3 * d;
fprintf('F(x0) = %g, d = -grad F(x0) = %s\n', F(x0), mat2str(d'));
fprintf('Lat cat vang: t* = %.8f, phi(t*) = %.8f (%d lan tinh F)\n', t3, p3, o3.nfev);
fprintf('fminbnd:      t* = %.8f, phi(t*) = %.8f\n', tr, pr);
fprintf('Diem moi x1 = x0 + t* d = %s, F(x1) = %.6f\n\n', mat2str(xnew', 6), F(xnew));

%% Vi du 4: gioi han do chinh xac cua phuong phap chi dung gia tri ham
%   f(x) = (x - 1)^2 + 1 tren [0, 3], x* = 1
fprintf('==================== VI DU 4 ====================\n');
f4  = @(x) (x - 1).^2 + 1;
df4 = @(x) 2*(x - 1);
fprintf('%8s | %18s | %18s\n', 'tol', 'sai so lat cat vang', 'sai so chia doi f''');
for t = [1e-4 1e-6 1e-8 1e-10 1e-12]
  xg4 = golden_section(f4, 0, 3, struct('tol', t));
  xd4 = bisection_min(f4, 0, 3, struct('tol', t, 'method', 'derivative', 'df', df4));
  fprintf('%8.0e | %18.2e | %18.2e\n', t, abs(xg4 - 1), abs(xd4 - 1));
end
fprintf('sqrt(eps) = %.2e\n\n', sqrt(eps));

%% Hinh 1: vi du 1 - cac khoang va diem thu; diem rong = diem dung lai tu buoc truoc
fig = figure('visible', 'off', 'Position', [100 100 900 470]);
hold on; grid on; box on;
K = 5;
for k = 0:K
  a = o1.a(k+1);  b = o1.b(k+1);  L = b - a;
  p1 = b - tau * L;  p2 = a + tau * L;
  plot([a b], [k k], 'k-', 'LineWidth', 3);
  if k == 0
    newPt = [p1 p2];  oldPt = [];
  elseif a ~= o1.a(k)                       % giu [x1, b]: x2 cu -> x1 moi
    newPt = p2;  oldPt = p1;
  else                                      % giu [a, x2]: x1 cu -> x2 moi
    newPt = p1;  oldPt = p2;
  end
  plot(newPt, k * ones(size(newPt)), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
  if ~isempty(oldPt)
    plot(oldPt, k, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
  end
  if k == 0
    text(1.55, 0.42, sprintf('L_0 = %.4f', L), 'FontSize', 9);
  else
    text(b + 0.03, k, sprintf('L_%d = %.4f', k, L), 'FontSize', 9);
  end
end
plot([x1 x1], [-0.5 K + 0.5], 'r--');
h1 = plot(nan, nan, 'bo', 'MarkerFaceColor', 'b');
h2 = plot(nan, nan, 'bo', 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
legend([h1 h2], {'diem tinh f moi', 'diem dung lai tu buoc truoc'}, 'Location', 'southeast');
set(gca, 'YDir', 'reverse');
axis([0 2 -0.5 K + 0.5]);
xlabel('x'); ylabel('buoc k');
title('Vi du 1: lat cat vang - moi buoc chi tinh f tai mot diem moi');
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_point_reuse.png'));
close(fig);

%% Hinh 2: vi du 2 - do dai khoang theo so lan tinh ham
fig = figure('visible', 'off');
Lg = og.b - og.a;  La = oa.b - oa.a;  Ld = od.b - od.a;
semilogy(2 + (0:numel(Lg) - 1), Lg, 'k^-', 'MarkerSize', 4); hold on; grid on;   % sau k buoc: 2 + k lan tinh f
semilogy(0:2:2*(numel(La) - 1), La, 'bo-', 'MarkerSize', 4);
semilogy(0:numel(Ld) - 1, Ld, 'rs-', 'MarkerSize', 4);
semilogy([0 45], [tol tol], 'k:');
text(22, 6e-7, 'tol = 10^{-6}', 'FontSize', 10);
xlim([0 45]);
xlabel('so lan tinh f (hoac f'')'); ylabel('do dai khoang L');
legend('Lat cat vang (he so 0.618 / lan tinh f)', 'Chia doi, gia tri ham (0.707 / lan)', ...
  'Chia doi, dao ham (0.5 / lan tinh f'')', 'Location', 'northeast');
title('Vi du 2: so sanh toc do hoi tu');
print(fig, '-dpng', '-r150', fullfile('figures', 'example2_comparison.png'));
close(fig);

%% Hinh 3: vi du 3 - tim kiem theo tia tren duong muc cua F
fig = figure('visible', 'off', 'Position', [100 100 1100 450]);
subplot(1, 2, 1); hold on; grid on; box on;
[X1, X2] = meshgrid(linspace(-0.5, 3.5, 300), linspace(0, 3.5, 300));
Z = (X1 - 2).^4 + (X1 - 2*X2).^2;
contour(X1, X2, Z, [0.1 0.37 1 3 8 16 30 52], 'LineColor', [0.6 0.6 0.6]);
plot([x0(1) x0(1) + 0.08 * d(1)], [x0(2) x0(2) + 0.08 * d(2)], 'b-', 'LineWidth', 1.2);
plot(x0(1), x0(2), 'ko', 'MarkerFaceColor', 'k');
plot(xnew(1), xnew(2), 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
plot(2, 1, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
text(x0(1) + 0.08, x0(2) + 0.12, 'x_0 = (0, 3)', 'FontSize', 10);
plot([xnew(1) 2.75], [xnew(2) + 0.08 2.15], 'k-');
text(2.2, 2.25, sprintf('x_1 = (%.3f, %.3f)', xnew(1), xnew(2)), 'FontSize', 10, 'BackgroundColor', 'w');
text(2.1, 0.75, 'x^* = (2, 1)', 'FontSize', 10, 'BackgroundColor', 'w');
text(0.9, 2.75, 'tia x_0 + t d', 'Color', 'b', 'FontSize', 10, 'BackgroundColor', 'w');
axis equal; axis([-0.5 3.5 0 3.5]);
xlabel('x_1'); ylabel('x_2');
title('Duong muc cua F va tia x_0 + t d');
subplot(1, 2, 2); hold on; grid on; box on;
tt = linspace(0, 0.12, 300);
pv = arrayfun(phi, tt);
plot(tt, pv, 'b-', 'LineWidth', 1.6);
plot(t3, p3, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(t3 - 0.012, p3 + 8, sprintf('t^* = %.4f', t3), 'FontSize', 10);
ylim([0 60]);
xlabel('t'); ylabel('\phi(t) = F(x_0 + t d)');
title('Ham mot bien \phi(t) can cuc tieu');
print(fig, '-dpng', '-r150', fullfile('figures', 'example3_line_search.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
