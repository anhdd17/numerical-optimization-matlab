% RUN_EXAMPLES  Loi giai Bai kiem tra so 1, De so 1 (phan A: cau 1-8, phan B).
% Dung lai chuong trinh cua cac cau trong bai tap lon (thu muc solutions/).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
S = fullfile('..', '..', 'solutions');
addpath(fullfile(S, 'q01-simplex'));               % simplex.m
addpath(fullfile(S, 'q04-kkt-conditions'));        % solve_nlp.m
addpath(fullfile(S, 'q05-lagrange-duality'));      % solve_lp.m
addpath(fullfile(S, 'q06-bisection-search'));      % bisection_min.m
addpath(fullfile(S, 'q07-golden-section-search')); % golden_section.m
addpath(fullfile(S, 'q08-gauss-seidel'));          % gauss_seidel_min.m, bracket_minimum.m
addpath(fullfile(S, 'q10-gradient-momentum'));     % gradient_descent.m, momentum_gradient.m
addpath(fullfile(S, 'q11-newton-method'));         % newton_method.m
addpath(fullfile(S, 'q12-quasi-newton'));          % quasi_newton.m, wolfe_line_search.m
if ~exist('figures', 'dir'), mkdir('figures'); end
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

%% Cau 1: Gauss-Seidel
%   F = 2x1^2 - 2x1x2 + 4x1x3 + 3x2^2 - 4x2x3 + 6x3^2 + 3x3 = 1/2 x'Hx + q'x
fprintf('==================== CAU 1 (Gauss-Seidel) ====================\n');
H1 = [4 -2 4; -2 6 -4; 4 -4 12];  q1 = [0; 0; 3];
f1 = @(x) 0.5 * x' * H1 * x + q1' * x;
fprintf('Tri rieng cua H: %s (xac dinh duong => F loi chat)\n', mat2str(eig(H1)', 6));
% Cuc tieu theo x_i khi co dinh cac bien con lai: dF/dx_i = 0
cm1 = @(x, i) -(q1(i) + H1(i, :) * x - H1(i, i) * x(i)) / H1(i, i);
[x1, F1, o1] = gauss_seidel_min(f1, [2; 3; 1], struct('tol', 1e-4, 'coordMin', cm1, 'display', true));
[x1b, ~, o1b] = gauss_seidel_min(f1, [2; 3; 1], struct('tol', 1e-4));
x1ex = -H1 \ q1;
fprintf('Gauss-Seidel (cong thuc): x* = %s, F* = %.6f, %d vong\n', mat2str(x1' + 0, 6), F1, o1.iterations);
fprintf('Gauss-Seidel (tim kiem mot chieu): x* = %s, %d vong\n', mat2str(x1b' + 0, 6), o1b.iterations);
fprintf('Nghiem dung -H\\q = %s = (1/3, -1/6, -5/12), F* = %.6f\n\n', mat2str(x1ex', 6), f1(x1ex));

%% Cau 2: B = (A1, A2, A3) co phai co so toi uu?
fprintf('==================== CAU 2 (co so) ====================\n');
A2 = [1 -1 2 -1; 2 1 -3 1; 1 1 1 1];  b2 = [2; 6; 7];  c2 = [2; 1; -1; -1];   % max c'x
B2 = A2(:, 1:3);
xB2 = B2 \ b2;  y2 = B2' \ c2(1:3);  d2 = c2' - y2' * A2;
fprintf('det B = %g (khac 0 => B la co so)\n', det(B2));
fprintf('x_B = B^-1 b = %s (>= 0 => co so chap nhan duoc), x = %s\n', mat2str(xB2', 6), mat2str([xB2; 0]', 6));
fprintf('y = B^-T c_B = %s, uoc luong c_j - y''A_j = %s\n', mat2str(y2', 6), mat2str(d2 + 0, 6));
fprintf('Bai toan max: moi c_j - y''A_j <= 0 => co so toi uu, F = %g\n', c2' * [xB2; 0]);
[xs2, fs2] = simplex(-c2, A2, b2, [1 2 3], struct('display', false));
fprintf('Don hinh (min -F) xuat phat tu B = {1,2,3}: x = %s, F = %g\n', mat2str(xs2' + 0, 6), -fs2);
[xl2, fl2] = solve_lp(-c2, [], [], A2, b2, zeros(4, 1));
fprintf('Doi chieu solver: x = %s, F = %g\n\n', mat2str(xl2', 6), -fl2);

%% Cau 3: Newton
fprintf('==================== CAU 3 (Newton) ====================\n');
% F = x1^4 + x2^2 + 4x1^2 + 3/2 x2^2 + 2x3^2 - x1x2 - 2x2x3 + x3x1  (x2^2 gop: 5/2 x2^2)
f3 = @(x) x(1)^4 + x(2)^2 + 4*x(1)^2 + 1.5*x(2)^2 + 2*x(3)^2 - x(1)*x(2) - 2*x(2)*x(3) + x(3)*x(1);
g3 = @(x) [4*x(1)^3 + 8*x(1) - x(2) + x(3); 5*x(2) - x(1) - 2*x(3); 4*x(3) - 2*x(2) + x(1)];
h3 = @(x) [12*x(1)^2 + 8, -1, 1; -1, 5, -2; 1, -2, 4];
fprintf('Tri rieng Hessian tai x = 0: %s (> 0; them 12x1^2 >= 0 => F loi chat)\n', mat2str(eig(h3([0;0;0]))', 6));
[x3, F3, o3] = newton_method(f3, g3, h3, [1; 0; 2], struct('tol', 1e-4, 'display', true));
x3r = x3;  x3r(abs(x3r) <= 1e-8) = 0;
fprintf('Newton: x* = %s, F* = %g, %d buoc\n\n', mat2str(x3r' + 0, 6), F3, o3.iterations);

%% Cau 4: cuc tieu dia phuong cua f tren [1, 5]
fprintf('==================== CAU 4 (mot bien) ====================\n');
f4 = @(x) 3*x + 2*x.^2 - 2*x.^3 - 2*x.^4 .* cos(x);
df4 = @(x) 3 + 4*x - 6*x.^2 - 8*x.^3 .* cos(x) + 2*x.^4 .* sin(x);
xx = linspace(1, 5, 4001);
sc = find(diff(sign(df4(xx))));
fprintf('f''(x) doi dau tren luoi tai x ~ %s; f''(1) = %.4f, f''(5) = %.4f\n', mat2str(xx(sc), 5), df4(1), df4(5));
[x4g, f4g, o4g] = golden_section(f4, 1, 5, struct('tol', 1e-4, 'display', true));
[x4b, f4b, o4b] = bisection_min(f4, 1, 2, struct('tol', 1e-4, 'method', 'derivative', 'df', df4));
fprintf('Lat cat vang tren [1,5]: x* = %.6f, f = %.6f (%d buoc)\n', x4g, f4g, o4g.iterations);
fprintf('Chia doi theo dau f'' tren [1,2]: x* = %.6f, f = %.6f (%d buoc)\n', x4b, f4b, o4b.iterations);
x4m = fzero(df4, [3 4.5]);
fprintf('Diem dung con lai x = %.6f la cuc dai dia phuong (f = %.4f); f(5) = %.4f (cuc tieu tai bien)\n\n', ...
  x4m, f4(x4m), f4(5));

figure('visible', 'off', 'Position', [100 100 1000 430]);
subplot(1, 2, 1);
h41 = plot(xx, f4(xx), 'b-', 'LineWidth', 1.5); hold on;
h42 = plot(x4g, f4g, 'ro', 'MarkerFaceColor', 'r');
h43 = plot(5, f4(5), 'ks', 'MarkerFaceColor', 'k');
h44 = plot(x4m, f4(x4m), 'g^', 'MarkerFaceColor', 'g');
xlabel('x'); ylabel('f(x)'); title('f tren [1, 5]'); grid on;
legend([h41, h42, h44, h43], {'f(x)', 'cuc tieu dia phuong x* = 1.369', 'cuc dai dia phuong', 'x = 5 (bien)'}, 'Location', 'southwest');
subplot(1, 2, 2);
xz = linspace(1, 2, 401);
plot(xz, f4(xz), 'b-', 'LineWidth', 1.5); hold on;
plot(x4g, f4g, 'ro', 'MarkerFaceColor', 'r');
text(x4g, f4g - 0.8, sprintf('x* = %.4f', x4g), 'HorizontalAlignment', 'center');
xlabel('x'); ylabel('f(x)'); title('Phong to [1, 2]'); grid on;
print('-dpng', '-r110', fullfile('figures', 'q4_function.png'));
close;

%% Cau 6: Gradient
fprintf('==================== CAU 6 (Gradient) ====================\n');
f6 = @(x) (x(2) - x(1)^2)^2 + 100*(1 - x(1))^2;
g6 = @(x) [-4*x(1)*(x(2) - x(1)^2) - 200*(1 - x(1)); 2*(x(2) - x(1)^2)];
x06 = [-0.5; 0.5];
[x6e, F6e, o6e] = gradient_descent(f6, g6, x06, struct('step', 'exact', 'tol', 1e-3, 'display', true));
[x6a, F6a, o6a] = gradient_descent(f6, g6, x06, struct('step', 'armijo', 'tol', 1e-3));
fprintf('Gradient tim kiem chinh xac: x* = %s, f = %.3e, %d buoc\n', mat2str(x6e', 6), F6e, o6e.iterations);
fprintf('Gradient Armijo:             x* = %s, f = %.3e, %d buoc\n', mat2str(x6a', 6), F6a, o6a.iterations);
fprintf('Nghiem dung x* = (1, 1), f* = 0\n\n');

figure('visible', 'off', 'Position', [100 100 1000 430]);
for sp = 1:2
  subplot(1, 2, sp);
  if sp == 1, ax = [-0.7 1.7 0.3 1.2]; lv = [1 5 20 50 100 200]; else, ax = [0.9 1.1 0.45 1.05]; lv = [0.01 0.05 0.2 0.5]; end
  [X1, X2] = meshgrid(linspace(ax(1), ax(2), 300), linspace(ax(3), ax(4), 300));
  Z = (X2 - X1.^2).^2 + 100*(1 - X1).^2;
  contour(X1, X2, Z, lv, 'LineColor', [0.7 0.7 0.7]); hold on;
  h61 = plot(o6a.X(1, :), o6a.X(2, :), 'b.-', 'MarkerSize', 6);
  h62 = plot(o6e.X(1, :), o6e.X(2, :), 'ro-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
  h63 = plot(1, 1, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y');
  axis(ax); xlabel('x_1'); ylabel('x_2'); grid on;
end
subplot(1, 2, 1);
legend([h62, h61, h63], {sprintf('tim kiem chinh xac (%d buoc)', o6e.iterations), ...
  sprintf('Armijo (%d buoc)', o6a.iterations), 'x* = (1, 1)'}, 'Location', 'northwest');
title('Gradient tu x^0 = (-0.5, 0.5)');
subplot(1, 2, 2); title('Phong to quanh x*');
print('-dpng', '-r110', fullfile('figures', 'q6_gradient_path.png'));
close;

%% Cau 7: Gradient quan tinh
fprintf('==================== CAU 7 (Gradient quan tinh) ====================\n');
% f = 1/2 x'Hx + q'x
H7 = [4 -2 -2; -2 10 1; -2 1 8];  q7 = [-2; 5; 8];
f7 = @(x) 0.5 * x' * H7 * x + q7' * x;  g7 = @(x) H7 * x + q7;
ev7 = eig(H7);  L7 = max(ev7);  mu7 = min(ev7);
a7 = 4 / (sqrt(L7) + sqrt(mu7))^2;  b7 = ((sqrt(L7) - sqrt(mu7)) / (sqrt(L7) + sqrt(mu7)))^2;
fprintf('Tri rieng H: %s; alpha = 4/(sqrt(L)+sqrt(mu))^2 = %.6f, beta = %.6f\n', mat2str(ev7', 6), a7, b7);
[x7, F7, o7] = momentum_gradient(f7, g7, [2; 1; 3], struct('alpha', a7, 'beta', b7, 'tol', 1e-4));
fprintf('%4s %-40s %14s %12s\n', 'k', 'x_k', 'f(x_k)', '||grad f||');
for k = [0:5, o7.iterations]
  fprintf('%4d %-40s %14.8g %12.3e\n', k, mat2str(o7.X(:, k + 1)', 6), o7.F(k + 1), o7.G(k + 1));
end
[~, ~, o7g] = gradient_descent(f7, g7, [2; 1; 3], struct('step', 'fixed', 'alpha', 2 / (L7 + mu7), 'tol', 1e-4));
x7ex = -H7 \ q7;
fprintf('Quan tinh: x* = %s, f* = %.6f, %d buoc (gradient buoc 2/(L+mu): %d buoc)\n', ...
  mat2str(x7', 6), F7, o7.iterations, o7g.iterations);
fprintf('Nghiem dung -H\\q = %s = (-2/9, -4/9, -1), f* = %.6f = -44/9\n\n', mat2str(x7ex', 6), f7(x7ex));

%% Cau 8: Tua Newton (BFGS)
fprintf('==================== CAU 8 (Tua Newton) ====================\n');
f8 = @(v) v(1)^4 + v(2)^4 + v(3)^4 + v(1)^2 + 5*v(2)^2 + 2*v(3)^2 - 2*v(1)*v(2) + v(1)*v(3) ...
          - 2*v(2)*v(3) + 2*v(1) - 5*v(2) + v(3);
g8 = @(v) [4*v(1)^3 + 2*v(1) - 2*v(2) + v(3) + 2;
           4*v(2)^3 + 10*v(2) - 2*v(1) - 2*v(3) - 5;
           4*v(3)^3 + 4*v(3) + v(1) - 2*v(2) + 1];
[x8, F8, o8] = quasi_newton(f8, g8, [2; 0; 1], struct('method', 'bfgs', 'tol', 1e-4, 'display', true));
[x8d, F8d, o8d] = quasi_newton(f8, g8, [2; 0; 1], struct('method', 'dfp', 'tol', 1e-4));
[x8r, F8r] = quasi_newton(f8, g8, [2; 0; 1], struct('tol', 1e-12));
fprintf('BFGS: x* = %s, F* = %.8f, %d buoc\n', mat2str(x8', 6), F8, o8.iterations);
fprintf('DFP:  x* = %s, F* = %.8f, %d buoc\n', mat2str(x8d', 6), F8d, o8d.iterations);
fprintf('BFGS tol 1e-12 (tham chieu): x* = %s, F* = %.10f\n', mat2str(x8r', 10), F8r);
H8 = 12 * diag(x8r.^2) + [2 -2 1; -2 10 -2; 1 -2 4];
fprintf('Tri rieng Hessian tai x*: %s (> 0); [2 -2 1; -2 10 -2; 1 -2 4] xac dinh duong => F loi chat\n\n', ...
  mat2str(eig(H8)', 6));

%% Cau 5: x = (0,0,9,0,8,0)
fprintf('==================== CAU 5 (phuong an) ====================\n');
A5 = [4 4 4 1 1 0; 8 6 4 3 0 1];  b5 = [44; 36];  c5 = [5; 1; 6; 2; 0; 0];   % max c'x
x5 = [0; 0; 9; 0; 8; 0];
fprintf('A x = %s (b = %s), x >= 0 => la phuong an\n', mat2str((A5 * x5)'), mat2str(b5'));
J5 = find(x5 > 0)';
fprintf('Cac cot ung voi x_j > 0: %s, hang = %d = so cot => phuong an cuc bien\n', mat2str(J5), rank(A5(:, J5)));
y5 = A5(:, J5)' \ c5(J5);  d5 = c5' - y5' * A5;
fprintf('y = %s, uoc luong c_j - y''A_j = %s (<= 0 => toi uu cho bai max), F = %g\n', ...
  mat2str(y5', 6), mat2str(d5 + 0, 6), c5' * x5);
[xs5, fs5] = simplex(-c5, A5, b5, [], struct('display', true));
fprintf('Don hinh (min -F) tu co so {x5, x6}: x = %s, F = %g\n', mat2str(xs5' + 0, 6), -fs5);
[xl5, fl5] = solve_lp(-c5, [], [], A5, b5, zeros(6, 1));
fprintf('Doi chieu solver: x = %s, F = %g\n\n', mat2str(xl5', 6), -fl5);

%% Phan B: min sum x_i ln(a + x_i), sum x_i = 1, x >= 0
fprintf('==================== PHAN B ====================\n');
% a = thang sinh cua hoc vien: thay gia tri cua ban vao day
a = 1;
fprintf('(a) g_i(t) = t ln(a+t): g'''' = (2a+t)/(a+t)^2 > 0 voi t >= 0 => f loi chat\n');
n = 2;
xP = [0.5; 0.5];
phi = @(t) t .* log(a + t) + (1 - t) .* log(a + 1 - t);
fprintf('(b) n = 2, a = %d: x* = (1/2, 1/2), f* = ln(a + 1/2) = %.10f\n', a, log(a + 0.5));
nuP = log(a + 0.5) + 0.5 / (a + 0.5);
fprintf('    nhan tu: nu* = ln(a+1/2) + 1/(2a+1) = %.10f, lambda* = (0, 0)\n', nuP);
% KKT
hk = @(t) log(a + t) + t ./ (a + t);
fprintf('    KKT: h(x_i) - nu - lambda_i = %s\n', mat2str((hk(xP) - nuP)', 3));
% Doi ngau: max g(nu, 0) theo nu (lambda = 0 vi dg/dlambda_i = -x_i <= 0 khi ...)
[nuD, negD] = golden_section(@(v) -dual_function_p(v, zeros(n, 1), a), nuP - 2, nuP + 2, struct('tol', 1e-10));
[gD, xD] = dual_function_p(nuD, zeros(n, 1), a);
fprintf('    Doi ngau: max_nu g(nu, 0) dat tai nu = %.8f, g = %.10f, x(nu) = %s\n', nuD, gD, mat2str(xD', 8));
fprintf('    Khoang cach doi ngau f* - g* = %.2e\n', log(a + 0.5) - gD);
fprintf('    g(nu*, lambda) voi lambda = (0.3, 0.7): %.6f <= f* (doi ngau yeu)\n', dual_function_p(nuP, [0.3; 0.7], a));
probB.f = @(x) sum(x .* log(a + x));
probB.g = @(x) -x;
probB.h = @(x) sum(x) - 1;
[xN, fN] = solve_nlp(probB, [0.9; 0.1]);
fprintf('    Doi chieu solver phi tuyen tu (0.9, 0.1): x = %s, f = %.10f\n', mat2str(xN', 8), fN);
fprintf('    Bang theo thang sinh a (n = 2, x* = (1/2, 1/2)):\n');
fprintf('    %3s %14s %14s %14s\n', 'a', 'f* = ln(a+1/2)', 'nu*', 'max g (so)');
for aa = 1:12
  nuA = log(aa + 0.5) + 1 / (2*aa + 1);
  [nuDA] = golden_section(@(v) -dual_function_p(v, [0; 0], aa), nuA - 2, nuA + 2, struct('tol', 1e-10));
  fprintf('    %3d %14.8f %14.8f %14.8f\n', aa, log(aa + 0.5), nuA, dual_function_p(nuDA, [0; 0], aa));
end
% n tong quat: x_i = 1/n
nG = 5;
probG.f = @(x) sum(x .* log(a + x));  probG.g = @(x) -x;  probG.h = @(x) sum(x) - 1;
[xG, fG] = solve_nlp(probG, (1:nG)' / sum(1:nG));
fprintf('    n = %d: solver x = %s, f = %.8f; ly thuyet x_i = 1/n, f* = ln(a + 1/n) = %.8f\n\n', ...
  nG, mat2str(xG', 6), fG, log(a + 1 / nG));

figure('visible', 'off', 'Position', [100 100 1000 430]);
subplot(1, 2, 1);
tt = linspace(0, 1, 401);
plot(tt, phi(tt), 'b-', 'LineWidth', 1.5); hold on;
plot(0.5, phi(0.5), 'ro', 'MarkerFaceColor', 'r');
text(0.5, phi(0.5) + 0.12, sprintf('x_1 = 1/2, f* = ln(%g) = %.4f', a + 0.5, phi(0.5)), 'HorizontalAlignment', 'center');
xlabel('x_1  (x_2 = 1 - x_1)'); ylabel('f'); title(sprintf('Bai toan goc, n = 2, a = %d', a)); grid on;
ylim([min(phi(tt)) - 0.05, max(phi(tt)) + 0.12]);
subplot(1, 2, 2);
nus = linspace(nuP - 1.5, nuP + 1.5, 301);
gs = arrayfun(@(v) dual_function_p(v, [0; 0], a), nus);
h91 = plot(nus, gs, 'b-', 'LineWidth', 1.5); hold on;
h92 = plot(nus, log(a + 0.5) * ones(size(nus)), 'k--');
h93 = plot(nuD, gD, 'ro', 'MarkerFaceColor', 'r');
xlabel('\nu'); ylabel('g(\nu, 0)'); title('Ham doi ngau'); grid on;
legend([h91, h92, h93], {'g(\nu, 0)', 'f* = ln(a + 1/2)', 'max g = f*'}, 'Location', 'south');
print('-dpng', '-r110', fullfile('figures', 'partB_primal_dual.png'));
close;
