% RUN_PRACTICE  (1) Kiem tra bo ham toolkit/ tren De so 1 (so voi loi giai trong
% exams/test1-paper1); (2) giai De luyen trong tai lieu on thi, doi chieu voi
% nghiem dung va solver.
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_practice
clear; clc;
addpath('toolkit');
addpath(fullfile('..', '..', 'solutions', 'q04-kkt-conditions'));   % solve_nlp.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% ===================== PHAN 1: bo ham tren De so 1 =====================
fprintf('############ DE SO 1 (kiem tra toolkit) ############\n');
fprintf('--- Cau 1 ---\n');
x = gs_quad([4 -2 4; -2 6 -4; 4 -4 12], [0; 0; 3], [2; 3; 1], 1e-4);
fprintf('--- Cau 2 ---\n');
basis_check([1 -1 2 -1; 2 1 -3 1; 1 1 1 1], [2; 6; 7], [2 1 -1 -1], [1 2 3], 'max');
fprintf('--- Cau 3 ---\n');
x = newton(@(x) [4*x(1)^3 + 8*x(1) - x(2) + x(3); 5*x(2) - x(1) - 2*x(3); 4*x(3) - 2*x(2) + x(1)], ...
           @(x) [12*x(1)^2 + 8, -1, 1; -1, 5, -2; 1, -2, 4], [1; 0; 2], 1e-4);
fprintf('--- Cau 4 ---\n');
x = golden(@(x) 3*x + 2*x.^2 - 2*x.^3 - 2*x.^4 .* cos(x), 1, 5, 1e-4);
fprintf('--- Cau 5 ---\n');
point_check([4 4 4 1 1 0; 8 6 4 3 0 1], [44; 36], [5 1 6 2 0 0], [0 0 9 0 8 0], 'max');
fprintf('--- Cau 6 ---\n');
x = grad_desc(@(x) (x(2) - x(1)^2)^2 + 100*(1 - x(1))^2, ...
              @(x) [-4*x(1)*(x(2) - x(1)^2) - 200*(1 - x(1)); 2*(x(2) - x(1)^2)], [-0.5; 0.5], 1e-3);
fprintf('--- Cau 7 ---\n');
H = [4 -2 -2; -2 10 1; -2 1 8];  q = [-2; 5; 8];  e = eig(H);  L = max(e);  m = min(e);
x = heavy_ball(@(x) H*x + q, [2; 1; 3], 4/(sqrt(L) + sqrt(m))^2, ((sqrt(L) - sqrt(m))/(sqrt(L) + sqrt(m)))^2, 1e-4);
fprintf('--- Cau 8 ---\n');
f8 = @(v) v(1)^4 + v(2)^4 + v(3)^4 + v(1)^2 + 5*v(2)^2 + 2*v(3)^2 - 2*v(1)*v(2) + v(1)*v(3) - 2*v(2)*v(3) + 2*v(1) - 5*v(2) + v(3);
g8 = @(v) [4*v(1)^3 + 2*v(1) - 2*v(2) + v(3) + 2; 4*v(2)^3 + 10*v(2) - 2*v(1) - 2*v(3) - 5; 4*v(3)^3 + 4*v(3) + v(1) - 2*v(2) + 1];
x = bfgs(f8, g8, [2; 0; 1], 1e-4);

%% ===================== PHAN 2: De luyen =====================
fprintf('\n############ DE LUYEN ############\n');
fprintf('--- Cau 1: Gauss-Seidel ---\n');
H1 = [6 2 0; 2 4 -2; 0 -2 8];  q1 = [1; 0; -4];
fprintf('Dinh thuc con chinh: %g %g %g\n', H1(1,1), det(H1(1:2,1:2)), det(H1));
x1 = gs_quad(H1, q1, [1; 1; 1], 1e-4);
x1b = gs_min(@(x) 0.5*x'*H1*x + q1'*x, [1; 1; 1], 1e-4);
fprintf('Nghiem dung -H\\q = %s, F* = %.8f\n', mat2str((-H1\q1)', 8), 0.5*(-H1\q1)'*H1*(-H1\q1) + q1'*(-H1\q1));

fprintf('--- Cau 2: co so ---\n');
basis_check([1 1 1 0; 1 -1 0 1], [2; 4], [2 1 1 0], [1 2], 'max');
basis_check([1 1 1 0; 1 -1 0 1], [2; 4], [2 1 1 0], [1 4], 'max');

fprintf('--- Cau 3: Newton ---\n');
g3 = @(x) [4*x(1)^3 + 2*x(1) + x(2) - 4; x(1) + 4*x(2) - x(3); -x(2) + 2*x(3) + 2];
h3 = @(x) [12*x(1)^2 + 2, 1, 0; 1, 4, -1; 0, -1, 2];
fprintf('Dinh thuc con chinh phan hang Hessian: %g %g %g\n', 2, det([2 1; 1 4]), det(h3([0;0;0])));
x3 = newton(g3, h3, [1; 1; 1], 1e-4);

fprintf('--- Cau 4: mot bien ---\n');
f4 = @(x) exp(-x) + x.^2 - 2*sin(3*x);
df4 = @(x) -exp(-x) + 2*x - 6*cos(3*x);
fprintf('f''(0) = %.4f, f''(2) = %.4f\n', df4(0), df4(2));
r1 = fzero(df4, [0.2 1]);  r2 = fzero(df4, [1.2 1.95]);
fprintf('Diem dung: %.6f (f = %.6f), %.6f (f = %.6f); f(2) = %.6f\n', r1, f4(r1), r2, f4(r2), f4(2));
x4 = golden(f4, 0, 2, 1e-4);
x4b = fminbnd(f4, 0, 2);
fprintf('fminbnd tren [0, 2]: %.6f\n', x4b);
xx = linspace(0, 2, 401);
fig = figure('visible', 'off', 'Position', [100 100 700 380]);
h1 = plot(xx, f4(xx), 'b-', 'LineWidth', 1.5); hold on;
h2 = plot(r1, f4(r1), 'ro', 'MarkerFaceColor', 'r');
h3p = plot(r2, f4(r2), 'g^', 'MarkerFaceColor', 'g');
h4 = plot(2, f4(2), 'ks', 'MarkerFaceColor', 'k');
grid on; xlabel('x'); ylabel('f(x)');
legend([h1, h2, h3p, h4], {'f(x) = e^{-x} + x^2 - 2 sin 3x', 'cuc tieu dia phuong', 'cuc dai dia phuong', 'x = 2 (bien)'}, 'Location', 'northwest');
print(fig, '-dpng', '-r110', fullfile('figures', 'practice_q4.png'));
close(fig);

fprintf('--- Cau 5: phuong an ---\n');
point_check([1 1 1 1 0; 2 1 3 0 1], [6; 9], [1 2 1 1 0], [1 2 1 2 2], 'max');
point_check([1 1 1 1 0; 2 1 3 0 1], [6; 9], [1 2 1 1 0], [0 6 0 0 3], 'max');

fprintf('--- Cau 6: gradient ---\n');
f6 = @(x) (x(1) - 1)^2 + 10*(x(2) - x(1)^2)^2;
g6 = @(x) [2*(x(1) - 1) - 40*x(1)*(x(2) - x(1)^2); 20*(x(2) - x(1)^2)];
x6 = grad_desc(f6, g6, [0; 0], 1e-3);

fprintf('--- Cau 7: gradient quan tinh ---\n');
H7 = [6 1 0; 1 4 -1; 0 -1 2];  q7 = [2; -1; 3];  e7 = eig(H7);  L7 = max(e7);  m7 = min(e7);
a7 = 4/(sqrt(L7) + sqrt(m7))^2;  b7 = ((sqrt(L7) - sqrt(m7))/(sqrt(L7) + sqrt(m7)))^2;
fprintf('Tri rieng %s, alpha = %.6f, beta = %.6f\n', mat2str(e7', 6), a7, b7);
x7 = heavy_ball(@(x) H7*x + q7, [1; 1; 1], a7, b7, 1e-4);
fprintf('Nghiem dung -H\\q = %s\n', mat2str((-H7\q7)', 8));

fprintf('--- Cau 8: tua Newton ---\n');
f8p = @(v) v(1)^4 + v(2)^4 + 2*v(1)^2 + 3*v(2)^2 + v(3)^2 - v(1)*v(2) + v(2)*v(3) - 2*v(1) + v(2) - 3*v(3);
g8p = @(v) [4*v(1)^3 + 4*v(1) - v(2) - 2; 4*v(2)^3 + 6*v(2) - v(1) + v(3) + 1; 2*v(3) + v(2) - 3];
x8 = bfgs(f8p, g8p, [1; 1; 1], 1e-4);
x8r = fminsearch(f8p, [0; 0; 0], optimset('TolX', 1e-10, 'TolFun', 1e-14, 'MaxFunEvals', 1e5, 'MaxIter', 1e5));
fprintf('fminsearch (doi chieu): %s, F = %.8f\n', mat2str(x8r', 6), f8p(x8r));

fprintf('--- Phan B ---\n');
% min (x1 - a/4)^2 + sum_{i>=2} x_i^2,  sum x = 1, x >= 0
for a = [2 8]
  for n = [2 4]
    c = [a/4; zeros(n - 1, 1)];
    if a <= 4, xt = [a/4 + (1 - a/4)/n; (1 - a/4)/n * ones(n - 1, 1)]; else, xt = [1; zeros(n - 1, 1)]; end
    nu = 2*(xt(1) - c(1));  lam = 2*(xt - c) - nu;
    pb.f = @(x) sum((x - c).^2);  pb.g = @(x) -x;  pb.h = @(x) sum(x) - 1;
    [xs, fs] = solve_nlp(pb, ones(n, 1) / n);
    % doi ngau: g(nu, lambda) = nu + sum(-mu_i c_i - mu_i^2/4), mu = nu + lambda
    mu = nu + lam;  gd = nu + sum(-mu .* c - mu.^2 / 4);
    fprintf('a = %d, n = %d: x* = %s, f* = %.6f, nu = %.4f, lambda = %s, g(nu,lambda) = %.6f; solver x = %s, f = %.6f\n', ...
      a, n, mat2str(xt', 4), sum((xt - c).^2), nu, mat2str(lam', 4), gd, mat2str(xs', 4), fs);
  end
end
