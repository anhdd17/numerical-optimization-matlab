% RUN_PRACTICE2  Giai De luyen so 2 (12 cau + tu luan) trong practice.md:
% chay phuong phap de yeu cau (ham cua bai tap lon), tinh tay tung buoc
% va doi chieu voi solver (glpk/linprog, sqp/fmincon qua solve_nlp).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_practice2
clear; clc;
exam_setup;
nodisp = struct('display', false);

%% Cau 1. Don hinh
fprintf('===== Cau 1: don hinh =====\n');
% max 4x1 + 3x2:  2x1 + x2 <= 10,  x1 + x2 <= 7,  x2 <= 5
A1 = [2 1 1 0 0; 1 1 0 1 0; 0 1 0 0 1];  b1 = [10; 7; 5];  c1 = [-4; -3; 0; 0; 0];
[x, fv, ef, o] = simplex(c1, A1, b1, [3 4 5]);
fprintf('x* = %s, F* = %g, %d buoc, exitflag %d\n', mat2str(x'), -fv, o.iterations, ef);
xs = solve_lp(c1(1:2), A1(:, 1:2), b1, [], [], [0; 0]);
fprintf('solve_lp: %s\n', mat2str(xs', 6));

%% Cau 2. Hai pha
fprintf('===== Cau 2: hai pha =====\n');
% min x1 - 2x2:  x1 + x2 >= 2,  -x1 + x2 <= 1
A2 = [1 1 -1 0; -1 1 0 1];  b2 = [2; 1];  c2 = [1; -2; 0; 0];
[x, fv, ef, o] = two_phase_simplex(c2, A2, b2, nodisp);
fprintf('so bien gia %d, w* = %g, exitflag %d (-3 = khong bi chan)\n', numel(o.artificial), o.phase1.w, ef);

%% Cau 3. Bai toan M
fprintf('===== Cau 3: bai toan M =====\n');
% min 2x1 + 3x2 + x3:  x1 + x2 + x3 = 5,  2x1 + x2 >= 4
A3 = [1 1 1 0; 2 1 0 -1];  b3 = [5; 4];  c3 = [2; 3; 1; 0];
[x, fv, ef, o] = big_m_simplex(c3, A3, b3, nodisp);
fprintf('so bien gia %d, x* = %s, F* = %g, exitflag %d\n', numel(o.artificial), mat2str(x'), fv, ef);
xs = solve_lp(c3(1:3), -[2 1 0], -4, [1 1 1], 5, [0; 0; 0]);
fprintf('solve_lp: %s\n', mat2str(xs', 6));

%% Cau 4. KKT
fprintf('===== Cau 4: KKT =====\n');
% min (x1-2)^2 + (x2-1)^2:  x1^2 - x2 <= 0,  x1 + x2 <= 2
p4.f = @(x) (x(1) - 2)^2 + (x(2) - 1)^2;  p4.g = @(x) [x(1)^2 - x(2); x(1) + x(2) - 2];
x4 = [1; 1];  gf = [2*(x4(1) - 2); 2*(x4(2) - 1)];  Jg = [2*x4(1), -1; 1, 1];
lam = -(Jg') \ gf;
fprintf('tai (1,1): g = %s, lambda = %s, rank Jg = %d\n', mat2str(p4.g(x4)'), mat2str(lam', 6), rank(Jg));
[xs, fs] = solve_nlp(p4, [0; 0]);
fprintf('solve_nlp: x = %s, f = %.6f\n', mat2str(xs', 6), fs);

%% Cau 5. Doi ngau Lagrange
fprintf('===== Cau 5: doi ngau =====\n');
% min x1^2 + 2x2^2:  x1 + x2 >= 3  ->  g(l) = 3l - 3l^2/8
l = linspace(0, 8, 8001);  gl = 3*l - 3*l.^2/8;  [gm, k] = max(gl);
fprintf('max g = %.6f tai lambda = %.4f\n', gm, l(k));
p5.f = @(x) x(1)^2 + 2*x(2)^2;  p5.g = @(x) 3 - x(1) - x(2);
[xs, fs] = solve_nlp(p5, [0; 0]);
fprintf('solve_nlp: x = %s, f = %.6f\n', mat2str(xs', 6), fs);

%% Cau 6. Chia doi theo dao ham
fprintf('===== Cau 6: chia doi =====\n');
f6 = @(x) x.^4 - 3*x.^2 + x;  df6 = @(x) 4*x.^3 - 6*x + 1;
a = 0; b = 2;
for k = 1:3
  c = (a + b)/2;  fprintf('buoc %d: c = %.4f, f''(c) = %+.6f', k, c, df6(c));
  if df6(c) > 0, b = c; else, a = c; end
  fprintf(' -> [%.4f, %.4f]\n', a, b);
end
[x6, ~, o6] = bisection_min(f6, 0, 2, struct('method', 'derivative', 'df', df6, 'tol', 1e-4));
fprintf('bisection_min: x* = %.6f (%d buoc); fzero f'' = %.6f; f'' tren [0,2] doi dau %d lan\n', ...
  x6, o6.iterations, fzero(df6, [0.5 2]), sum(diff(sign(df6(linspace(0, 2, 2001)))) ~= 0));

%% Cau 7. Hooke-Jeeves
fprintf('===== Cau 7: Hooke-Jeeves =====\n');
H7 = [4 1; 1 2];  q7 = [-3; -1];  f7 = @(x) 0.5*x'*H7*x + q7'*x;
fprintf('F(0,0) = %g, F(0.5,0) = %g, F(0.5,0.5) = %g, F(0.5,-0.5) = %g\n', f7([0;0]), f7([0.5;0]), f7([0.5;0.5]), f7([0.5;-0.5]));
[x7, F7, o7] = hooke_jeeves(f7, [0; 0], struct('delta', 0.5, 'tol', 1e-6));
fprintf('hooke_jeeves: x = %s, F = %.6f, %d lan tinh F; -H\\q = %s\n', mat2str(x7', 6), F7, o7.nfev, mat2str((-H7\q7)', 6));

%% Cau 8. Do nuoc (water-filling)
fprintf('===== Cau 8: water-filling =====\n');
[x8, nu8, lam8, lev8] = water_filling([0.5; 1; 2], 1);
fprintf('x* = %s, nu = %.4f, muc nuoc = %.4f, lambda = %s, f* = %.6f\n', mat2str(x8', 6), nu8, lev8, mat2str(lam8', 6), -sum(log([0.5; 1; 2] + x8)));

%% Cau 9. Frank-Wolfe
fprintf('===== Cau 9: Frank-Wolfe =====\n');
f9 = @(x) (x(1) - 2)^2 + (x(2) - 2)^2;  g9 = @(x) [2*(x(1) - 2); 2*(x(2) - 2)];
V = [0 4 0; 0 0 2];  x = [0; 0];
for k = 1:3
  [~, j] = min(g9(x)' * V);  y = V(:, j);  d = y - x;
  gam = min(1, max(0, -(g9(x)'*d) / (2*(d'*d))));
  fprintf('buoc %d: y = %s, gap = %.4f, gamma = %.4f', k, mat2str(y'), g9(x)'*(x - y), gam);
  x = x + gam*d;  fprintf(' -> x = %s, f = %.4f\n', mat2str(x', 6), f9(x));
end
% FW hoi tu cham (zic-zac giua hai dinh, sai so ~ 1/k): dem so buoc den khe <= tol
for tol9 = [1e-2 1e-3 1e-4]
  [x9, f9v, o9] = frank_wolfe(f9, g9, [1 2], 4, [], [], [0; 0], [0; 0], struct('tol', tol9, 'maxIter', 1e5));
  fprintf('frank_wolfe tol %.0e: x = %s, f = %.6f, %d buoc\n', tol9, mat2str(x9', 6), f9v, o9.iterations);
end
[x9, f9v, o9] = frank_wolfe(f9, g9, [1 2], 4, [], [], [0; 0], [0; 0], struct('tol', 1e-4, 'variant', 'away', 'maxIter', 1e5));
fprintf('frank_wolfe away tol 1e-4: x = %s, f = %.6f, %d buoc; nghiem dung (1.6, 1.2), f* = 0.8\n', mat2str(x9', 6), f9v, o9.iterations);

%% Cau 10. Rosen
fprintf('===== Cau 10: Rosen =====\n');
% min (x1-2)^2 + (x2-2)^2:  x1 + x2 <= 2,  -x1 <= 0,  -x2 <= 0,  x0 = (2, 0)
g10 = @(x) [2*(x(1) - 2); 2*(x(2) - 2)];  x0 = [2; 0];
N = [1 1; 0 -1];  u = -(N*N') \ (N*g10(x0));  P = eye(2) - N'*((N*N')\N);
fprintf('W = {1, 3}: P = %s, u = %s\n', mat2str(P), mat2str(u'));
N = [1 1];  P = eye(2) - N'*((N*N')\N);  d = -P*g10(x0);
fprintf('W = {1}: P = %s, d = %s\n', mat2str(P), mat2str(d'));
p10.f = @(x) (x(1) - 2)^2 + (x(2) - 2)^2;  p10.grad = g10;
p10.g = @(x) [x(1) + x(2) - 2; -x(1); -x(2)];  p10.jacg = @(x) [1 1; -1 0; 0 -1];
[x10, f10, o10] = gradient_projection(p10, x0, struct('lineSearch', 'exact', 'display', true));
fprintf('gradient_projection: x = %s, f = %.6f\n', mat2str(x10', 6), f10);

%% Cau 11. Newton co rang buoc dang thuc
fprintf('===== Cau 11: Newton Ax = b =====\n');
f11 = @(x) exp(x(1)) + exp(2*x(2)) + x(3)^2;
g11 = @(x) [exp(x(1)); 2*exp(2*x(2)); 2*x(3)];
h11 = @(x) diag([exp(x(1)), 4*exp(2*x(2)), 2]);
A11 = [1 1 1];  b11 = 1;  x0 = [1; 0; 0];
K = [h11(x0), A11'; A11, 0];  s = K \ [-g11(x0); 0];  dx = s(1:3);
fprintf('buoc Newton tai x0: dx = %s, w = %.6f, lambda^2/2 = %.6f\n', mat2str(dx', 6), s(4), dx'*h11(x0)*dx/2);
[x11, nu11, o11] = newton_eq(f11, g11, h11, A11, b11, x0, struct('tol', 1e-12));
fprintf('newton_eq: x* = %s, nu = %.6f, f* = %.6f, %d buoc\n', mat2str(x11', 6), nu11, f11(x11), o11.iterations);
fprintf('KKT: grad f + nu = %s\n', mat2str((g11(x11) + nu11)', 3));

%% Cau 12. Beale
fprintf('===== Cau 12: Beale =====\n');
% min x1^2 + 2x2^2 - 2x1 - 6x2:  x1 + x2 + x3 = 2,  x >= 0
c12 = [-2; -6; 0];  D12 = diag([2 4 0]);
fprintf('buoc 1 (x2 vao): theta1 = %g (x3 = 2 - x2 = 0), theta2 = %g (4x2 - 6 = 0)\n', 2, 6/4);
[x12, f12, ef12, o12] = beale_qp(c12, D12, [1 1 1], 2, struct('B0', 3, 'display', true));
fprintf('beale_qp: x* = %s, f* = %.6f, exitflag %d\n', mat2str(x12', 6), f12, ef12);
p12.f = @(x) x(1)^2 + 2*x(2)^2 - 2*x(1) - 6*x(2);  p12.g = @(x) [x(1) + x(2) - 2; -x(1); -x(2)];
[xs, fs] = solve_nlp(p12, [0; 0]);
fprintf('solve_nlp: x = %s, f = %.6f\n', mat2str(xs', 6), fs);

%% Phan B. Ham chan + doi ngau
fprintf('===== Phan B =====\n');
% min x1^2 + x2^2  s.t.  a - x1 - x2 <= 0
for a = [2 7]
  for t = [1 10 100]
    s = (t*a + sqrt(t^2*a^2 + 4*t)) / (4*t);  lt = 1/(t*(2*s - a));
    fprintf('a = %d, t = %3d: x*(t) = (%.6f, %.6f), lambda(t) = %.6f, f0 - g = %.2e (f0 - g - 1/t = %.1e)\n', ...
      a, t, s, s, lt, 2*s^2 - (a*lt - lt^2/2), 2*s^2 - (a*lt - lt^2/2) - 1/t);
  end
  pb.f0 = @(x) x(1)^2 + x(2)^2;  pb.grad0 = @(x) 2*x;  pb.hess0 = @(x) 2*eye(2);
  pb.fi = @(x) a - x(1) - x(2);  pb.jac = @(x) [-1 -1];
  [xb, fb, ob] = barrier_method(pb, [a; a], struct('t0', 1, 'mu', 10, 'tol', 1e-6));
  fprintf('barrier_method a = %d: x = %s, f = %.6f (a^2/2 = %g), %d buoc ngoai\n', a, mat2str(xb', 8), fb, a^2/2, ob.outer);
end
fprintf('ly thuyet: ceil(log(m/(eps t0))/log mu) = %d lan tang t, cong lan dinh tam dau tai t0 => 7 buoc ngoai\n', ceil(log(1/1e-6)/log(10)));
