% CHECK_UNCONSTRAINED  Tinh dap an cac cau hoi mau cua topics/line-search.md
% (cau 6, 7) va topics/unconstrained.md (cau 8-12).
%   >> check_unconstrained
exam_setup;

%% Cau 6. Chia doi
fprintf('===== Cau 6 =====\n');
L0 = 2;  dl = 1e-3;  ep = 0.01;
fprintf('6.1 so buoc chia doi (gia tri ham): log2((L0-2d)/(eps-2d)) = %.4f\n', log2((L0 - 2*dl)/(ep - 2*dl)));
fh = @(x) x.^2 - sin(x);  dfh = @(x) 2*x - cos(x);
a = 0; b = 1;
for k = 1:2
  c = (a + b)/2;  fprintf('6.2 buoc %d: c = %.4f, f''(c) = %.6f\n', k, c, dfh(c));
  if dfh(c) > 0, b = c; else, a = c; end
end
fprintf('    khoang sau 2 buoc [%.4f, %.4f]; x* = %.6f (fzero f'')\n', a, b, fzero(dfh, [0 1]));
f6 = @(x) x.^4 - 14*x.^3 + 60*x.^2 - 70*x;
[x6, ~, o6] = bisection_min(f6, 0, 2, struct('tol', 1e-5));
fprintf('6.3 f = x^4-14x^3+60x^2-70x tren [0,2]: x* = %.6f, %d buoc, %d lan tinh f\n', x6, o6.iterations, o6.nfev);
c = 1;  d = 1e-3;
fprintf('    buoc 1: f(1-d) = %.6f, f(1+d) = %.6f\n', f6(c - d), f6(c + d));

%% Cau 7. Lat cat vang
fprintf('===== Cau 7 =====\n');
tau = (sqrt(5) - 1)/2;
fprintf('7.1 so buoc [0,3] -> 1e-3: %.4f\n', log(1e-3/3)/log(tau));
x1 = 1 - tau;  x2 = tau;
fprintf('7.2 x1 = %.6f f = %.6f; x2 = %.6f f = %.6f; diem moi = %.6f\n', x1, fh(x1), x2, fh(x2), x2 - tau*x2);
[x7, ~, o7] = golden_section(@(x) exp(x) - 2*x, 0, 2, struct('tol', 1e-4));
fprintf('7.3 e^x - 2x tren [0,2]: x* = %.6f (%d buoc), ln 2 = %.6f\n', x7, o7.iterations, log(2));

%% Cau 8. Gauss-Seidel
fprintf('===== Cau 8 =====\n');
H8 = [2 1; 1 2];  q8 = [-3; 0];
x = gs_quad(H8, q8, [0; 0], 1e-6);
fprintf('8.1 nghiem dung %s\n', mat2str((-H8\q8)'));

%% Cau 9. Hooke-Jeeves
fprintf('===== Cau 9 =====\n');
F9 = @(x) (x(1) - 2)^4 + (x(1) - 2*x(2))^2;
fprintf('9.1 F(0,3) = %g, F(0.5,3) = %g, F(0.5,3.5) = %g, F(0.5,2.5) = %g, F(1,2) = %g\n', ...
  F9([0;3]), F9([0.5;3]), F9([0.5;3.5]), F9([0.5;2.5]), F9([1;2]));
[x9, F9v, o9] = hooke_jeeves(F9, [0; 3], struct('delta', 0.5, 'tol', 1e-6));
fprintf('    hooke_jeeves Delta0 = 0.5: x = %s, F = %g, %d lan tinh F\n', mat2str(x9', 6), F9v, o9.nfev);

%% Cau 10. Gradient
fprintf('===== Cau 10 =====\n');
H10 = diag([2 8]);  x0 = [1; 1];  g0 = H10*x0;
al = (g0'*g0)/(g0'*H10*g0);
fprintf('10.1 alpha = %.6f (= 17/130 = %.6f), x1 = %s\n', al, 17/130, mat2str((x0 - al*g0)', 6));
L = 8; mu = 2;
fprintf('10.2 2/L = %.4f, 2/(L+mu) = %.4f, ti so (k-1)/(k+1) = %.4f\n', 2/L, 2/(L+mu), (L/mu - 1)/(L/mu + 1));
fprintf('10.3 heavy-ball alpha = %.6f, beta = %.6f, ti so = %.4f\n', 4/(sqrt(L)+sqrt(mu))^2, ((sqrt(L)-sqrt(mu))/(sqrt(L)+sqrt(mu)))^2, (sqrt(L/mu)-1)/(sqrt(L/mu)+1));
[~, ~, oa] = gradient_descent(@(x) 0.5*x'*H10*x, @(x) H10*x, x0, struct('step', 'fixed', 'alpha', 0.2, 'tol', 1e-6));
[~, ~, ob] = momentum_gradient(@(x) 0.5*x'*H10*x, @(x) H10*x, x0, struct('alpha', 4/(sqrt(L)+sqrt(mu))^2, 'beta', 1/9, 'tol', 1e-6));
fprintf('     so buoc den ||g|| <= 1e-6: gradient 0.2: %d, heavy-ball: %d\n', oa.iterations, ob.iterations);

%% Cau 11. Newton
fprintf('===== Cau 11 =====\n');
x = 1;  for k = 1:3, x = x - 4*x^3/(12*x^2); fprintf('11.1 f = x^4: x_%d = %.6f\n', k, x); end
f11 = @(x) x(1)^4 + x(1)*x(2) + (1 + x(2))^2;
g11 = @(x) [4*x(1)^3 + x(2); x(1) + 2*(1 + x(2))];
h11 = @(x) [12*x(1)^2, 1; 1, 2];
x0 = [0; 0];  d = -h11(x0)\g11(x0);
fprintf('11.2 H(0) eig = %s, d = %s, g''d = %g, f(x0) = %g, f(x0+d) = %g\n', mat2str(eig(h11(x0))', 4), mat2str(d'), g11(x0)'*d, f11(x0), f11(x0 + d));
[x11, f11v, o11] = newton_method(f11, g11, h11, x0, struct('variant', 'modified', 'tol', 1e-10));
fprintf('     Newton cai bien: x* = %s, f* = %.6f, %d buoc\n', mat2str(x11', 6), f11v, o11.iterations);
x11r = fminsearch(f11, [0; 0], optimset('TolX', 1e-10, 'TolFun', 1e-12));
fprintf('     fminsearch: %s\n', mat2str(x11r', 6));
H3 = [4 1; 1 2];  q3 = [-1; 2];  x0 = [5; -7];
fprintf('11.3 bac hai: 1 buoc Newton tu (5,-7): %s, -H\\q = %s\n', mat2str((x0 - H3\(H3*x0 + q3))', 6), mat2str((-H3\q3)', 6));

%% Cau 12. Tua Newton
fprintf('===== Cau 12 =====\n');
s = [1; 0];  y = [2; 1];  r = 1/(y'*s);  I = eye(2);
Hb = (I - r*(s*y'))*(I - r*(y*s')) + r*(s*s');
Hd = I - (y*y')/(y'*y) + (s*s')/(y'*s);
fprintf('12.1 BFGS H1 = %s, H1*y = %s\n', mat2str(Hb, 6), mat2str((Hb*y)'));
fprintf('     DFP   H1 = %s, H1*y = %s\n', mat2str(Hd, 6), mat2str((Hd*y)'));
A12 = [2 0 0; 0 4 1; 0 1 3];  b12 = [1; 1; 1];
[x12, ~, o12] = quasi_newton(@(x) 0.5*x'*A12*x - b12'*x, @(x) A12*x - b12, [0; 0; 0], struct('lineSearch', 'exact', 'tol', 1e-10));
fprintf('12.2 bac hai n = 3, tim kiem chinh xac: %d buoc, ||H - A^-1|| = %.2e\n', o12.iterations, norm(o12.H - inv(A12)));
