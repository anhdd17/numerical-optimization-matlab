% CHECK_CONSTRAINED  Tinh dap an cho cac cau hoi mau cua topics/constrained-methods.md
% (cau 14 Frank-Wolfe, cau 15 Rosen / Murtagh-Sargent, cau 16 Newton co rang buoc
% dang thuc, cau 17 ham chan). Dung lai chuong trinh cua solutions/q14..q17.
%   >> cd exams/study-guide
%   >> check_constrained
exam_setup;

%% ===================== CAU 14: FRANK-WOLFE =====================
fprintf('\n############ CAU 14 ############\n');
% Mau 14.1 + 14.2: min (x1-2.5)^2 + 4(x2-1.25)^2, x1 + 2x2 <= 4, x >= 0, x0 = (1, 0)
f1 = @(x) (x(1) - 2.5)^2 + 4*(x(2) - 1.25)^2;
g1 = @(x) [2*(x(1) - 2.5); 8*(x(2) - 1.25)];
A1 = [1 2];  b1 = 4;  lb1 = [0; 0];
[~, ~, o] = frank_wolfe(f1, g1, A1, b1, [], [], lb1, [1; 0], struct('display', true, 'maxIter', 3));
fprintf('y0 = %s, G0 = %g, gamma0 = %g, x1 = %s\n', mat2str(o.Y(:, 1)'), o.gap(1), o.gammas(1), mat2str(o.X(:, 2)'));
[xq, fq] = solve_qp(diag([2 8]), [-5; -10], A1, b1, [], [], lb1);
fprintf('solve_qp: x* = %s, f* = %.6f (hang so 2.5^2 + 4*1.25^2 = %g)\n', mat2str(xq', 6), fq + 2.5^2 + 4*1.25^2, 2.5^2 + 4*1.25^2);
[xa, fa, oa] = frank_wolfe(f1, g1, A1, b1, [], [], lb1, [1; 0], struct('variant', 'away', 'tol', 1e-8));
fprintf('FW buoc lui: x = %s, f = %.10f, %d buoc\n', mat2str(xa', 8), fa, oa.iterations);
% Mau 14.3: min (x1-1)^2 + (x2-2)^2, x1 + x2 <= 2, x >= 0, x0 = (0, 0)
f3 = @(x) (x(1) - 1)^2 + (x(2) - 2)^2;
g3 = @(x) [2*(x(1) - 1); 2*(x(2) - 2)];
[x3, f3v, o3] = frank_wolfe(f3, g3, [1 1], 2, [], [], [0; 0], [0; 0], struct('display', true, 'tol', 1e-12));
fprintf('Mau 14.3: x = %s, f = %.10f, so buoc = %d, gap = %s\n', mat2str(x3', 8), f3v, o3.iterations, mat2str(o3.gap, 4));
fprintf('gamma_k = 2/(k+2): k = 0..4 -> %s\n', mat2str(2 ./ ((0:4) + 2), 4));

%% ===================== CAU 15: ROSEN =====================
fprintf('\n############ CAU 15 ############\n');
P1.f    = @(x) 2*x(1)^2 + 2*x(2)^2 - 2*x(1)*x(2) - 4*x(1) - 6*x(2);
P1.grad = @(x) [4*x(1) - 2*x(2) - 4; 4*x(2) - 2*x(1) - 6];
P1.g    = @(x) [x(1) + x(2) - 2; x(1) + 5*x(2) - 5; -x(1); -x(2)];
P1.jacg = @(x) [1 1; 1 5; -1 0; 0 -1];
[xr, fr, orr] = gradient_projection(P1, [0; 0], struct('lineSearch', 'exact', 'display', true));
fprintf('Rosen: x* = %s, f* = %.8f, lambda = %s, so buoc = %d\n', mat2str(xr', 8), fr, mat2str(orr.lambda', 6), orr.iterations);
fprintf('35/31 = %.8f, 24/31 = %.8f, -222/31 = %.8f, 32/31 = %.8f\n', 35/31, 24/31, -222/31, 32/31);
% Buoc 0 bang tay: W = {3, 4}
N = [-1 0; 0 -1];  gr = P1.grad([0; 0]);
u = -(N*N') \ (N*gr);
fprintf('Buoc 0: grad = %s, u = %s\n', mat2str(gr'), mat2str(u'));
N = [-1 0];  P = eye(2) - N'*((N*N')\N);  d = -P*gr;
fprintf('Bo (4): P = %s, d = %s\n', mat2str(P), mat2str(d'));
% alpha_max: rang buoc (1), (2)
amax = min([(2 - 0) / ([1 1]*d), (5 - 0) / ([1 5]*d)]);
fprintf('alpha_max = %.6f; cuc tieu theo tia tai %.6f\n', amax, fminbnd(@(a) P1.f(a*d), 0, 1));
% Buoc 1: x1 = (0,1), W = {2, 3}
x1 = [0; 1];  gr = P1.grad(x1);  N = [1 5; -1 0];
u = -(N*N') \ (N*gr);
fprintf('Buoc 1: grad = %s, u = %s\n', mat2str(gr'), mat2str(u', 6));
N = [1 5];  P = eye(2) - N'*((N*N')\N);  d = -P*gr;
fprintf('W = {2}: 26*P = %s, 26*d = %s\n', mat2str(26*P, 6), mat2str(26*d', 6));
% Murtagh-Sargent
[xm, fm, om] = gradient_projection(P1, [0; 0], struct('lineSearch', 'exact', 'metric', 'bfgs'));
fprintf('Murtagh-Sargent: x = %s, so buoc = %d\n', mat2str(xm', 8), om.iterations);
% Cong thuc M-S voi H = I trung voi Rosen
H = eye(2);  PH = H - H*N'*((N*H*N')\(N*H));
fprintf('||P_H(H = I) - P|| = %g\n', norm(PH - P));

%% ===================== CAU 16: NEWTON CO RANG BUOC DANG THUC =====================
fprintf('\n############ CAU 16 ############\n');
% Mau 16.1: min 1/2||x||^2, x1 + x2 + x3 = 3, x1 - x3 = 1
fA = @(x) 0.5 * (x' * x);  gA = @(x) x;  hA = @(x) eye(3);
AA = [1 1 1; 1 0 -1];  bA = [3; 1];
[xA, nuA, oA] = newton_eq(fA, gA, hA, AA, bA, [2; 0; 1], struct('display', true));
fprintf('x* = %s, nu = %s, f* = %g, lambda^2/2 tai x0 = %g, f(x0) - f* = %g\n', mat2str(xA', 6), mat2str(nuA', 6), fA(xA), oA.dec(1), fA([2;0;1]) - fA(xA));
[xB, nuB, oB] = newton_eq(fA, gA, hA, AA, bA, [0; 0; 0], struct('display', true));
fprintf('Xuat phat 0: x = %s sau %d buoc, ||Ax0 - b|| = %.4f\n', mat2str(xB', 6), oB.iterations, norm(bA));
% Mau 16.2: min x1^4 + x2^2, x1 + x2 = 2, x0 = (1, 1)
fC = @(x) x(1)^4 + x(2)^2;  gC = @(x) [4*x(1)^3; 2*x(2)];  hC = @(x) diag([12*x(1)^2, 2]);
AC = [1 1];  bC = 2;  x0 = [1; 1];
K = [hC(x0) AC'; AC 0];  sol = K \ [-gC(x0); 0];
dx = sol(1:2);  w = sol(3);  lam2 = dx' * hC(x0) * dx;
fprintf('Buoc Newton tai (1,1): dx = %s, w = %.6f, lambda^2 = %.6f, lambda^2/2 = %.6f\n', mat2str(dx', 6), w, lam2, lam2 / 2);
fprintf('-1/7 = %.6f, 16/7 = %.6f, 2/7 = %.6f\n', -1/7, 16/7, 2/7);
[xC, nuC, oC] = newton_eq(fC, gC, hC, AC, bC, x0, struct('display', true));
fprintf('x* = %s, nu = %.8f, f* = %.8f, %d buoc\n', mat2str(xC', 8), nuC, fC(xC), oC.iterations);
fprintf('Nghiem 2x^3 + x - 2 = 0: %s\n', mat2str(roots([2 0 1 -2])', 8));
% Mau 16.3: phan du nguyen thuy sau buoc t: (1 - t) ||Ax - b||
fprintf('(1 - 0.5) * %.4f = %.4f\n', norm(bA), 0.5 * norm(bA));

%% ===================== CAU 17: HAM CHAN =====================
fprintf('\n############ CAU 17 ############\n');
Q.f0    = P1.f;  Q.grad0 = P1.grad;  Q.hess0 = @(x) [4 -2; -2 4];
Q.fi    = P1.g;  Q.jac   = P1.jacg;
[xb, fb, ob] = barrier_method(Q, [0.5; 0.5], struct('display', true));
fprintf('Ham chan: x = %s, f = %.9f, so buoc ngoai = %d, tong Newton = %d, lambda = %s\n', ...
  mat2str(xb', 8), fb, ob.outer, sum(ob.newton), mat2str(ob.lambda', 5));
fprintf('ceil(log10(4/1e-8)) = %d\n', ceil(log10(4 / 1e-8)));
fprintf('Khe m/t: %s\n', mat2str(ob.gap(1:4), 4));
% Mau 17.2: min x s.t. 1 - x <= 0: x*(t) = 1 + 1/t
for t = [1 10 100]
  xt = fminbnd(@(x) t*x - log(x - 1), 1 + 1e-12, 10, optimset('TolX', 1e-12));
  fprintf('t = %g: x*(t) = %.8f (1 + 1/t = %.8f), lambda = -1/(t f1) = %.6f\n', t, xt, 1 + 1/t, -1 / (t * (1 - xt)));
end
% Mau 17.3: cau 13 vi du 1, xuat phat khong chap nhan duoc
R.f0    = @(x) (x(1) - 3)^2 + (x(2) - 2)^2;
R.grad0 = @(x) [2*(x(1) - 3); 2*(x(2) - 2)];
R.hess0 = @(x) 2*eye(2);
R.fi    = @(x) [x(1)^2 + x(2)^2 - 5; x(1) + 2*x(2) - 4];
R.jac   = @(x) [2*x(1), 2*x(2); 1, 2];
R.hessc = @(x, w) 2*w(1)*eye(2);
[xR, fR, oR] = barrier_method(R, [3; 3]);
fprintf('Cau 13 vi du 1: x = %s, f = %.8f, lambda = %s, pha 1 s = %.4f\n', mat2str(xR', 8), fR, mat2str(oR.lambda', 6), oR.phase1.s);
