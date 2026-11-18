% CHECK_KKT_DUALITY  Kiem chung so lieu cua topics/kkt-duality-convex.md
% (cau 4, 5, 13: chinh quy, KKT, doi ngau, quy hoach loi).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> check_kkt_duality
exam_setup;
o = struct('display', false);

%% ================= CAU 4 =================
fprintf('############ CAU 4 ############\n');
% Q4.1: min (x1-3)^2 + (x2-3)^2,  g1 = x1 + x2 - 4 <= 0,  g2 = -x1 <= 0
p.gradf = @(x) [2*(x(1) - 3); 2*(x(2) - 3)];
p.g = @(x) [x(1) + x(2) - 4; -x(1)];
p.jacg = @(x) [1 1; -1 0];
r = kkt_check(p, [2; 2], o);
fprintf('Q4.1 tai (2,2): KKT = %d, lambda = %s, LICQ = %d\n', r.KKT, mat2str(r.lambda', 4), r.LICQ);
r = kkt_check(p, [0; 4], o);
lam = [1 -1; 1 0] \ (-p.gradf([0; 4]));   % giai he dung voi g1, g2 tich cuc
fprintf('Q4.1 tai (0,4): KKT = %d, LICQ = %d, nghiem he dung lambda = %s\n', r.KKT, r.LICQ, mat2str(lam', 4));
q.f = @(x) (x(1) - 3)^2 + (x(2) - 3)^2;  q.g = p.g;
[xs, fs] = solve_nlp(q, [0; 0]);
fprintf('Q4.1 solver: x = %s, f = %.6f\n', mat2str(xs', 6), fs);

% Q4.2: g1 = -x1, g2 = -x2, g3 = -x1 - 2x2 tai (0,0), f = x1 + x2
p2.gradf = @(x) [1; 1];
p2.g = @(x) [-x(1); -x(2); -x(1) - 2*x(2)];
p2.jacg = @(x) [-1 0; 0 -1; -1 -2];
r = kkt_check(p2, [0; 0], o);
fprintf('Q4.2: LICQ = %d, MFCQ = %d, KKT = %d, nhan tu duy nhat = %d\n', r.LICQ, r.MFCQ, r.KKT, r.uniqueMultipliers);

% Q4.3: min x1^2 + x2^2 + x3^2,  x1 + 2x2 + 3x3 = 14
p3.gradf = @(x) 2*x;  p3.h = @(x) x(1) + 2*x(2) + 3*x(3) - 14;  p3.jach = @(x) [1 2 3];
r = kkt_check(p3, [1; 2; 3], o);
fprintf('Q4.3: KKT = %d, mu = %.6f, f = %g\n', r.KKT, r.mu, 14);

% E4.1: min (x1-2)^2 + (x2-3)^2,  g1 = x1 + x2 - 3 <= 0,  g2 = x1 - 2x2 <= 0
e.gradf = @(x) [2*(x(1) - 2); 2*(x(2) - 3)];
e.g = @(x) [x(1) + x(2) - 3; x(1) - 2*x(2)];
e.jacg = @(x) [1 1; 1 -2];
r = kkt_check(e, [1; 2], o);
fprintf('E4.1 tai (1,2): KKT = %d, lambda = %s, LICQ = %d\n', r.KKT, mat2str(r.lambda', 4), r.LICQ);
xb = [1 1; 1 -2] \ [3; 0];
lb = [1 1; 1 -2]' \ (-e.gradf(xb));
fprintf('E4.1 ca hai tich cuc: x = %s, lambda = %s\n', mat2str(xb', 4), mat2str(lb', 4));
q.f = @(x) (x(1) - 2)^2 + (x(2) - 3)^2;  q.g = e.g;
[xs, fs] = solve_nlp(q, [0; 0]);
fprintf('E4.1 solver: x = %s, f = %.6f\n', mat2str(xs', 6), fs);

%% ================= CAU 5 =================
fprintf('############ CAU 5 ############\n');
% Q5.1: min x1^2 + x2^2,  g = 2 - x1 - x2 <= 0
res = dual_qp(2*eye(2), [0; 0], [-1 -1], -2);
fprintf('Q5.1: q(1) = %.4f, q(2) = %.4f, lambda* = %.4f, d* = %.4f, p* = %.4f, x* = %s\n', ...
  res.q(1, []), res.q(2, []), res.lambda, res.d, res.p, mat2str(res.xP', 4));

% Q5.2: min 2x1 + 3x2,  x1 + x2 >= 4,  x1 + 3x2 >= 6,  x >= 0  (them bien bu)
r = lp_dual([2 3 0 0]', [1 1 -1 0; 1 3 0 -1], [4; 6], false);
fprintf('Q5.2: x* = %s, p* = %g; y* = %s, d* = %g; bu = %s\n', mat2str(r.x', 4), r.p, mat2str(r.y', 4), r.d, mat2str(r.compl', 3));

% Q5.3: min x s.t. x^2 <= 0: q(lambda) = -1/(4 lambda)
fprintf('Q5.3: q(1) = %.4f, q(10) = %.4f, q(1000) = %.6f (-> 0 = p*, khong dat)\n', -1/4, -1/40, -1/4000);

% E5.1: min x1^2 + x2^2 - 2x1 - 4x2,  x1 + x2 <= 1
res = dual_qp(2*eye(2), [-2; -4], [1 1], 1);
fprintf('E5.1: lambda* = %.4f, d* = %.4f, p* = %.4f, x(lambda*) = %s, q(0) = %.4f, q(1) = %.4f\n', ...
  res.lambda, res.d, res.p, mat2str(res.xD', 4), res.q(0, []), res.q(1, []));

% E5.2: max 3x1 + 2x2,  x1 + x2 <= 4,  x1 + 3x2 <= 6,  x1 <= 3.5,  x >= 0
%       = min -3x1 - 2x2 voi bien bu x3, x4, x5
r = lp_dual([-3 -2 0 0 0]', [1 1 1 0 0; 1 3 0 1 0; 1 0 0 0 1], [4; 6; 3.5], false);
fprintf('E5.2: x* = %s, z* = %g; y (bai max) = %s, w* = %g\n', mat2str(r.x(1:2)', 4), -r.p, mat2str(-r.y', 4), -r.d);
[xl, fl] = solve_lp([4; 6; 3.5], -[1 1; 1 3; 1 0]', -[3; 2], [], [], zeros(3, 1));
fprintf('E5.2: giai truc tiep bai doi ngau min 4y1+6y2+3.5y3: y = %s, w = %g\n', mat2str(xl', 4), fl);

%% ================= CAU 13 =================
fprintf('############ CAU 13 ############\n');
% Q13.1: g1 = x1^2 + x2^2 - 1 <= 0, g2 = 1 - x1 <= 0 (affine); f = x2
s.g = @(x) [x(1)^2 + x(2)^2 - 1; 1 - x(1)];
[ok1, xh1, s1] = slater_check(s, [0; 0]);
[ok2, xh2, s2] = slater_check(s, [0; 0], [false; true]);
fprintf('Q13.1: Slater chat: ok = %d (s* = %.2e); Slater tinh: ok = %d (s* = %.2e)\n', ok1, s1, ok2, s2);
p4.gradf = @(x) [0; 1];  p4.g = s.g;  p4.jacg = @(x) [2*x(1) 2*x(2); -1 0];
r = kkt_check(p4, [1; 0], o);
fprintf('Q13.1 tai (1,0): KKT = %d, LICQ = %d, ||grad L|| (nhan tu khong am gan nhat) = %.4f\n', r.KKT, r.LICQ, r.stationarity);

% Q13.2: water filling alpha = (1, 2, 4), P = 3
[x, nu, lambda, L] = water_filling([1; 2; 4], 3);
fprintf('Q13.2: x = %s, nu = %.6f, lambda = %s, L = %.4f\n', mat2str(x', 4), nu, mat2str(lambda', 4), L);

% Q13.3: min x1^2 + x2^2,  g = 1 - x1 - x2 <= 0, diem (0.5, 0.5)
p5.gradf = @(x) 2*x;  p5.g = @(x) 1 - x(1) - x(2);  p5.jacg = @(x) [-1 -1];
r = kkt_check(p5, [0.5; 0.5], o);
fprintf('Q13.3: KKT = %d, lambda = %.4f\n', r.KKT, r.lambda);

% E13.1: min -sum log(alpha_i + x_i), sum x = 1, x >= 0, alpha = (1, 1.5, 3)
al = [1; 1.5; 3];
[x, nu, lambda, L] = water_filling(al, 1);
f13 = -sum(log(al + x));
mu = nu - lambda;
qd = -nu * 1 + sum(log(mu) + 1 - al .* mu);
fprintf('E13.1: x = %s, nu = %.6f (= 4/7 = %.6f), lambda = %s (5/21 = %.6f), L = %.4f\n', ...
  mat2str(x', 4), nu, 4/7, mat2str(lambda', 6), 5/21, L);
fprintf('E13.1: f* = %.6f, q(nu*, lambda*) = %.6f\n', f13, qd);
w.f = @(x) -sum(log(al + x));  w.g = @(x) -x;  w.h = @(x) sum(x) - 1;
[xs, fs] = solve_nlp(w, [1; 1; 1] / 3);
fprintf('E13.1 solver: x = %s, f = %.6f\n', mat2str(xs', 6), fs);

% E13.3: min sum (1/2 d_i x_i^2 - c_i x_i), sum x = 1, x >= 0; d = (1,1,1), c = (2, 1.5, 0)
c = [2; 1.5; 0];
v.f = @(x) 0.5 * sum(x.^2) - c' * x;  v.g = @(x) -x;  v.h = @(x) sum(x) - 1;
[xs, fs] = solve_nlp(v, [1; 1; 1] / 3);
xt = [0.75; 0.25; 0];  nut = -1.25;  lt = xt - c - nut;
fprintf('E13.3: ly thuyet x = %s, nu = %.4f, lambda = %s, f = %.6f; solver x = %s, f = %.6f\n', ...
  mat2str(xt', 4), nut, mat2str(lt', 4), v.f(xt), mat2str(xs', 6), fs);
