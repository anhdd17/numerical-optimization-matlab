function [x, fval, out] = barrier_method(prob, x0, opts)
% BARRIER_METHOD  Giai bai toan quy hoach loi tong quat
%
%       min f0(x)   s.t.  f_i(x) <= 0 (i = 1..m, f_i loi),   A x = b
%
%   bang thuat toan ham chan Newton (log barrier; Boyd & Vandenberghe, ch. 11).
%   Vong ngoai: voi t tang dan (t <- mu*t), giai bai toan dinh tam
%
%       min  t f0(x) + phi(x),   phi(x) = -sum_i log(-f_i(x)),   s.t.  A x = b
%
%   bang Newton co rang buoc dang thuc (newton_eq.m, cau 16), xuat phat tu
%   nghiem cua buoc truoc. Nghiem x*(t) nam tren "duong trung tam"; diem
%   lambda_i = -1/(t f_i(x*(t))) chap nhan duoc cho bai toan doi ngau va
%   khe doi ngau bang m/t, nen dung khi m/t <= tol.
%   Neu x0 khong chap nhan chat (f_i(x0) < 0) hoac A x0 ~= b, giai truoc
%   bai toan pha 1  min s  s.t. f_i(x) <= s, A x = b, s >= -1  cung bang ham chan.
%
%   [x, fval, out] = barrier_method(prob, x0)
%   [x, fval, out] = barrier_method(prob, x0, opts)
%
% Dau vao:
%   prob - struct:
%       .f0, .grad0, .hess0 : ham muc tieu, gradient, Hessian
%       .fi    : x -> [f_1(x); ...; f_m(x)]          (m x 1)
%       .jac   : x -> ma tran Jacobi cua fi          (m x n)
%       .hessc : (x, w) -> sum_i w_i * hess f_i(x)   (n x n); bo trong neu moi f_i affine
%       .A, .b : rang buoc dang thuc (co the bo trong)
%   x0   - diem xuat phat (bat ky; neu khong chap nhan chat thi chay pha 1)
%   opts - (tuy chon) struct:
%       .t0      : t ban dau (mac dinh 1)
%       .mu      : he so tang t (mac dinh 10)
%       .tol     : sai so khe doi ngau m/t (mac dinh 1e-8)
%       .tolInner: sai so lambda^2/2 cua Newton dinh tam (mac dinh 1e-10)
%       .display : in tung buoc ngoai (mac dinh false)
%       .stopNeg : (dung noi bo cho pha 1) dung ngay khi f0(x) < 0
%
% Dau ra:
%   x, fval - nghiem gan dung va gia tri f0
%   out     - struct: .outer (so buoc ngoai), .newton (so buoc Newton moi buoc
%             ngoai), .t, .gap (m/t), .X (cac diem trung tam), .F, .lambda (nhan tu
%             cua f_i), .nu (nhan tu cua A x = b), .phase1 (ket qua pha 1, neu co),
%             .exitflag (1: dat tol; -2: bai toan khong co diem chap nhan chat)
%
% Can newton_eq.m (cau 16) trong path.

  if nargin < 3, opts = struct(); end
  t       = get_opt(opts, 't0', 1);
  mu      = get_opt(opts, 'mu', 10);
  tol     = get_opt(opts, 'tol', 1e-8);
  tolIn   = get_opt(opts, 'tolInner', 1e-10);
  display = get_opt(opts, 'display', false);
  stopNeg = get_opt(opts, 'stopNeg', false);
  if exist('newton_eq', 'file') ~= 2
    error('barrier_method:path', 'Can newton_eq.m (cau 16): addpath(''../q16-newton-equality-constrained'').');
  end

  x = x0(:);
  n = numel(x);
  if isfield(prob, 'A') && ~isempty(prob.A)
    A = prob.A;  b = prob.b(:);
  else
    A = zeros(0, n);  b = zeros(0, 1);
  end
  if ~isfield(prob, 'hessc') || isempty(prob.hessc)
    prob.hessc = @(x, w) zeros(n);
  end
  m = numel(prob.fi(x));
  out = struct('outer', 0, 'newton', [], 't', [], 'gap', [], 'X', zeros(n, 0), 'F', [], ...
               'lambda', zeros(m, 1), 'nu', zeros(size(A, 1), 1), 'phase1', [], 'exitflag', 0);

  % ----- pha 1: tim diem chap nhan chat -----
  if any(prob.fi(x) >= 0) || norm(A * x - b) > 1e-9 * max(1, norm(b))
    s0 = max(max(prob.fi(x)), 0) + 1;
    P1.f0 = @(z) z(end);
    P1.grad0 = @(z) [zeros(n, 1); 1];
    P1.hess0 = @(z) zeros(n + 1);
    P1.fi = @(z) [prob.fi(z(1:n)) - z(end); -1 - z(end)];
    P1.jac = @(z) [prob.jac(z(1:n)), -ones(m, 1); zeros(1, n), -1];
    P1.hessc = @(z, w) blkdiag(prob.hessc(z(1:n), w(1:m)), 0);
    P1.A = [A, zeros(size(A, 1), 1)];  P1.b = b;
    if norm(A * x - b) > 1e-9 * max(1, norm(b))
      x = x - A' * ((A * A') \ (A * x - b));     % chieu len {A x = b}; s du lon de chap nhan chat
      s0 = max(max(prob.fi(x)), 0) + 1;
    end
    [z, s, o1] = barrier_method(P1, [x; s0], struct('mu', mu, 'tol', 1e-8, 'stopNeg', true));
    out.phase1 = struct('s', s, 'x', z(1:n), 'outer', o1.outer, 'newton', sum(o1.newton));
    if display
      fprintf('Pha 1: s* = %.4g sau %d buoc ngoai, %d buoc Newton\n', s, o1.outer, sum(o1.newton));
    end
    if s >= 0
      out.exitflag = -2;                         % khong co diem chap nhan chat
      x = z(1:n);  fval = prob.f0(x);
      return;
    end
    x = z(1:n);
  end

  % ----- vong ngoai -----
  if display
    fprintf('%4s %10s %16s %12s %8s\n', 'k', 't', 'f0(x*(t))', 'm/t', 'Newton');
  end
  for k = 1:200
    Ft = @(z) t * prob.f0(z) - sum(log(max(-prob.fi(z), 0)));       % = +Inf ngoai mien
    Gt = @(z) t * prob.grad0(z) + prob.jac(z)' * (1 ./ -prob.fi(z));
    Ht = @(z) hess_t(prob, z, t);
    [x, nu, o] = newton_eq(Ft, Gt, Ht, A, b, x, struct('start', 'feasible', 'tol', tolIn, 'maxIter', 200));
    if o.exitflag ~= 1
      warning('barrier_method:centering', 'Buoc dinh tam t = %g ket thuc voi exitflag = %d.', t, o.exitflag);
    end
    out.outer = k;
    out.newton(end+1) = o.iterations;
    out.t(end+1) = t;  out.gap(end+1) = m / t;
    out.X(:, end+1) = x;  out.F(end+1) = prob.f0(x);
    out.lambda = -1 ./ (t * prob.fi(x));
    out.nu = nu / t;
    if display
      fprintf('%4d %10.4g %16.10g %12.4e %8d\n', k, t, out.F(end), m / t, o.iterations);
    end
    if stopNeg && out.F(end) < 0
      break;
    end
    if m / t <= tol
      out.exitflag = 1;
      break;
    end
    t = mu * t;
  end
  fval = prob.f0(x);
end

% =====================================================================
function H = hess_t(prob, x, t)
% Hessian cua t f0 + phi
  fv = prob.fi(x);  J = prob.jac(x);
  d = 1 ./ (-fv);
  H = t * prob.hess0(x) + J' * diag(d.^2) * J + prob.hessc(x, d);
end

function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
