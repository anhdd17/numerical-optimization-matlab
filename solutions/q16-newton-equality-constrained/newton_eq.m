function [x, nu, out] = newton_eq(f, grad, hess, A, b, x0, opts)
% NEWTON_EQ  Giai bai toan quy hoach loi voi rang buoc tuyen tinh (dang thuc)
%
%       min f(x)   s.t.  A x = b        (f loi, kha vi hai lan; A: p x n, hang p)
%
%   bang thuat toan Newton (Boyd & Vandenberghe, Convex Optimization, ch. 10).
%   Buoc Newton (dx, w) la nghiem cua he KKT
%
%       [ H   A' ] [ dx ]     [ -g        ]        H = hess f(x), g = grad f(x)
%       [ A   0  ] [ w  ]  =  [ -(A x - b) ]
%
%   - Diem xuat phat CHAP NHAN DUOC (A x0 = b): ve phai thu hai bang 0, moi
%     diem lap deu chap nhan duoc; tim kiem quay lui tren f; dung khi
%     lambda^2/2 <= tol, voi lambda^2 = dx' H dx (giam luong Newton).
%   - Diem xuat phat BAT KY (chi can x0 thuoc mien xac dinh cua f): Newton
%     cho he r(x, nu) = (g + A' nu, A x - b) = 0 (Newton "nguyen-doi ngau");
%     tim kiem quay lui tren ||r||; dung khi ||A x - b|| <= tolFeas va ||r|| <= tol.
%
%   [x, nu, out] = newton_eq(f, grad, hess, A, b, x0)
%   [x, nu, out] = newton_eq(f, grad, hess, A, b, x0, opts)
%
% Dau vao:
%   f, grad, hess - ham, gradient, Hessian (handle); f(x) = Inf ngoai mien xac dinh
%   A, b          - rang buoc A x = b
%   x0            - diem xuat phat (thuoc mien xac dinh cua f)
%   opts          - (tuy chon) struct:
%       .start   : 'auto' (mac dinh: 'feasible' neu ||A x0 - b|| <= tolFeas),
%                  'feasible', 'infeasible'
%       .nu0     : nhan tu ban dau cho 'infeasible' (mac dinh 0)
%       .tol     : sai so (mac dinh 1e-10)
%       .tolFeas : sai so cua A x = b (mac dinh 1e-10 * max(1, ||b||))
%       .alpha, .beta : tham so tim kiem quay lui (mac dinh 0.01, 0.5)
%       .maxIter : so buoc toi da (mac dinh 100)
%       .display : in tung buoc (mac dinh false)
%
% Dau ra:
%   x, nu - nghiem va nhan tu Lagrange cua A x = b (grad f(x) + A' nu = 0)
%   out   - struct: .iterations, .X, .F, .dec (lambda^2/2, ban chap nhan duoc),
%           .rnorm (||r||), .feas (||A x - b||), .ts (cac buoc t), .start,
%           .exitflag (1: dat sai so; 0: het maxIter; -1: he KKT suy bien / x0 ngoai mien)

  if nargin < 7, opts = struct(); end
  tolFeas = get_opt(opts, 'tolFeas', 1e-10 * max(1, norm(b)));
  start   = lower(get_opt(opts, 'start', 'auto'));
  tol     = get_opt(opts, 'tol', 1e-10);
  alpha   = get_opt(opts, 'alpha', 0.01);
  beta    = get_opt(opts, 'beta', 0.5);
  maxIter = get_opt(opts, 'maxIter', 100);
  display = get_opt(opts, 'display', false);

  x = x0(:);
  [p, n] = size(A);
  b = b(:);
  if strcmp(start, 'auto')
    if norm(A * x - b) <= tolFeas, start = 'feasible'; else, start = 'infeasible'; end
  end
  if strcmp(start, 'feasible') && norm(A * x - b) > tolFeas
    error('newton_eq:x0', 'start = ''feasible'' can A*x0 = b (||A x0 - b|| = %.2e).', norm(A * x - b));
  end
  nu = get_opt(opts, 'nu0', zeros(p, 1));
  fx = f(x);
  out = struct('iterations', 0, 'X', x, 'F', fx, 'dec', [], 'rnorm', [], 'feas', norm(A * x - b), ...
               'ts', [], 'start', start, 'exitflag', 0);
  if ~isfinite(fx)
    out.exitflag = -1;
    warning('newton_eq:domain', 'x0 nam ngoai mien xac dinh cua f.');
    return;
  end
  rfun = @(x, nu) [grad(x) + A' * nu; A * x - b];

  if display
    if strcmp(start, 'feasible')
      fprintf('%4s %16s %14s %10s\n', 'k', 'f(x_k)', 'lambda^2/2', 't_k');
    else
      fprintf('%4s %16s %12s %12s %10s\n', 'k', 'f(x_k)', '||r||', '||Ax - b||', 't_k');
    end
  end

  for k = 0:maxIter
    g = grad(x);  H = hess(x);
    K = [H, A'; A, zeros(p)];
    if rcond(K) < 1e-15
      out.exitflag = -1;
      break;
    end
    if strcmp(start, 'feasible')
      sol = K \ [-g; zeros(p, 1)];
      dx = sol(1:n);  w = sol(n+1:end);
      dec = dx' * H * dx / 2;                  % lambda(x)^2 / 2
      out.dec(end+1) = dec;
      nu = w;                                  % uoc luong nhan tu
      if dec <= tol || k == maxIter
        if display, fprintf('%4d %16.10g %14.4e\n', k, fx, dec); end
        out.exitflag = double(dec <= tol);
        break;
      end
      t = 1;                                   % quay lui tren f, giu x trong mien xac dinh
      while true
        ft = f(x + t * dx);
        if isfinite(ft) && ft <= fx + alpha * t * (g' * dx), break; end
        t = beta * t;
        if t < 1e-16, break; end
      end
      if display, fprintf('%4d %16.10g %14.4e %10.4g\n', k, fx, dec, t); end
      x = x + t * dx;  fx = ft;
    else
      r = rfun(x, nu);
      rn = norm(r);
      out.rnorm(end+1) = rn;
      feas = norm(A * x - b);
      if (feas <= tolFeas && rn <= tol) || k == maxIter
        if display, fprintf('%4d %16.10g %12.4e %12.4e\n', k, fx, rn, feas); end
        out.exitflag = double(feas <= tolFeas && rn <= tol);
        break;
      end
      sol = -(K \ r);
      dx = sol(1:n);  dnu = sol(n+1:end);
      t = 1;                                   % quay lui tren ||r||
      while true
        xt = x + t * dx;
        if isfinite(f(xt)) && norm(rfun(xt, nu + t * dnu)) <= (1 - alpha * t) * rn, break; end
        t = beta * t;
        if t < 1e-16, break; end
      end
      if display, fprintf('%4d %16.10g %12.4e %12.4e %10.4g\n', k, fx, rn, feas, t); end
      x = x + t * dx;  nu = nu + t * dnu;  fx = f(x);
    end
    out.ts(end+1) = t;
    out.X(:, end+1) = x;  out.F(end+1) = fx;  out.feas(end+1) = norm(A * x - b);
    out.iterations = k + 1;
  end
end

% =====================================================================
function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
