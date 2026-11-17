function [x, fval, out] = gradient_projection(prob, x0, opts)
% GRADIENT_PROJECTION  Giai bai toan quy hoach loi
%
%       min f(x)   s.t.  g_i(x) <= 0 (i = 1..m),   h_j(x) = 0 (j = 1..p)
%
%   bang phuong phap chieu gradient cua Rosen (1960, 1961) va bien the metric
%   bien thien cua Murtagh & Sargent (1969, 1973). Tai diem chap nhan duoc x_k,
%   voi tap lam viec W (rang buoc dang thuc + bat dang thuc tich cuc), N la ma
%   tran Jacobi cua cac rang buoc trong W (moi dong la mot gradient):
%
%       P = H - H N'(N H N')^{-1} N H          (H = I: Rosen; H ~ nghich dao Hessian: M-S)
%       d = -P grad f(x),   u = -(N H N')^{-1} N H grad f(x)   (uoc luong nhan tu)
%
%   - d ~= 0: di theo d (nam trong khong gian tiep xuc cua W), sau do keo ve
%     mat W bang buoc Newton (phuc hoi tinh chap nhan duoc) neu rang buoc cong.
%   - d = 0 va u_i >= 0 voi moi bat dang thuc trong W: x la diem KKT, dung.
%   - d = 0 va co u_i < 0: bo rang buoc do khoi W va lap lai.
%
%   [x, fval, out] = gradient_projection(prob, x0)
%   [x, fval, out] = gradient_projection(prob, x0, opts)
%
% Dau vao:
%   prob - struct: .f, .grad (gradient n x 1),
%                  .g, .jacg (m x 1, m x n), .h, .jach (p x 1, p x n) - co the bo trong
%   x0   - diem xuat phat CHAP NHAN DUOC
%   opts - (tuy chon) struct:
%       .metric     : 'identity' (Rosen, mac dinh) | 'bfgs' (Murtagh-Sargent)
%       .lineSearch : 'exact'  - cuc tieu f(x + a d) tren [0, a_max] (chia doi theo dau
%                                dao ham, bisection_min.m cau 6); chi dung khi moi rang
%                                buoc trong W la affine (duong di thang)
%                     'armijo' - quay lui, co phuc hoi (mac dinh)
%       .tol        : dung khi ||d|| <= tol va u >= -tol (mac dinh 1e-8)
%       .tolActive  : g_i(x) >= -tolActive coi la tich cuc (mac dinh 1e-8)
%       .maxIter    : so buoc toi da (mac dinh 500)
%       .display    : in tung buoc (mac dinh false)
%
% Dau ra:
%   x, fval - nghiem gan dung va gia tri f
%   out     - struct: .iterations, .X, .F, .W (tap lam viec cuoi), .lambda (nhan tu
%             cua g, m x 1), .mu (nhan tu cua h), .nfev, .ngev, .drops (so lan bo
%             rang buoc), .exitflag (1: KKT; 0: het maxIter; -1: khong tim duoc buoc;
%             2: buoc qua nho so voi sai so lam tron - dung o gioi han do chinh xac)

  if nargin < 3, opts = struct(); end
  metric  = lower(get_opt(opts, 'metric', 'identity'));
  lsearch = lower(get_opt(opts, 'lineSearch', 'armijo'));
  tol     = get_opt(opts, 'tol', 1e-8);
  tolA    = get_opt(opts, 'tolActive', 1e-8);
  maxIter = get_opt(opts, 'maxIter', 500);
  display = get_opt(opts, 'display', false);

  x = x0(:);
  n = numel(x);
  [gv, ~] = eval_c(prob, 'g', 'jacg', x, n);
  [hv, ~] = eval_c(prob, 'h', 'jach', x, n);
  m = numel(gv);  p = numel(hv);
  if any(gv > 1e-8) || any(abs(hv) > 1e-8)
    error('gradient_projection:x0', 'x0 phai chap nhan duoc.');
  end
  H = eye(n);
  fx = prob.f(x);  gf = prob.grad(x);
  out = struct('iterations', 0, 'X', x, 'F', fx, 'W', [], 'lambda', zeros(m, 1), ...
               'mu', zeros(p, 1), 'nfev', 1, 'ngev', 1, 'drops', 0, 'exitflag', 0);
  if display
    fprintf('%4s %-24s %14s %-8s %-10s %-22s %10s\n', 'k', 'x_k', 'f(x_k)', 'bo', 'W', 'u (nhan tu)', 'alpha');
  end

  W = find(gv >= -tolA)';                      % tap lam viec ban dau
  for k = 1:maxIter
    % ----- huong chieu, bo rang buoc neu can -----
    dropped = [];
    while true
      [gv, Jg] = eval_c(prob, 'g', 'jacg', x, n);
      [~, Jh] = eval_c(prob, 'h', 'jach', x, n);
      N = [Jh; Jg(W, :)];
      if isempty(N)
        u = zeros(0, 1);  d = -H * gf;
      else
        M = N * H * N';
        u = -(M \ (N * H * gf));
        d = -H * gf - H * N' * u;
      end
      uIneq = u(p+1:end);
      if norm(d) > tol * max(1, norm(gf)) || isempty(uIneq) || min(uIneq) >= -tol
        break;
      end
      [~, j] = min(uIneq);                     % bo rang buoc co nhan tu am nhat
      dropped(end+1) = W(j);
      W(j) = [];
      out.drops = out.drops + 1;
    end
    if display
      if isempty(dropped), sd = ''; else, sd = sprintf('bo %s', mat2str(dropped)); end
      fprintf('%4d %-24s %14.8g %-8s %-10s %-22s', k - 1, mat2str(x', 6), fx, sd, mat2str(W), mat2str(u', 4));
    end
    if norm(d) <= tol * max(1, norm(gf))
      out.exitflag = 1;
      if display, fprintf('\n'); end
      break;
    end

    % ----- tim do dai buoc -----
    inact = setdiff(1:m, W);
    if strcmp(lsearch, 'exact')
      [a, xn, nf] = exact_step(prob, x, d, inact, n);
      out.nfev = out.nfev + nf;
    else
      [a, xn, nf] = armijo_step(prob, x, fx, gf, d, W, p, n, get_opt(opts, 'alpha0', 1));
      out.nfev = out.nfev + nf;
    end
    if a == 0
      out.exitflag = -1;
      if display, fprintf('  (khong tim duoc buoc)\n'); end
      break;
    end
    if display, fprintf(' %10.4g\n', a); end

    % ----- cap nhat -----
    gfn = prob.grad(xn);  out.ngev = out.ngev + 1;
    if strcmp(metric, 'bfgs')
      % cap nhat BFGS cho nghich dao Hessian cua ham Lagrange (nhan tu u giu nguyen)
      [~, Jgn] = eval_c(prob, 'g', 'jacg', xn, n);
      [~, Jhn] = eval_c(prob, 'h', 'jach', xn, n);
      Nn = [Jhn; Jgn(W, :)];
      s = xn - x;
      if isempty(N), y = gfn - gf; else, y = (gfn + Nn' * u) - (gf + N' * u); end
      if s' * y > 1e-12 * norm(s) * norm(y)
        rho = 1 / (s' * y);  V = eye(n) - rho * (s * y');
        H = V * H * V' + rho * (s * s');
      end
    end
    fn = prob.f(xn);
    if norm(xn - x) <= 1e-14 * max(1, norm(x)) && abs(fn - fx) <= 4 * eps * max(1, abs(fx))
      out.exitflag = 2;                        % khong con tien duoc: gioi han do chinh xac may
      x = xn;  fx = fn;  gf = gfn;
      break;
    end
    x = xn;  fx = fn;  gf = gfn;
    gv = eval_c(prob, 'g', 'jacg', x, n);
    W = union(W, find(gv >= -tolA)');          % them rang buoc moi tro thanh tich cuc
    W = W(:)';
    out.X(:, end+1) = x;  out.F(end+1) = fx;
    out.iterations = k;
  end

  % nhan tu cuoi cung (theo metric I)
  [~, Jg] = eval_c(prob, 'g', 'jacg', x, n);
  [~, Jh] = eval_c(prob, 'h', 'jach', x, n);
  N = [Jh; Jg(W, :)];
  if ~isempty(N)
    u = -((N * N') \ (N * gf));
    out.mu = u(1:p);
    out.lambda(W) = u(p+1:end);
  end
  out.W = W;
  fval = fx;
end

% =====================================================================
function [a, xn, nf] = exact_step(prob, x, d, inact, n)
% Buoc chinh xac tren duong thang x + a d (rang buoc trong W la affine).
  nf = 0;
  dphi = @(t) prob.grad(x + t * d)' * d;
  feas = @(t) all(gi(prob, x + t * d, inact, n) <= 1e-12);
  % mo rong den khi gap rang buoc chan hoac dao ham doi dau
  lo = 0;  t = 1;
  while feas(t) && dphi(t) < 0 && t < 1e8
    lo = t;  t = 2 * t;
  end
  if ~feas(t)
    hi = t;                                    % tim bien: a_max thuoc [lo, hi]
    for it = 1:200
      mid = (lo + hi) / 2;
      if feas(mid), lo = mid; else, hi = mid; end
      if hi - lo <= 1e-15 * max(1, hi), break; end
    end
    amax = lo;
  else
    amax = t;
  end
  if dphi(amax) <= 0
    a = amax;
  else
    [a, ~, o] = bisection_min(@(t) prob.f(x + t * d), 0, amax, ...
      struct('method', 'derivative', 'df', dphi, 'tol', 1e-14 * amax, 'maxIter', 200));
    nf = nf + o.nfev;
  end
  xn = x + a * d;
end

function [a, xn, nf] = armijo_step(prob, x, fx, gf, d, W, p, n, a0)
% Buoc Armijo tren "duong di" a -> R(x + a d), R = phuc hoi ve mat {g_W = 0, h = 0}:
%   1) neu x + a0 d (sau phuc hoi) vi pham rang buoc ngoai W: tim bien a_max bang
%      chia doi (rang buoc moi tro thanh tich cuc - tuong tu a_max cua Rosen);
%   2) neu a0 tot: thu nhan doi khi f con giam (mo rong);
%   3) quay lui (chia doi) den khi thoa Armijo.
  nf = 0;
  slope = gf' * d;                             % < 0
  a = a0;
  [ok, z] = feasible_at(prob, x, d, a, W, p, n);
  if ~ok
    lo = 0;  hi = a;                           % lo chap nhan duoc, hi khong
    for it = 1:100
      mid = (lo + hi) / 2;
      if feasible_at(prob, x, d, mid, W, p, n), lo = mid; else, hi = mid; end
      if hi - lo <= 1e-14 * hi, break; end
    end
    a = lo;
    [~, z] = feasible_at(prob, x, d, a, W, p, n);
  else
    fz = prob.f(z);  nf = nf + 1;
    for it = 1:30                              % mo rong
      [ok2, z2] = feasible_at(prob, x, d, 2 * a, W, p, n);
      if ~ok2, break; end
      f2 = prob.f(z2);  nf = nf + 1;
      if f2 >= fz, break; end
      a = 2 * a;  z = z2;  fz = f2;
    end
  end
  for it = 1:60                                % quay lui
    fz = prob.f(z);  nf = nf + 1;
    if fz <= fx + 1e-4 * a * slope
      xn = z;
      return;
    end
    a = a / 2;
    [ok, z] = feasible_at(prob, x, d, a, W, p, n);
    if ~ok, break; end
  end
  a = 0;  xn = x;
end

function [ok, z] = feasible_at(prob, x, d, a, W, p, n)
  [z, ok] = restore(prob, x + a * d, W, p, n);
  ok = ok && all(eval_c(prob, 'g', 'jacg', z, n) <= 1e-10);
end

function [z, ok] = restore(prob, z, W, p, n)
% Buoc Newton (chieu theo chuan nho nhat) dua z ve mat g_W(z) = 0, h(z) = 0.
  ok = true;
  for it = 1:30
    [gv, Jg] = eval_c(prob, 'g', 'jacg', z, n);
    [hv, Jh] = eval_c(prob, 'h', 'jach', z, n);
    r = [hv; gv(W)];
    if isempty(r) || norm(r) <= 1e-12
      return;
    end
    N = [Jh; Jg(W, :)];
    z = z - N' * ((N * N') \ r);
  end
  ok = norm(r) <= 1e-9;
end

function v = gi(prob, x, idx, n)
  v = eval_c(prob, 'g', 'jacg', x, n);
  v = v(idx);
end

function [v, J] = eval_c(prob, fname, jname, x, n)
  if isfield(prob, fname) && ~isempty(prob.(fname))
    v = prob.(fname)(x);  v = v(:);
    if nargout > 1
      J = reshape(prob.(jname)(x), numel(v), n);
    end
  else
    v = zeros(0, 1);  J = zeros(0, n);
  end
end

function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
