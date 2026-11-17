function [x, fval, out] = frank_wolfe(f, grad, A, b, Aeq, beq, lb, x0, opts)
% FRANK_WOLFE  Giai bai toan quy hoach loi rang buoc tuyen tinh
%
%       min f(x)   s.t.  A x <= b,  Aeq x = beq,  x >= lb
%
%   (f loi kha vi, mien chap nhan duoc D la da dien BI CHAN) bang giai thuat
%   Frank-Wolfe (gradient co dieu kien). Moi buoc:
%       y_k = argmin { grad f(x_k)' y : y thuoc D }      (quy hoach tuyen tinh)
%       d_k = y_k - x_k,   x_{k+1} = x_k + gamma_k d_k,   gamma_k thuoc [0, 1].
%   Khe Frank-Wolfe  gap_k = grad f(x_k)'(x_k - y_k) >= f(x_k) - f*  dung lam
%   tieu chuan dung (chung chi sai so).
%
%   [x, fval, out] = frank_wolfe(f, grad, A, b, Aeq, beq, lb)
%   [x, fval, out] = frank_wolfe(f, grad, A, b, Aeq, beq, lb, x0, opts)
%
% Dau vao:
%   f, grad          - ham muc tieu va gradient (handle)
%   A, b, Aeq, beq   - rang buoc tuyen tinh (bo trong [] neu khong co)
%   lb               - can duoi cua x (vector, co the -inf; [] = -inf)
%   x0               - diem xuat phat CHAP NHAN DUOC; [] = mot dinh cua D
%                      (nghiem co so cua bai toan tuyen tinh min grad f(0)'y)
%   opts             - (tuy chon) struct:
%       .step    : 'exact'    - gamma_k = argmin f(x_k + gamma d_k) tren [0, gamma_max]
%                               (chia doi theo dau dao ham, bisection_min.m cau 6) (mac dinh)
%                  'standard' - gamma_k = 2/(k + 2)  (k = 0, 1, ...)
%       .variant : 'standard' (mac dinh) hoac 'away' (Frank-Wolfe co buoc lui ra
%                  xa dinh xau nhat, Guelat & Marcotte 1986; x0 nen la mot dinh -
%                  neu khong, x0 duoc dung nhu mot phan tu cua tap S ban dau)
%       .tol     : dung khi gap_k <= tol (mac dinh 1e-6)
%       .maxIter : so buoc toi da (mac dinh 1000)
%       .display : in tung buoc (mac dinh false)
%
% Dau ra:
%   x, fval - nghiem gan dung va gia tri f
%   out     - struct: .iterations, .X (day diem), .F, .gap (khe FW theo buoc),
%             .Y (cac dinh y_k), .gammas, .awaySteps, .dropSteps (bien the 'away'),
%             .nLP (so bai toan tuyen tinh da giai), .exitflag (1: gap <= tol; 0: het maxIter)
%
% Can solve_lp.m (cau 5) va, voi step = 'exact', bisection_min.m (cau 6) trong path.

  if nargin < 8, x0 = []; end
  if nargin < 9, opts = struct(); end
  step    = lower(get_opt(opts, 'step', 'exact'));
  variant = lower(get_opt(opts, 'variant', 'standard'));
  tol     = get_opt(opts, 'tol', 1e-6);
  maxIter = get_opt(opts, 'maxIter', 1000);
  display = get_opt(opts, 'display', false);
  if exist('solve_lp', 'file') ~= 2 || (strcmp(step, 'exact') && exist('bisection_min', 'file') ~= 2)
    error('frank_wolfe:path', ['Can solve_lp.m (cau 5) va bisection_min.m (cau 6): ' ...
      'addpath(''../q05-lagrange-duality'', ''../q06-bisection-search'').']);
  end

  if strcmp(variant, 'away') && ~strcmp(step, 'exact')
    error('frank_wolfe:away', 'Bien the ''away'' can tim kiem chinh xac (opts.step = ''exact'').');
  end
  lmo = @(c) lp_vertex(c, A, b, Aeq, beq, lb);       % "oracle" tuyen tinh
  nLP = 0;
  if isempty(x0)
    x0 = lmo(grad(zeros(numel(lb), 1)));
    nLP = 1;
  end
  x = x0(:);
  n = numel(x);
  % Bien the 'away': x = sum w_j S(:, j), S la tap dinh dang dung
  S = x;  w = 1;

  fx = f(x);
  out = struct('iterations', 0, 'X', x, 'F', fx, 'gap', [], 'Y', zeros(n, 0), ...
               'gammas', [], 'awaySteps', 0, 'dropSteps', 0, 'nLP', 0, 'exitflag', 0);
  % cap phat truoc (them cot bang end+1 trong vong lap ton O(k^2) khi k lon)
  X = zeros(n, maxIter + 1);  X(:, 1) = x;
  Y = zeros(n, maxIter + 1);
  F = zeros(1, maxIter + 1);  F(1) = fx;
  GAP = zeros(1, maxIter + 1);  GAM = zeros(1, maxIter);
  if display
    fprintf('%4s %-26s %-20s %14s %12s %10s %s\n', 'k', 'x_k', 'y_k', 'f(x_k)', 'gap_k', 'gamma_k', '');
  end

  for k = 0:maxIter
    g = grad(x);
    y = lmo(g);  nLP = nLP + 1;
    gap = g' * (x - y);
    GAP(k + 1) = gap;
    Y(:, k + 1) = y;
    if gap <= tol || k == maxIter
      if display
        fprintf('%4d %-26s %-20s %14.8g %12.4e\n', k, mat2str(x', 6), mat2str(y', 4), fx, gap);
      end
      out.exitflag = double(gap <= tol);
      break;
    end

    % ----- chon huong -----
    kind = 'FW';
    d = y - x;  gmax = 1;
    if strcmp(variant, 'away')
      [~, jv] = max(g' * S);                         % dinh "xau nhat" dang dung
      dA = x - S(:, jv);
      if g' * dA < g' * d                            % huong lui giam nhanh hon
        kind = 'A';
        d = dA;
        gmax = w(jv) / (1 - w(jv));
        if ~isfinite(gmax), gmax = 1e12; end         % chi con mot dinh (khong xay ra)
      end
    end

    % ----- do dai buoc -----
    if strcmp(step, 'standard')
      gam = min(2 / (k + 2), gmax);
    else
      dphi = @(t) grad(x + t * d)' * d;
      if dphi(gmax) <= 0
        gam = gmax;
      else
        [gam, ~, ~] = bisection_min(@(t) f(x + t * d), 0, gmax, ...
          struct('method', 'derivative', 'df', dphi, 'tol', 1e-14 * gmax, 'maxIter', 200));
      end
    end

    % ----- cap nhat trong so (bien the 'away') -----
    if strcmp(variant, 'away')
      if strcmp(kind, 'FW')
        w = (1 - gam) * w;
        j = find(sum(abs(S - y), 1) < 1e-9, 1);
        if isempty(j)
          S(:, end+1) = y;  w(end+1) = gam;
        else
          w(j) = w(j) + gam;
        end
      else
        out.awaySteps = out.awaySteps + 1;
        w = (1 + gam) * w;
        w(jv) = w(jv) - gam;
      end
      keep = w > 1e-12;
      if any(~keep), out.dropSteps = out.dropSteps + nnz(~keep); end
      S = S(:, keep);  w = w(keep) / sum(w(keep));
    end

    x = x + gam * d;
    fx = f(x);
    X(:, k + 2) = x;  F(k + 2) = fx;  GAM(k + 1) = gam;
    out.iterations = k + 1;
    if display
      fprintf('%4d %-26s %-20s %14.8g %12.4e %10.4g %s\n', k, mat2str(X(:, k + 1)', 6), ...
        mat2str(y', 4), F(k + 1), gap, gam, kind);
    end
  end
  K = out.iterations;
  out.X = X(:, 1:K + 1);  out.F = F(1:K + 1);  out.gap = GAP(1:K + 1);
  out.Y = Y(:, 1:K + 1);  out.gammas = GAM(1:K);
  out.nLP = nLP;
  fval = fx;
end

% =====================================================================
function y = lp_vertex(c, A, b, Aeq, beq, lb)
% Dinh cua D cuc tieu ham tuyen tinh c'y (solve_lp.m cua cau 5)
  [y, ~, status] = solve_lp(c, A, b, Aeq, beq, lb);
  if status ~= 1
    error('frank_wolfe:lp', 'Bai toan tuyen tinh con khong giai duoc (D rong hoac khong bi chan).');
  end
  y = y(:);
  y(abs(y) < 1e-12) = 0;
end

function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
