function res = kkt_check(prob, x, opts)
% KKT_CHECK  Kiem tra cac dieu kien chinh quy (LICQ, MFCQ) va dieu kien KKT
%            tai mot diem x cua bai toan
%
%   min f(x)   s.t.  g_i(x) <= 0 (i = 1..m),   h_j(x) = 0 (j = 1..p)
%
%   res = kkt_check(prob, x)
%   res = kkt_check(prob, x, opts)
%
% Dau vao:
%   prob - struct gom cac ham (handle), bo trong neu khong co rang buoc:
%          .gradf : x -> grad f(x)                (n x 1)
%          .g     : x -> [g_1(x); ...; g_m(x)]    (m x 1)
%          .jacg  : x -> ma tran Jacobi cua g     (m x n), dong i = grad g_i'
%          .h     : x -> [h_1(x); ...; h_p(x)]    (p x 1)
%          .jach  : x -> ma tran Jacobi cua h     (p x n)
%   x    - diem can kiem tra (n x 1)
%   opts - (tuy chon) .tol (mac dinh 1e-6), .display (mac dinh true)
%
% Dau ra: struct res gom
%   .feasible   - x thoa man moi rang buoc (sai so tol)
%   .active     - chi so cac rang buoc bat dang thuc tich cuc g_i(x) = 0
%   .LICQ       - cac gradient tich cuc {grad g_i, i tich cuc} U {grad h_j} doc lap tuyen tinh
%   .MFCQ       - grad h_j doc lap tuyen tinh va ton tai d: grad h_j'*d = 0,
%                 grad g_i'*d < 0 voi moi i tich cuc
%   .lambda     - nhan tu cua g (m x 1, bang 0 voi rang buoc khong tich cuc)
%   .mu         - nhan tu cua h (p x 1)
%   .stationarity - || grad f + sum lambda_i grad g_i + sum mu_j grad h_j ||
%   .KKT        - x la diem KKT (chap nhan duoc va ton tai nhan tu)
%   .uniqueMultipliers - nhan tu (neu co) la duy nhat

  if nargin < 3, opts = struct(); end
  tol     = get_opt(opts, 'tol', 1e-6);
  display = get_opt(opts, 'display', true);

  x  = x(:);
  n  = numel(x);
  gf = prob.gradf(x);
  gf = gf(:);
  [gv, Jg] = eval_constr(prob, 'g', 'jacg', x, n);
  [hv, Jh] = eval_constr(prob, 'h', 'jach', x, n);
  m = numel(gv);
  p = numel(hv);

  % ---------- Tinh chap nhan duoc va tap tich cuc ----------
  res.feasible = all(gv <= tol) && all(abs(hv) <= tol);
  act = find(abs(gv) <= tol)';
  res.active = act;
  JA = Jg(act, :);                          % gradient cac rang buoc tich cuc
  K  = [JA; Jh]';                           % n x (|A| + p)

  % ---------- LICQ ----------
  res.LICQ = isempty(K) || rank(K, tol) == size(K, 2);

  % ---------- MFCQ: max t  s.t.  JA*d + t <= 0, Jh*d = 0, -1 <= d <= 1, 0 <= t <= 1 ----------
  hIndep = isempty(Jh) || rank(Jh, tol) == p;
  if isempty(act)
    res.MFCQ = hIndep;
  elseif ~hIndep
    res.MFCQ = false;
  else
    na = numel(act);
    cLP = [zeros(n, 1); -1];
    Ain = [JA, ones(na, 1)];
    Aeq = [Jh, zeros(p, 1)];
    lb  = [-ones(n, 1); 0];
    ub  = [ ones(n, 1); 1];
    z = solve_lp(cLP, Ain, zeros(na, 1), Aeq, zeros(p, 1), lb, ub);
    res.MFCQ = z(end) > tol;
  end

  % ---------- Nhan tu: min || gf + JA'*lam + Jh'*mu ||, lam >= 0 ----------
  % Tach mu = mu_plus - mu_minus de dung lsqnonneg (co trong MATLAB va Octave).
  C = [JA', Jh', -Jh'];
  if isempty(C)
    v = zeros(0, 1);
  else
    v = lsqnonneg(C, -gf);
  end
  na = numel(act);
  lamA = reshape(v(1:na), na, 1);
  muv  = reshape(v(na+1:na+p) - v(na+p+1:end), p, 1);
  res.lambda = zeros(m, 1);
  res.lambda(act) = lamA;
  res.mu = muv;
  r = gf + JA' * lamA + Jh' * muv;
  res.stationarity = norm(r);
  res.KKT = res.feasible && res.stationarity <= tol * max(1, norm(gf));
  res.uniqueMultipliers = res.KKT && (isempty(K) || rank(K, tol) == size(K, 2));

  if display
    print_report(res, x, gv, hv, m, p);
  end
end

% =====================================================================
function [v, J] = eval_constr(prob, fname, jname, x, n)
  if isfield(prob, fname) && ~isempty(prob.(fname))
    v = prob.(fname)(x);
    v = v(:);
    J = prob.(jname)(x);
    J = reshape(J, numel(v), n);
  else
    v = zeros(0, 1);
    J = zeros(0, n);
  end
end

function z = solve_lp(c, Ain, bin, Aeq, beq, lb, ub)
% min c'*z s.t. Ain*z <= bin, Aeq*z = beq, lb <= z <= ub (glpk hoac linprog).
  if exist('OCTAVE_VERSION', 'builtin') ~= 0
    A = [Ain; Aeq];
    b = [bin; beq];
    ctype = [repmat('U', 1, size(Ain, 1)), repmat('S', 1, size(Aeq, 1))];
    z = glpk(c, A, b, lb, ub, ctype, repmat('C', 1, numel(c)), 1);
  else
    z = linprog(c, Ain, bin, Aeq, beq, lb, ub, optimoptions('linprog', 'Display', 'none'));
  end
end

function print_report(res, x, gv, hv, m, p)
  yn = {'khong', 'co'};
  fprintf('  Diem x = %s\n', mat2str(x', 5));
  if m > 0, fprintf('  g(x) = %s\n', mat2str(gv' + 0, 5)); end
  if p > 0, fprintf('  h(x) = %s\n', mat2str(hv' + 0, 5)); end
  fprintf('  Chap nhan duoc: %s;  rang buoc bat dang thuc tich cuc: %s\n', ...
    yn{res.feasible + 1}, mat2str(res.active));
  fprintf('  LICQ: %s;  MFCQ: %s\n', yn{res.LICQ + 1}, yn{res.MFCQ + 1});
  if m > 0, fprintf('  lambda = %s\n', mat2str(res.lambda', 5)); end
  if p > 0, fprintf('  mu     = %s\n', mat2str(res.mu', 5)); end
  fprintf('  ||grad_x L|| = %.2e  =>  diem KKT: %s', res.stationarity, yn{res.KKT + 1});
  if res.KKT && ~res.uniqueMultipliers
    fprintf(' (nhan tu KHONG duy nhat)');
  end
  fprintf('\n\n');
end

function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
