function [x, fval, exitflag, out] = big_m_simplex(c, A, b, opts)
% BIG_M_SIMPLEX  Giai bai toan quy hoach tuyen tinh dang chinh tac bang
%                phuong phap ham phat (bai toan M).
%
%   min  f(x) = c'*x
%   s.t. A*x = b,  x >= 0
%
%   Them bien gia x^a cho cac dong chua co cot don vi va giai bai toan M:
%   min  c'*x + M*sum(x^a)   s.t.  A*x + E*x^a = b,  x, x^a >= 0
%
%   [x, fval, exitflag, out] = big_m_simplex(c, A, b)
%   [x, fval, exitflag, out] = big_m_simplex(c, A, b, opts)
%
% Dau vao:
%   c, A, b - du lieu bai toan (A: m x n)
%   opts    - (tuy chon) struct gom cac truong:
%             .M       : [] (mac dinh) - M la so lon tuy y: moi uoc luong
%                        Delta_j = alpha_j*M + beta_j duoc luu thanh hai
%                        thanh phan va so sanh theo thu tu tu dien;
%                        so duong   - dung gia tri M cu the, goi simplex.m
%                        cua cau 1 (can addpath('../q01-simplex'))
%             .rule    : 'dantzig' (mac dinh) hoac 'bland'
%             .display : true (mac dinh) in bang don hinh moi buoc
%             .maxIter : so buoc lap toi da (mac dinh 100)
%             .tol     : sai so so hoc (mac dinh 1e-10)
%
% Dau ra:
%   x        - phuong an toi uu cua bai toan goc (n x 1)
%   fval     - gia tri toi uu c'*x
%   exitflag -  1 : tim duoc phuong an toi uu
%               0 : vuot qua so buoc lap toi da
%              -2 : bai toan goc khong co phuong an chap nhan duoc
%                   (bai toan M toi uu nhung con bien gia duong)
%              -3 : bai toan M khong bi chan duoi => bai toan goc khong co
%                   phuong an toi uu
%              -4 : phat hien xoay vong
%   out      - struct: .artificial (chi so bien gia), .xa (gia tri bien gia),
%              .iterations, .basis, .alpha, .beta (uoc luong cuoi),
%              .tableau, .xHist (day phuong an, chi phan x goc)

  if nargin < 4, opts = struct(); end
  Mval    = get_opt(opts, 'M', []);
  rule    = get_opt(opts, 'rule', 'dantzig');
  display = get_opt(opts, 'display', true);
  maxIter = get_opt(opts, 'maxIter', 100);
  tol     = get_opt(opts, 'tol', 1e-10);

  c = c(:);
  b = b(:);
  [m, n] = size(A);
  if numel(c) ~= n || numel(b) ~= m
    error('big_m_simplex:dim', 'Kich thuoc c, A, b khong tuong thich.');
  end

  % Dua ve b >= 0
  neg = b < 0;
  A(neg, :) = -A(neg, :);
  b(neg)    = -b(neg);

  % ---------- Them bien gia cho cac dong chua co cot don vi ----------
  B = zeros(1, m);
  for i = 1:m
    e = zeros(m, 1);
    e(i) = 1;
    j = find(all(abs(A - e) <= tol, 1), 1);
    if ~isempty(j), B(i) = j; end
  end
  rowsArt = find(B == 0);
  k = numel(rowsArt);
  art = n + (1:k);
  Aart = zeros(m, k);
  Aart(sub2ind([m, k], rowsArt, 1:k)) = 1;
  B(rowsArt) = art;
  AM = [A, Aart];
  N  = n + k;
  cM = [zeros(n, 1); ones(k, 1)];           % he so cua M trong ham muc tieu
  c0 = [c; zeros(k, 1)];                    % phan khong chua M

  % ---------- Truong hop M la mot so cu the: dung simplex.m cua cau 1 ----------
  if ~isempty(Mval)
    if exist('simplex', 'file') ~= 2
      error('big_m_simplex:path', ...
        'Can ham simplex.m cua cau 1: addpath(''../q01-simplex'').');
    end
    [xx, ~, exitflag, o] = simplex(c0 + Mval * cM, AM, b, B, opts);
    [x, fval, exitflag, out] = finish(xx, exitflag, n, art, c, tol, display);
    out.iterations = o.iterations;
    out.basis      = o.basis;
    out.tableau    = o.tableau;
    out.xHist      = o.xHist(:, 1:n);
    return;
  end

  % ---------- M la so lon tuy y: uoc luong Delta = alpha*M + beta ----------
  T = [AM, b];                              % B = cac cot don vi, nen B^-1 = I
  T = AM(:, B) \ T;
  exitflag = 0;
  iter     = 0;
  history  = {mat2str(sort(B))};
  xHist    = basic_solution(T, B, N)';
  xx       = [];

  while true
    alpha = cM(B)' * T(:, 1:N) - cM';       % he so cua M trong Delta_j
    beta  = c0(B)' * T(:, 1:N) - c0';       % phan khong chua M
    fM    = cM(B)' * T(:, end);
    f0    = c0(B)' * T(:, end);
    alpha(abs(alpha) <= tol) = 0;
    beta(abs(beta) <= tol)   = 0;

    if display
      print_tableau(T, B, cM, c0, alpha, beta, fM, f0, iter, n);
    end

    % Delta_j > 0 theo thu tu tu dien: alpha_j > 0, hoac alpha_j = 0 va beta_j > 0
    cand = find(alpha > 0 | (alpha == 0 & beta > 0));
    if isempty(cand)
      exitflag = 1;
      if display, fprintf('=> Moi Delta_j <= 0: bai toan M dat toi uu.\n\n'); end
      break;
    end
    if iter >= maxIter
      if display, fprintf('=> Vuot qua so buoc lap toi da (%d).\n\n', maxIter); end
      break;
    end

    % Chon bien vao co so
    if strcmpi(rule, 'bland')
      s = cand(1);
    else
      s = cand(1);                          % Delta lon nhat theo thu tu tu dien
      for j = cand(2:end)
        if alpha(j) > alpha(s) || (alpha(j) == alpha(s) && beta(j) > beta(s))
          s = j;
        end
      end
    end

    col = T(:, s);
    if all(col <= tol)
      exitflag = -3;
      if display
        fprintf('=> Delta_%d > 0 nhung moi z_i%d <= 0: bai toan M khong bi chan duoi.\n\n', s, s);
      end
      break;
    end

    % Chon bien ra theo ty so nho nhat
    ratio = inf(m, 1);
    pos = col > tol;
    ratio(pos) = T(pos, end) ./ col(pos);
    theta = min(ratio);
    ties = find(abs(ratio - theta) <= tol * max(1, abs(theta)));
    if strcmpi(rule, 'bland')
      [~, kk] = min(B(ties));
      r = ties(kk);
    else
      r = ties(1);
    end

    if display
      fprintf('   Bien vao: x%d (Delta = %s), bien ra: x%d (theta = %g), phan tu truc: %g\n\n', ...
        s, m_str(alpha(s), beta(s)), B(r), theta, T(r, s));
    end

    T(r, :) = T(r, :) / T(r, s);
    for i = [1:r-1, r+1:m]
      T(i, :) = T(i, :) - T(i, s) * T(r, :);
    end
    B(r) = s;
    iter = iter + 1;
    xHist(end+1, :) = basic_solution(T, B, N)';

    key = mat2str(sort(B));
    if any(strcmp(history, key))
      exitflag = -4;
      if display
        fprintf('=> Co so %s lap lai: xoay vong. Hay dung opts.rule = ''bland''.\n\n', key);
      end
      break;
    end
    history{end+1} = key;
  end

  xx = basic_solution(T, B, N);
  [x, fval, exitflag, out] = finish(xx, exitflag, n, art, c, tol, display);
  out.iterations = iter;
  out.basis      = B;
  out.alpha      = alpha;
  out.beta       = beta;
  out.tableau    = T;
  out.xHist      = xHist(:, 1:n);
end

% =====================================================================
function [x, fval, exitflag, out] = finish(xx, exitflag, n, art, c, tol, display)
% Ket luan cho bai toan goc tu phuong an cua bai toan M.
  xx(abs(xx) <= tol) = 0;                   % bo sai so lam tron
  x  = xx(1:n);
  xa = xx(art);
  fval = c' * x;
  if exitflag == 1 && any(xa > tol)
    exitflag = -2;
    fval = NaN;
    if display
      fprintf('=> Bai toan M toi uu nhung con bien gia duong: bai toan goc khong co phuong an chap nhan duoc.\n\n');
    end
  elseif exitflag == -3
    fval = -Inf;
  end
  out = struct('artificial', art, 'xa', xa);
end

function x = basic_solution(T, B, N)
  x = zeros(N, 1);
  x(B) = T(:, end);
end

function s = m_str(a, b)
% Bieu dien a*M + b dang chuoi, vi du "7M + 4", "-M", "0.2".
  if a == 0
    s = sprintf('%.4g', b);
    return;
  end
  if a == 1
    s = 'M';
  elseif a == -1
    s = '-M';
  else
    s = sprintf('%.4gM', a);
  end
  if b > 0
    s = sprintf('%s+%.4g', s, b);
  elseif b < 0
    s = sprintf('%s%.4g', s, b);
  end
end

function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end

function print_tableau(T, B, cM, c0, alpha, beta, fM, f0, iter, n)
  [m, n1] = size(T);
  N = n1 - 1;
  w = 10;
  line = repmat('-', 1, 29 + w * N);
  name = @(j) sprintf('x%d', j);
  fprintf('--- Bang don hinh (bai toan M), buoc lap %d ---\n', iter);
  fprintf('%8s %8s %10s |', 'c_B', 'Co so', 'Phuong an');
  for j = 1:N
    if j > n, fprintf('%*s', w, [name(j) '*']); else, fprintf('%*s', w, name(j)); end
  end
  fprintf('\n%8s %8s %10s |', '', '', '');
  for j = 1:N, fprintf('%*s', w, m_str(cM(j), c0(j))); end
  fprintf('\n%s\n', line);
  for i = 1:m
    fprintf('%8s %8s %10.4g |', m_str(cM(B(i)), c0(B(i))), name(B(i)), T(i, end));
    fprintf('%*.4g', [w * ones(1, N); T(i, 1:N)]);
    fprintf('\n');
  end
  fprintf('%s\n', line);
  fprintf('%8s %8s %10.4g |', '', 'Delta(M)', fM);
  fprintf('%*.4g', [w * ones(1, N); alpha]);
  fprintf('\n%8s %8s %10.4g |', '', 'Delta', f0);
  fprintf('%*.4g', [w * ones(1, N); beta]);
  fprintf('\n   (bien co dau * la bien gia; dong Delta(M) la he so cua M)\n\n');
end
