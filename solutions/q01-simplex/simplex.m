function [x, fval, exitflag, out] = simplex(c, A, b, B, opts)
% SIMPLEX  Giai bai toan quy hoach tuyen tinh dang chinh tac bang
%          thuat toan don hinh (bang don hinh).
%
%   min  f(x) = c'*x
%   s.t. A*x = b,  x >= 0
%
%   [x, fval, exitflag, out] = simplex(c, A, b)
%   [x, fval, exitflag, out] = simplex(c, A, b, B)
%   [x, fval, exitflag, out] = simplex(c, A, b, B, opts)
%
% Dau vao:
%   c    - vector he so ham muc tieu (n x 1)
%   A    - ma tran rang buoc (m x n), rank(A) = m
%   b    - vector ve phai (m x 1)
%   B    - (tuy chon) chi so cac bien co so ban dau (1 x m). Neu bo trong,
%          chuong trinh tu tim m cot don vi trong A lam co so xuat phat.
%   opts - (tuy chon) struct gom cac truong:
%          .rule    : 'dantzig' (mac dinh) hoac 'bland' (chong xoay vong)
%          .display : true (mac dinh) in bang don hinh moi buoc
%          .maxIter : so buoc lap toi da (mac dinh 100)
%          .tol     : sai so so hoc (mac dinh 1e-10)
%
% Dau ra:
%   x        - phuong an toi uu (n x 1)
%   fval     - gia tri toi uu c'*x
%   exitflag -  1 : tim duoc phuong an toi uu
%               0 : vuot qua so buoc lap toi da
%              -3 : ham muc tieu khong bi chan duoi
%              -4 : phat hien xoay vong (lap lai co so)
%   out      - struct: .iterations, .basis, .delta, .tableau, .xHist,
%              .direction (huong giam vo han khi exitflag = -3)

  % ---------- Xu ly tham so ----------
  if nargin < 4, B = []; end
  if nargin < 5, opts = struct(); end
  rule    = get_opt(opts, 'rule', 'dantzig');
  display = get_opt(opts, 'display', true);
  maxIter = get_opt(opts, 'maxIter', 100);
  tol     = get_opt(opts, 'tol', 1e-10);

  c = c(:);
  b = b(:);
  [m, n] = size(A);
  if numel(c) ~= n || numel(b) ~= m
    error('simplex:dim', 'Kich thuoc c, A, b khong tuong thich.');
  end

  % Dua ve b >= 0 (nhan -1 vao cac rang buoc co b_i < 0)
  neg = b < 0;
  A(neg, :) = -A(neg, :);
  b(neg)    = -b(neg);

  % Tim co so xuat phat neu chua cho
  if isempty(B)
    B = find_identity_basis(A, tol);
    if isempty(B)
      error('simplex:nobasis', ...
        ['Khong tim thay co so don vi trong A. Hay cung cap B ' ...
         'hoac dung thuat toan don hinh hai pha (cau 2).']);
    end
  end
  B = B(:)';
  if numel(B) ~= m
    error('simplex:basis', 'Co so B phai co dung m = %d phan tu.', m);
  end

  % ---------- Lap bang don hinh ban dau: T = [B^-1*A, B^-1*b] ----------
  AB = A(:, B);
  if rank(AB) < m
    error('simplex:singular', 'Ma tran co so A(:,B) suy bien.');
  end
  T = AB \ [A, b];
  if any(T(:, end) < -tol)
    error('simplex:infeasible', ...
      'Co so B khong cho phuong an chap nhan duoc (B^-1*b co phan tu am).');
  end

  exitflag = 0;
  iter     = 0;
  history  = {mat2str(sort(B))};
  xHist    = basic_solution(T, B, n)';
  out      = struct('direction', []);

  % ---------- Vong lap don hinh ----------
  while true
    cB    = c(B);
    delta = (cB' * T(:, 1:n)) - c';   % uoc luong Delta_j = c_B'*B^-1*A_j - c_j
    fcur  = cB' * T(:, end);

    if display
      print_tableau(T, B, c, delta, fcur, iter);
    end

    % Dau hieu toi uu: Delta_j <= 0 voi moi j
    cand = find(delta > tol);
    if isempty(cand)
      exitflag = 1;
      if display, fprintf('=> Dat dieu kien toi uu: moi Delta_j <= 0.\n\n'); end
      break;
    end
    if iter >= maxIter
      exitflag = 0;
      if display, fprintf('=> Vuot qua so buoc lap toi da (%d).\n\n', maxIter); end
      break;
    end

    % Chon bien vao co so (cot xoay s)
    if strcmpi(rule, 'bland')
      s = cand(1);                    % chi so nho nhat co Delta_j > 0
    else
      [~, k] = max(delta(cand));      % Delta_j lon nhat (quy tac Dantzig)
      s = cand(k);
    end

    % Dau hieu khong bi chan: cot z_s <= 0
    col = T(:, s);
    if all(col <= tol)
      exitflag = -3;
      d = zeros(n, 1);
      d(s) = 1;
      d(B) = -col;
      out.direction = d;
      if display
        fprintf('=> Delta_%d > 0 nhung moi z_i%d <= 0: ham muc tieu khong bi chan duoi.\n\n', s, s);
      end
      break;
    end

    % Chon bien ra khoi co so (dong xoay r) theo ty so nho nhat
    ratio = inf(m, 1);
    pos = col > tol;
    ratio(pos) = T(pos, end) ./ col(pos);
    theta = min(ratio);
    ties = find(abs(ratio - theta) <= tol * max(1, abs(theta)));
    if strcmpi(rule, 'bland')
      [~, k] = min(B(ties));          % bien co chi so nho nhat roi co so
      r = ties(k);
    else
      r = ties(1);
    end

    if display
      fprintf('   Bien vao: x%d (Delta = %g), bien ra: x%d (theta = %g), phan tu truc: %g\n\n', ...
        s, delta(s), B(r), theta, T(r, s));
    end

    % Phep xoay quanh phan tu T(r,s)
    T(r, :) = T(r, :) / T(r, s);
    for i = [1:r-1, r+1:m]
      T(i, :) = T(i, :) - T(i, s) * T(r, :);
    end
    B(r) = s;
    iter = iter + 1;
    xHist(end+1, :) = basic_solution(T, B, n)';

    % Phat hien xoay vong (chi xay ra khi suy bien va dung quy tac Dantzig)
    key = mat2str(sort(B));
    if any(strcmp(history, key))
      exitflag = -4;
      if display
        fprintf('=> Co so %s lap lai: thuat toan bi xoay vong. Hay dung opts.rule = ''bland''.\n\n', key);
      end
      break;
    end
    history{end+1} = key;
  end

  x    = basic_solution(T, B, n);
  fval = c' * x;

  out.iterations = iter;
  out.basis      = B;
  out.delta      = delta;
  out.tableau    = T;
  out.xHist      = xHist;
end

% =====================================================================
function x = basic_solution(T, B, n)
  x = zeros(n, 1);
  x(B) = T(:, end);
end

function B = find_identity_basis(A, tol)
  [m, ~] = size(A);
  B = zeros(1, m);
  for i = 1:m
    e = zeros(m, 1);
    e(i) = 1;
    j = find(all(abs(A - e) <= tol, 1), 1);
    if isempty(j)
      B = [];
      return;
    end
    B(i) = j;
  end
end

function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end

function print_tableau(T, B, c, delta, f, iter)
  [m, n1] = size(T);
  n = n1 - 1;
  fprintf('--- Bang don hinh, buoc lap %d ---\n', iter);
  fprintf('%6s %8s %10s |', 'c_B', 'Co so', 'Phuong an');
  for j = 1:n, fprintf('%9s', sprintf('x%d', j)); end
  fprintf('\n');
  fprintf('%6s %8s %10s |', '', '', '');
  for j = 1:n, fprintf('%9.4g', c(j)); end
  fprintf('\n');
  fprintf('%s\n', repmat('-', 1, 27 + 9 * n));
  for i = 1:m
    fprintf('%6.4g %8s %10.4g |', c(B(i)), sprintf('x%d', B(i)), T(i, end));
    fprintf('%9.4g', T(i, 1:n));
    fprintf('\n');
  end
  fprintf('%s\n', repmat('-', 1, 27 + 9 * n));
  fprintf('%6s %8s %10.4g |', '', 'Delta', f);
  fprintf('%9.4g', delta);
  fprintf('\n\n');
end
