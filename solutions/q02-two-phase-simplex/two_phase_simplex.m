function [x, fval, exitflag, out] = two_phase_simplex(c, A, b, opts)
% TWO_PHASE_SIMPLEX  Giai bai toan quy hoach tuyen tinh dang chinh tac bang
%                    thuat toan don hinh hai pha.
%
%   min  f(x) = c'*x
%   s.t. A*x = b,  x >= 0
%
%   [x, fval, exitflag, out] = two_phase_simplex(c, A, b)
%   [x, fval, exitflag, out] = two_phase_simplex(c, A, b, opts)
%
%   Khong can co so xuat phat: pha 1 tu tim mot phuong an cuc bien ban dau
%   (hoac ket luan bai toan vo nghiem), pha 2 giai bai toan goc tu do.
%   Dung ham simplex.m cua cau 1 (solutions/q01-simplex) cho ca hai pha.
%
% Dau vao:
%   c, A, b - du lieu bai toan (A: m x n, khong can rank(A) = m)
%   opts    - (tuy chon) struct, chuyen tiep cho simplex.m:
%             .rule ('dantzig' | 'bland'), .display, .maxIter, .tol
%
% Dau ra:
%   x        - phuong an toi uu (n x 1)
%   fval     - gia tri toi uu c'*x
%   exitflag -  1 : tim duoc phuong an toi uu
%               0 : vuot qua so buoc lap toi da
%              -2 : bai toan khong co phuong an chap nhan duoc
%              -3 : ham muc tieu khong bi chan duoi
%              -4 : phat hien xoay vong
%   out      - struct:
%              .artificial   : chi so cac bien gia da them (cot n+1, n+2, ...)
%              .phase1       : out cua simplex o pha 1 (them .w = gia tri toi uu)
%              .phase2       : out cua simplex o pha 2 (rong neu dung o pha 1)
%              .removedRows  : cac dong (rang buoc) thua bi loai sau pha 1
%              .iterations   : tong so buoc lap hai pha
%              .xHist        : day phuong an (chi phan x goc) qua ca hai pha

  if nargin < 4, opts = struct(); end
  if exist('simplex', 'file') ~= 2
    error('two_phase_simplex:path', ...
      'Can ham simplex.m cua cau 1: addpath(''../q01-simplex'').');
  end
  tol     = get_opt(opts, 'tol', 1e-10);
  display = get_opt(opts, 'display', true);

  c = c(:);
  b = b(:);
  [m, n] = size(A);
  if numel(c) ~= n || numel(b) ~= m
    error('two_phase_simplex:dim', 'Kich thuoc c, A, b khong tuong thich.');
  end

  % Dua ve b >= 0
  neg = b < 0;
  A(neg, :) = -A(neg, :);
  b(neg)    = -b(neg);

  % ---------- Pha 1: them bien gia cho cac dong chua co cot don vi ----------
  B = zeros(1, m);
  for i = 1:m
    e = zeros(m, 1);
    e(i) = 1;
    j = find(all(abs(A - e) <= tol, 1), 1);
    if ~isempty(j), B(i) = j; end
  end
  rowsArt = reshape(find(B == 0), 1, []);    % 1 x k ca khi k = 0 (find tren 1 x 1 cho 0 x 0)
  k = numel(rowsArt);
  art = n + (1:k);
  Aart = zeros(m, k);
  Aart(sub2ind([m, k], rowsArt, 1:k)) = 1;
  B(rowsArt) = art;

  A1 = [A, Aart];
  c1 = [zeros(n, 1); ones(k, 1)];            % w = tong cac bien gia

  if display
    fprintf('########## PHA 1: min w = %s ##########\n', art_sum_str(art));
  end
  [~, w, flag1, out1] = simplex(c1, A1, b, B, opts);
  out1.w = w;

  out = struct('artificial', art, 'phase1', out1, 'phase2', [], ...
               'removedRows', [], 'iterations', out1.iterations, ...
               'xHist', out1.xHist(:, 1:n));

  if flag1 ~= 1                               % pha 1 luon bi chan duoi boi 0
    x = zeros(n, 1); fval = NaN; exitflag = flag1;
    return;
  end
  if w > tol
    x = zeros(n, 1); fval = NaN; exitflag = -2;
    if display
      fprintf('=> w* = %g > 0: bai toan goc khong co phuong an chap nhan duoc.\n\n', w);
    end
    return;
  end

  % ---------- Dua cac bien gia con lai (bang 0) ra khoi co so ----------
  T = out1.tableau;                           % T = B^-1 * [A1, b]
  B = out1.basis;
  keep = true(1, m);
  for r = find(B > n)
    j = find(abs(T(r, 1:n)) > tol & ~ismember(1:n, B), 1);
    if isempty(j)
      keep(r) = false;                        % dong r toan 0: rang buoc thua
      if display
        fprintf('   Bien gia x%d o lai co so, dong cua no bang 0 tren cac cot goc: loai rang buoc thua.\n', B(r));
      end
    else
      if display
        fprintf('   Xoay de dua bien gia x%d ra, x%d vao co so (theta = 0).\n', B(r), j);
      end
      T(r, :) = T(r, :) / T(r, j);
      for i = [1:r-1, r+1:m]
        T(i, :) = T(i, :) - T(i, j) * T(r, :);
      end
      B(r) = j;
    end
  end
  out.removedRows = find(~keep);

  % ---------- Pha 2: bai toan goc, xuat phat tu co so tim duoc ----------
  % He T(keep, 1:n) x = T(keep, end) tuong duong A x = b va co san cot don vi o B.
  A2 = T(keep, 1:n);
  b2 = T(keep, end);
  B2 = B(keep);
  if display
    fprintf('\n########## PHA 2: bai toan goc, co so xuat phat %s ##########\n', ...
      mat2str(B2));
  end
  [x, fval, exitflag, out2] = simplex(c, A2, b2, B2, opts);

  out.phase2     = out2;
  out.iterations = out.iterations + out2.iterations;
  out.xHist      = [out.xHist; out2.xHist];
end

% =====================================================================
function s = art_sum_str(art)
  s = strjoin(arrayfun(@(j) sprintf('x%d', j), art, 'UniformOutput', false), ' + ');
end

function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
