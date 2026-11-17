function [x, fval, exitflag, out] = beale_qp(c, D, A, b, opts)
% BEALE_QP  Giai bai toan quy hoach toan phuong dang chinh tac
%
%       min f(x) = c'x + 1/2 x'Dx   s.t.  A x = b,  x >= 0      (D doi xung, nua xac dinh duong)
%
%   bang giai thuat don hinh Beale (1955, 1959). Moi buoc, cac bien co so duoc
%   bieu dien qua cac bien phi co so; f tro thanh ham bac hai cua cac bien phi
%   co so (bang 0 tai diem hien tai). Chon bien vao x_k co dao ham rieng
%   df/dx_k < 0 (hoac bien tu do u co df/du ~= 0) va tang no cho den khi:
%     (a) mot bien co so giam ve 0  (theta1, nhu don hinh)  => xoay co so;
%     (b) df/dx_k = 0               (theta2 = cuc tieu cua f theo huong)
%         => them BIEN TU DO u = 1/2 df/dx_k (khong rang buoc dau) va mot rang
%            buoc dinh nghia u; x_k vao co so tren dong moi, u phi co so (= 0).
%   Bien tu do nao vao co so thi bi loai cung dong dinh nghia cua no.
%   Dung khi df/dx_j >= 0 voi moi bien phi co so x_j va df/du = 0 voi moi bien tu do.
%
%   [x, fval, exitflag, out] = beale_qp(c, D, A, b)
%   [x, fval, exitflag, out] = beale_qp(c, D, A, b, opts)
%
% Dau vao:
%   c, D, A, b - du lieu bai toan (A: p x n)
%   opts       - (tuy chon) struct:
%       .B0      : co so xuat phat (chi so cot, A(:, B0) kha nghich, A(:,B0)\b >= 0).
%                  Bo trong: tim bang pha 1 cua two_phase_simplex.m (cau 2).
%       .tol     : sai so (mac dinh 1e-10)
%       .maxIter : so buoc toi da (mac dinh 500)
%       .display : in tung buoc (mac dinh false)
%
% Dau ra:
%   x, fval  - nghiem toi uu va gia tri f
%   exitflag - 1: toi uu; 0: het maxIter; -2: khong co phuong an; -3: f khong bi chan duoi
%   out      - struct: .iterations, .X (day diem), .F, .events (cell: 'xoay', 'tu do',
%              'loai u'), .nFree (so bien tu do da tao), .basis (co so cuoi, chi so
%              trong 1..n), .lambda (dao ham rieng df/dx_j tai nghiem = nhan tu cua x_j >= 0)

  if nargin < 5, opts = struct(); end
  tol     = get_opt(opts, 'tol', 1e-10);
  maxIter = get_opt(opts, 'maxIter', 500);
  display = get_opt(opts, 'display', false);
  c = c(:);  b = b(:);
  [p, n] = size(A);
  D = (D + D') / 2;
  out = struct('iterations', 0, 'X', zeros(n, 0), 'F', [], 'events', {{}}, 'nFree', 0, ...
               'basis', [], 'lambda', zeros(n, 1));

  % ----- co so xuat phat -----
  B = get_opt(opts, 'B0', []);
  if isempty(B)
    if exist('two_phase_simplex', 'file') ~= 2
      error('beale_qp:path', 'Can two_phase_simplex.m (cau 2): addpath(''../q01-simplex'', ''../q02-two-phase-simplex'').');
    end
    [~, ~, flag1, o1] = two_phase_simplex(zeros(n, 1), A, b, struct('display', false));
    if flag1 == -2
      x = zeros(n, 1);  fval = NaN;  exitflag = -2;
      return;
    end
    keep = setdiff(1:p, o1.removedRows);
    A = A(keep, :);  b = b(keep);  p = numel(keep);
    B = o1.phase2.basis;
  end
  B = B(:)';

  % He mo rong: M z = r, z = (x, u); cot n+t la bien tu do u_t (chi xuat hien o dong dinh nghia)
  M = A;  r = b;
  freeRow = zeros(1, 0);                     % freeRow(t) = dong dinh nghia cua bien tu do thu t
  freeId = zeros(1, 0);                      % freeId(t) = so thu tu (u1, u2, ...) khi tao
  exitflag = 0;
  if display
    fprintf('%4s %-28s %12s %-8s %-8s %10s %10s  %s\n', 'k', 'x_k', 'f(x_k)', 'vao', 'ra', 'theta1', 'theta2', 'buoc');
  end

  for k = 0:maxIter
    nz = size(M, 2);
    MB = M(:, B);
    z = zeros(nz, 1);  z(B) = MB \ r;
    z(B(B <= n & abs(z(B)) < tol)) = 0;
    x = z(1:n);
    fx = c' * x + 0.5 * x' * D * x;
    out.X(:, end+1) = x;  out.F(end+1) = fx;
    grad = c + D * x;

    % dao ham rieng cua f theo moi bien phi co so (cac bien phi co so khac giu = 0)
    N = setdiff(1:nz, B);
    G = zeros(1, numel(N));  Hd = G;  DX = zeros(n, numel(N));  DZB = zeros(numel(B), numel(N));
    for jj = 1:numel(N)
      j = N(jj);
      dzB = -(MB \ M(:, j));
      dz = zeros(nz, 1);  dz(j) = 1;  dz(B) = dzB;
      DX(:, jj) = dz(1:n);  DZB(:, jj) = dzB;
      G(jj) = dz(1:n)' * grad;
      Hd(jj) = dz(1:n)' * D * dz(1:n);
    end
    isFreeN = N > n;
    score = zeros(size(G));
    score(~isFreeN) = max(-G(~isFreeN), 0);  % bien x: can df/dx_j < 0
    score(isFreeN) = abs(G(isFreeN));        % bien tu do: can df/du ~= 0
    [best, jj] = max(score);
    if isempty(best) || best <= tol * max(1, norm(grad))
      exitflag = 1;
      if display, fprintf('%4d %-28s %12.8g   toi uu\n', k, mat2str(x' + 0, 6), fx); end
      break;
    end
    if k == maxIter, break; end

    j = N(jj);
    sigma = 1;
    if isFreeN(jj), sigma = -sign(G(jj)); end
    dzB = sigma * DZB(:, jj);
    % (a) ti so: bien co so x_i (khong am) giam ve 0
    theta1 = Inf;  ileave = 0;
    for ii = 1:numel(B)
      if B(ii) <= n && dzB(ii) < -tol
        th = z(B(ii)) / (-dzB(ii));
        if ileave == 0 || th < theta1 - tol || (abs(th - theta1) <= tol && B(ii) < B(ileave))
          theta1 = th;  ileave = ii;
        end
      end
    end
    % (b) cuc tieu cua f theo huong
    if Hd(jj) > tol, theta2 = abs(G(jj)) / Hd(jj); else, theta2 = Inf; end
    if isinf(theta1) && isinf(theta2)
      exitflag = -3;
      if display, fprintf('%4d %-28s %12.8g   f khong bi chan duoi\n', k, mat2str(x' + 0, 6), fx); end
      break;
    end

    vname = var_name(j, n, freeId);
    if theta1 <= theta2
      lname = var_name(B(ileave), n, freeId);
      B(ileave) = j;
      ev = 'xoay';
    else
      lname = '-';
      % bien tu do moi: u = 1/2 df/dz_j = 1/2 dx_j'(c + D x)  <=>  1/2 dx_j'D x - u = -1/2 dx_j'c
      dxj = DX(:, jj);
      M = [M, zeros(size(M, 1), 1); [0.5 * dxj' * D, zeros(1, nz - n), -1]];
      r = [r; -0.5 * dxj' * c];
      out.nFree = out.nFree + 1;
      freeRow(end+1) = size(M, 1);
      freeId(end+1) = out.nFree;
      B = [B, j];                                % z_j vao co so tren dong moi, u = 0 phi co so
      ev = sprintf('tao u%d', out.nFree);
    end
    % bien tu do nao nam trong co so thi loai cung dong dinh nghia cua no
    % (u khong rang buoc dau nen dong do khong han che gi cac bien con lai)
    jf = B(B > n);
    for jdel = sort(jf, 'descend')
      t = jdel - n;
      rowDel = freeRow(t);
      ev = sprintf('%s, loai u%d', ev, freeId(t));
      M(rowDel, :) = [];  r(rowDel) = [];
      M(:, jdel) = [];
      B(B == jdel) = [];
      B(B > jdel) = B(B > jdel) - 1;
      freeRow(t) = [];  freeId(t) = [];
      freeRow(freeRow > rowDel) = freeRow(freeRow > rowDel) - 1;
    end
    out.events{end+1} = ev;
    if display
      fprintf('%4d %-28s %12.8g %-8s %-8s %10.4g %10.4g  %s\n', k, mat2str(x' + 0, 6), fx, vname, lname, ...
        theta1, theta2, ev);
    end
    out.iterations = k + 1;
  end
  fval = fx;
  out.basis = B(B <= n);
  out.lambda(N(~isFreeN)) = G(~isFreeN);
end

% =====================================================================
function s = var_name(j, n, freeId)
  if j <= n, s = sprintf('x%d', j); else, s = sprintf('u%d', freeId(j - n)); end
end

function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
