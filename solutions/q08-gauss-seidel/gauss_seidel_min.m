function [x, fval, out] = gauss_seidel_min(f, x0, opts)
% GAUSS_SEIDEL_MIN  Tim cuc tieu cua ham n bien f bang giai thuat Gauss-Seidel
%                   (phuong phap toa do tuan hoan - cyclic coordinate descent).
%
%   Moi vong lap k: lan luot voi i = 1, ..., n, cuc tieu f theo bien x_i,
%   cac bien con lai giu nguyen (gia tri moi nhat):
%       t* = argmin_t f(x + t e_i),   x <- x + t* e_i.
%   Bai toan mot chieu duoc giai bang: tim khoang chua cuc tieu
%   (bracket_minimum.m) + lat cat vang (golden_section.m cua cau 7).
%
%   [x, fval, out] = gauss_seidel_min(f, x0)
%   [x, fval, out] = gauss_seidel_min(f, x0, opts)
%
% Dau vao:
%   f    - ham muc tieu (handle), nhan vector cot n x 1
%   x0   - diem xuat phat (n x 1)
%   opts - (tuy chon) struct:
%          .tol     : dung khi ||x_k - x_{k-1}|| <= tol (mac dinh 1e-6)
%          .ftol    : hoac khi f(x_{k-1}) - f(x_k) <= ftol*(1 + |f|) (mac dinh 1e-12)
%          .maxIter : so vong lap toi da (mac dinh 500)
%          .lsTol   : do chinh xac tim kiem mot chieu (mac dinh 1e-8)
%          .h       : buoc ban dau khi tim khoang (mac dinh 1)
%          .coordMin: (tuy chon) handle @(x, i) tra ve gia tri toi uu cua x_i
%                     khi biet cong thuc (vi du ham bac hai); khi do khong
%                     dung tim kiem mot chieu
%          .display : in tung vong lap (mac dinh false)
%
% Dau ra:
%   x, fval - diem cuc tieu tim duoc va gia tri f
%   out     - struct: .iterations (so vong), .nfev (so lan tinh f),
%             .X (n x (K+1): diem sau moi vong), .path (n x (nK+1): diem sau
%             moi buoc toa do, de ve duong zic-zac), .F (f sau moi vong),
%             .exitflag (1: dat tol/ftol, 0: het maxIter)

  if nargin < 3, opts = struct(); end
  tol      = get_opt(opts, 'tol', 1e-6);
  ftol     = get_opt(opts, 'ftol', 1e-12);
  maxIter  = get_opt(opts, 'maxIter', 500);
  lsTol    = get_opt(opts, 'lsTol', 1e-8);
  h        = get_opt(opts, 'h', 1);
  coordMin = get_opt(opts, 'coordMin', []);
  display  = get_opt(opts, 'display', false);
  if isempty(coordMin) && exist('golden_section', 'file') ~= 2
    error('gauss_seidel_min:path', ...
      'Can golden_section.m cua cau 7: addpath(''../q07-golden-section-search'').');
  end

  x = x0(:);
  n = numel(x);
  fx = f(x);
  out = struct('iterations', 0, 'nfev', 1, 'X', x, 'path', x, 'F', fx, 'exitflag', 0);
  if display
    fprintf('%5s %-40s %14s %12s\n', 'k', 'x_k', 'f(x_k)', '||dx||');
    fprintf('%5d %-40s %14.8g\n', 0, mat2str(x', 6), fx);
  end

  for k = 1:maxIter
    xOld = x;  fOld = fx;
    for i = 1:n
      if ~isempty(coordMin)
        x(i) = coordMin(x, i);
      else
        ei = zeros(n, 1);  ei(i) = 1;
        phi = @(t) f(x + t * ei);
        [lo, hi, nb] = bracket_minimum(phi, h);
        [t, ~, o] = golden_section(phi, lo, hi, struct('tol', lsTol));
        out.nfev = out.nfev + nb + o.nfev;
        if phi(t) < fx                       % chi nhan buoc neu f giam
          x(i) = x(i) + t;
        end
        out.nfev = out.nfev + 1;
      end
      fx = f(x);
      out.nfev = out.nfev + 1;
      out.path(:, end+1) = x;
    end
    out.X(:, end+1) = x;
    out.F(end+1) = fx;
    out.iterations = k;
    dx = norm(x - xOld);
    if display
      fprintf('%5d %-40s %14.8g %12.3e\n', k, mat2str(x', 6), fx, dx);
    end
    if dx <= tol || (fOld - fx) <= ftol * (1 + abs(fx))
      out.exitflag = 1;
      break;
    end
  end
  fval = fx;
end

% =====================================================================
function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
