function [x, fval, out] = quasi_newton(f, grad, x0, opts)
% QUASI_NEWTON  Tim cuc tieu cua ham n bien kha vi f bang thuat toan tua Newton:
%
%       d_k = -H_k * grad f(x_k),   x_{k+1} = x_k + alpha_k d_k,
%
%   trong do H_k xap xi nghich dao cua Hessian, cap nhat hang hai tu
%   s_k = x_{k+1} - x_k, y_k = grad f(x_{k+1}) - grad f(x_k) sao cho
%   H_{k+1} y_k = s_k (phuong trinh cat tuyen).
%
%   [x, fval, out] = quasi_newton(f, grad, x0)
%   [x, fval, out] = quasi_newton(f, grad, x0, opts)
%
% Dau vao:
%   f, grad - ham muc tieu va gradient (handle), nhan vector cot n x 1
%   x0      - diem xuat phat (n x 1)
%   opts    - (tuy chon) struct:
%       .method     : 'bfgs' (mac dinh) hoac 'dfp'
%       .lineSearch : 'wolfe' - dieu kien Wolfe manh, wolfe_line_search.m (mac dinh)
%                     'exact' - cuc tieu phi(a) = f(x + a d): tim khoang bang
%                               bracket_minimum.m (cau 8), roi chia doi theo dau
%                               phi'(a) = grad f(x + a d)'d bang bisection_min.m
%                               (cau 6) den do chinh xac may
%       .H0         : ma tran ban dau (mac dinh I)
%       .scaleH0    : truoc lan cap nhat dau tien thay H0 boi (s'y / y'y) I
%                     (mac dinh true neu khong cho .H0)
%       .tol        : dung khi ||grad f(x_k)|| <= tol (mac dinh 1e-8)
%       .maxIter    : so buoc lap toi da (mac dinh 1000)
%       .c1, .c2    : tham so Wolfe (mac dinh 1e-4, 0.9)
%       .display    : in tung buoc (mac dinh false)
%
% Dau ra:
%   x, fval - diem tim duoc va gia tri f
%   out     - struct: .iterations, .nfev, .ngev, .X, .F, .G (||grad f||),
%             .alphas, .S, .Y (cac vector s_k, y_k), .H (ma tran H cuoi cung),
%             .skipped (so lan bo qua cap nhat vi s'y <= 0),
%             .exitflag (1: dat tol; 0: het maxIter; -1: f khong huu han)

  if nargin < 4, opts = struct(); end
  method  = lower(get_opt(opts, 'method', 'bfgs'));
  lsearch = lower(get_opt(opts, 'lineSearch', 'wolfe'));
  tol     = get_opt(opts, 'tol', 1e-8);
  maxIter = get_opt(opts, 'maxIter', 1000);
  c1      = get_opt(opts, 'c1', 1e-4);
  c2      = get_opt(opts, 'c2', 0.9);
  display = get_opt(opts, 'display', false);
  if ~any(strcmp(method, {'bfgs', 'dfp'}))
    error('quasi_newton:method', 'opts.method phai la ''bfgs'' hoac ''dfp''.');
  end
  if strcmp(lsearch, 'exact') && (exist('bisection_min', 'file') ~= 2 || exist('bracket_minimum', 'file') ~= 2)
    error('quasi_newton:path', ['Can bisection_min.m (cau 6) va bracket_minimum.m (cau 8): ' ...
      'addpath(''../q06-bisection-search'', ''../q08-gauss-seidel'').']);
  end

  x = x0(:);
  n = numel(x);
  if isfield(opts, 'H0') && ~isempty(opts.H0)
    H = opts.H0;
    scale = get_opt(opts, 'scaleH0', false);
  else
    H = eye(n);
    scale = get_opt(opts, 'scaleH0', true);
  end
  fx = f(x);  g = grad(x);
  out = struct('iterations', 0, 'nfev', 1, 'ngev', 1, 'X', x, 'F', fx, 'G', norm(g), ...
               'alphas', [], 'S', zeros(n, 0), 'Y', zeros(n, 0), 'H', H, 'skipped', 0, 'exitflag', 0);
  if display
    fprintf('%4s %-30s %16s %12s %10s\n', 'k', 'x_k', 'f(x_k)', '||grad f||', 'alpha_k');
    fprintf('%4d %-30s %16.10g %12.4e\n', 0, mat2str(x', 8), fx, norm(g));
  end

  for k = 1:maxIter
    if norm(g) <= tol
      out.exitflag = 1;
      break;
    end
    d = -H * g;
    if g' * d >= 0                             % H mat xac dinh duong (sai so): khoi dong lai
      H = eye(n);  d = -g;
    end

    if strcmp(lsearch, 'wolfe')
      [a, fn, gn, nf, ng] = wolfe_line_search(f, grad, x, d, fx, g, 1, c1, c2);
      out.nfev = out.nfev + nf;  out.ngev = out.ngev + ng;
    else
      phi = @(t) f(x + t * d);
      [lo, hi, nb] = bracket_minimum(phi, 1);
      lo = max(lo, 0);                         % chi xet a >= 0 (d la huong giam)
      dphi = @(t) grad(x + t * d)' * d;
      [a, fn, o] = bisection_min(phi, lo, hi, struct('method', 'derivative', 'df', dphi, ...
                                  'tol', 4 * eps * max(1, hi), 'maxIter', 200));
      gn = grad(x + a * d);
      out.nfev = out.nfev + nb + o.nfev;  out.ngev = out.ngev + o.ndfev + 1;
    end

    s = a * d;
    xn = x + s;
    y = gn - g;
    sy = s' * y;
    if sy > 1e-12 * norm(s) * norm(y)
      if scale && k == 1
        H = (sy / (y' * y)) * eye(n);          % uoc luong do cong ban dau
      end
      if strcmp(method, 'bfgs')
        rho = 1 / sy;
        V = eye(n) - rho * (s * y');
        H = V * H * V' + rho * (s * s');
      else
        Hy = H * y;
        H = H - (Hy * Hy') / (y' * Hy) + (s * s') / sy;
      end
    else
      out.skipped = out.skipped + 1;           % s'y <= 0: giu nguyen H de H xac dinh duong
    end

    x = xn;  fx = fn;  g = gn;
    out.X(:, end+1) = x;  out.F(end+1) = fx;  out.G(end+1) = norm(g);
    out.alphas(end+1) = a;  out.S(:, end+1) = s;  out.Y(:, end+1) = y;
    out.iterations = k;
    if display
      fprintf('%4d %-30s %16.10g %12.4e %10.4g\n', k, mat2str(x', 8), fx, norm(g), a);
    end
    if ~isfinite(fx)
      out.exitflag = -1;
      break;
    end
  end
  if out.exitflag == 0 && norm(g) <= tol
    out.exitflag = 1;
  end
  out.H = H;
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
