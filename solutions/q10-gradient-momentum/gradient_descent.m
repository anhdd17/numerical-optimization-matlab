function [x, fval, out] = gradient_descent(f, grad, x0, opts)
% GRADIENT_DESCENT  Tim cuc tieu cua ham n bien kha vi f bang thuat toan gradient:
%                   x_{k+1} = x_k - alpha_k * grad f(x_k).
%
%   [x, fval, out] = gradient_descent(f, grad, x0)
%   [x, fval, out] = gradient_descent(f, grad, x0, opts)
%
% Dau vao:
%   f, grad - ham muc tieu va gradient (handle), nhan vector cot n x 1
%   x0      - diem xuat phat (n x 1)
%   opts    - (tuy chon) struct:
%             .step    : cach chon do dai buoc
%                        'fixed'  - alpha_k = opts.alpha (mac dinh)
%                        'armijo' - quay lui: bat dau tu opts.alpha, nhan rho
%                                   cho den khi f(x - a g) <= f(x) - c a ||g||^2
%                        'exact'  - tim kiem chinh xac theo tia: cuc tieu
%                                   phi(a) = f(x - a g) bang bracket_minimum.m
%                                   (cau 8) + golden_section.m (cau 7)
%             .alpha   : buoc co dinh / buoc thu dau tien (mac dinh 1)
%             .rho, .c : tham so Armijo (mac dinh 0.5, 1e-4)
%             .tol     : dung khi ||grad f(x_k)|| <= tol (mac dinh 1e-6)
%             .maxIter : so buoc lap toi da (mac dinh 10000)
%             .display : in tung buoc (mac dinh false)
%
% Dau ra:
%   x, fval - diem tim duoc va gia tri f
%   out     - struct: .iterations, .nfev, .ngev (so lan tinh f, grad f),
%             .X (n x (K+1) day diem), .F, .G (f va ||grad f|| theo buoc),
%             .alphas (cac buoc da dung), .exitflag (1: dat tol; 0: het maxIter;
%             -1: phan ky - f khong con huu han)

  if nargin < 4, opts = struct(); end
  step    = get_opt(opts, 'step', 'fixed');
  alpha0  = get_opt(opts, 'alpha', 1);
  rho     = get_opt(opts, 'rho', 0.5);
  c       = get_opt(opts, 'c', 1e-4);
  tol     = get_opt(opts, 'tol', 1e-6);
  maxIter = get_opt(opts, 'maxIter', 10000);
  display = get_opt(opts, 'display', false);
  if strcmpi(step, 'exact') && (exist('golden_section', 'file') ~= 2 || exist('bracket_minimum', 'file') ~= 2)
    error('gradient_descent:path', ['Can golden_section.m (cau 7) va bracket_minimum.m (cau 8): ' ...
      'addpath(''../q07-golden-section-search'', ''../q08-gauss-seidel'').']);
  end

  x = x0(:);
  fx = f(x);  g = grad(x);
  out = struct('iterations', 0, 'nfev', 1, 'ngev', 1, 'X', x, 'F', fx, ...
               'G', norm(g), 'alphas', [], 'exitflag', 0);
  if display
    fprintf('%5s %-28s %14s %12s %10s\n', 'k', 'x_k', 'f(x_k)', '||grad f||', 'alpha_k');
    fprintf('%5d %-28s %14.8g %12.4e\n', 0, mat2str(x', 6), fx, norm(g));
  end

  for k = 1:maxIter
    if norm(g) <= tol
      out.exitflag = 1;
      break;
    end
    d = -g;
    switch lower(step)
      case 'fixed'
        a = alpha0;
      case 'armijo'
        a = alpha0;
        fa = f(x + a * d);  out.nfev = out.nfev + 1;
        while fa > fx - c * a * (g' * g)
          a = rho * a;
          fa = f(x + a * d);  out.nfev = out.nfev + 1;
        end
      case 'exact'
        phi = @(t) f(x + t * d);
        [lo, hi, nb] = bracket_minimum(phi, alpha0);
        lo = max(lo, 0);                       % chi xet a >= 0 (d la huong giam)
        [a, ~, o] = golden_section(phi, lo, hi, struct('tol', 1e-10 * max(1, hi)));
        out.nfev = out.nfev + nb + o.nfev;
      otherwise
        error('gradient_descent:step', 'opts.step phai la ''fixed'', ''armijo'' hoac ''exact''.');
    end
    x = x + a * d;
    fx = f(x);  g = grad(x);
    out.nfev = out.nfev + 1;  out.ngev = out.ngev + 1;
    out.X(:, end+1) = x;  out.F(end+1) = fx;  out.G(end+1) = norm(g);
    out.alphas(end+1) = a;
    out.iterations = k;
    if display
      fprintf('%5d %-28s %14.8g %12.4e %10.4g\n', k, mat2str(x', 6), fx, norm(g), a);
    end
    if ~isfinite(fx)
      out.exitflag = -1;
      break;
    end
  end
  if out.exitflag == 0 && norm(g) <= tol
    out.exitflag = 1;
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
