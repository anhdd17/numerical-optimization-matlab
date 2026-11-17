function [x, fval, out] = momentum_gradient(f, grad, x0, opts)
% MOMENTUM_GRADIENT  Tim cuc tieu cua ham n bien kha vi f bang thuat toan gradient
%                    quan tinh (heavy-ball, Polyak 1964):
%
%       x_{k+1} = x_k - alpha * grad f(x_k) + beta * (x_k - x_{k-1}),   x_{-1} = x_0.
%
%   Tuong duong voi dang "van toc":  v_{k+1} = beta*v_k - alpha*grad f(x_k),
%   x_{k+1} = x_k + v_{k+1},  v_0 = 0.
%
%   [x, fval, out] = momentum_gradient(f, grad, x0, opts)
%
% Dau vao:
%   f, grad - ham muc tieu va gradient (handle)
%   x0      - diem xuat phat (n x 1)
%   opts    - struct:
%             .alpha   : buoc (bat buoc)
%             .beta    : he so quan tinh, 0 <= beta < 1 (mac dinh 0.9);
%                        beta = 0 la thuat toan gradient buoc co dinh
%             .tol     : dung khi ||grad f(x_k)|| <= tol (mac dinh 1e-6)
%             .maxIter : so buoc lap toi da (mac dinh 10000)
%
% Dau ra:
%   x, fval - diem tim duoc va gia tri f
%   out     - struct: .iterations, .ngev, .X, .F, .G, .exitflag
%             (1: dat tol; 0: het maxIter; -1: phan ky)

  if nargin < 4 || ~isfield(opts, 'alpha')
    error('momentum_gradient:alpha', 'Can opts.alpha.');
  end
  alpha   = opts.alpha;
  beta    = get_opt(opts, 'beta', 0.9);
  tol     = get_opt(opts, 'tol', 1e-6);
  maxIter = get_opt(opts, 'maxIter', 10000);

  x = x0(:);
  xPrev = x;
  g = grad(x);
  out = struct('iterations', 0, 'ngev', 1, 'X', x, 'F', f(x), 'G', norm(g), 'exitflag', 0);
  for k = 1:maxIter
    if norm(g) <= tol
      out.exitflag = 1;
      break;
    end
    xNew = x - alpha * g + beta * (x - xPrev);
    xPrev = x;
    x = xNew;
    g = grad(x);
    out.ngev = out.ngev + 1;
    out.X(:, end+1) = x;  out.F(end+1) = f(x);  out.G(end+1) = norm(g);
    out.iterations = k;
    if ~all(isfinite(x))
      out.exitflag = -1;
      break;
    end
  end
  if out.exitflag == 0 && norm(g) <= tol
    out.exitflag = 1;
  end
  fval = f(x);
end

% =====================================================================
function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
