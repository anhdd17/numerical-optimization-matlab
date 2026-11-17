function [x, fval, out] = newton_method(f, grad, hess, x0, opts)
% NEWTON_METHOD  Tim cuc tieu cua ham n bien f (kha vi hai lan) bang thuat toan Newton:
%
%       H_k d_k = -grad f(x_k),   x_{k+1} = x_k + alpha_k d_k,   H_k = hess f(x_k).
%
%   [x, fval, out] = newton_method(f, grad, hess, x0)
%   [x, fval, out] = newton_method(f, grad, hess, x0, opts)
%
% Dau vao:
%   f, grad, hess - ham, gradient (n x 1) va ma tran Hessian (n x n) (handle)
%   x0            - diem xuat phat (n x 1)
%   opts          - (tuy chon) struct:
%       .variant : 'pure'     - Newton thuan tuy, alpha_k = 1 (mac dinh)
%                  'damped'   - Newton giam buoc: alpha_k chon theo Armijo
%                               (thu 1, 1/2, 1/4, ...)
%                  'modified' - Newton cai bien: neu H_k khong xac dinh duong,
%                               thay bang H_k + tau*I (tang tau den khi phan
%                               tich Cholesky thanh cong), cong voi Armijo
%       .tol     : dung khi ||grad f(x_k)|| <= tol (mac dinh 1e-10)
%       .maxIter : so buoc lap toi da (mac dinh 100)
%       .c, .rho : tham so Armijo (mac dinh 1e-4, 0.5)
%       .display : in tung buoc (mac dinh false)
%
% Dau ra:
%   x, fval - diem tim duoc va gia tri f
%   out     - struct: .iterations, .X, .F, .G (||grad f||), .alphas, .taus,
%             .exitflag (1: dat tol; 0: het maxIter; -1: phan ky / H suy bien;
%             -2: 'damped' gap H khong xac dinh duong, d = -H\g khong la huong giam)

  if nargin < 5, opts = struct(); end
  variant = get_opt(opts, 'variant', 'pure');
  tol     = get_opt(opts, 'tol', 1e-10);
  maxIter = get_opt(opts, 'maxIter', 100);
  c       = get_opt(opts, 'c', 1e-4);
  rho     = get_opt(opts, 'rho', 0.5);
  display = get_opt(opts, 'display', false);

  x = x0(:);
  n = numel(x);
  fx = f(x);  g = grad(x);
  out = struct('iterations', 0, 'X', x, 'F', fx, 'G', norm(g), 'alphas', [], ...
               'taus', [], 'exitflag', 0);
  if display
    fprintf('%4s %-30s %16s %12s %8s %8s\n', 'k', 'x_k', 'f(x_k)', '||grad f||', 'alpha', 'tau');
    fprintf('%4d %-30s %16.10g %12.4e\n', 0, mat2str(x', 8), fx, norm(g));
  end

  for k = 1:maxIter
    if norm(g) <= tol
      out.exitflag = 1;
      break;
    end
    H = hess(x);
    tau = 0;
    if strcmpi(variant, 'modified')
      % Tang tau den khi H + tau*I xac dinh duong (Cholesky thanh cong)
      beta = 1e-3;
      if min(diag(H)) > 0, tau = 0; else, tau = -min(diag(H)) + beta; end
      while true
        [~, p] = chol(H + tau * eye(n));
        if p == 0, break; end
        tau = max(2 * tau, beta);
      end
      H = H + tau * eye(n);
    end
    if rcond(H) < 1e-14
      out.exitflag = -1;                     % Hessian suy bien: khong giai duoc
      break;
    end
    d = -H \ g;
    if ~strcmpi(variant, 'pure') && g' * d >= 0
      out.exitflag = -2;                     % d khong phai huong giam: Armijo bat luc
      break;
    end

    alpha = 1;
    if ~strcmpi(variant, 'pure')
      fa = f(x + d);
      while fa > fx + c * alpha * (g' * d) && alpha > 1e-12
        alpha = rho * alpha;
        fa = f(x + alpha * d);
      end
    end

    x = x + alpha * d;
    fx = f(x);  g = grad(x);
    out.X(:, end+1) = x;  out.F(end+1) = fx;  out.G(end+1) = norm(g);
    out.alphas(end+1) = alpha;  out.taus(end+1) = tau;
    out.iterations = k;
    if display
      fprintf('%4d %-30s %16.10g %12.4e %8.4g %8.3g\n', k, mat2str(x', 8), fx, norm(g), alpha, tau);
    end
    if ~all(isfinite(x)) || ~isfinite(fx)
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
