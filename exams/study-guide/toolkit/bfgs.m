function x = bfgs(f, g, x0, tol)
% BFGS  Tua Newton: d = -H*grad, buoc Armijo (thu 1, 1/2, ...), cap nhat
%   H+ = (I - r s y') H (I - r y s') + r s s',  r = 1/(y's)  (bo qua neu y's <= 0).
%   Dung khi ||grad f|| <= tol.
  x = x0(:);  n = numel(x);  H = eye(n);  gx = g(x);
  fprintf('%4s  %-40s %14s %12s\n', 'k', 'x_k', 'f(x_k)', '||grad||');
  for k = 0:1000
    fprintf('%4d  %-40s %14.8g %12.3e\n', k, mat2str(x', 6), f(x), norm(gx));
    if norm(gx) <= tol, break; end
    d = -H * gx;  a = 1;
    while f(x + a * d) > f(x) + 1e-4 * a * (gx' * d), a = a / 2; end
    s = a * d;  xn = x + s;  gn = g(xn);  y = gn - gx;
    if y' * s > 1e-12
      r = 1 / (y' * s);  I = eye(n);
      H = (I - r * (s * y')) * H * (I - r * (y * s')) + r * (s * s');
    end
    x = xn;  gx = gn;
  end
end
