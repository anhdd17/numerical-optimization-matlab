function x = newton(g, h, x0, tol)
% NEWTON  Newton thuan tuy: giai H(x_k) d = -grad f(x_k), x_{k+1} = x_k + d.
%   Dung khi ||grad f|| <= tol.
  x = x0(:);
  fprintf('%4s  %-40s %12s\n', 'k', 'x_k', '||grad||');
  for k = 0:100
    d = g(x);
    fprintf('%4d  %-40s %12.3e\n', k, mat2str(x', 6), norm(d));
    if norm(d) <= tol, break; end
    x = x - h(x) \ d;
  end
end
