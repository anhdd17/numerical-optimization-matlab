function x = gs_min(f, x0, tol)
% GS_MIN  Gauss-Seidel (toa do tuan hoan) cho ham bat ky: cuc tieu lan luot
%   theo tung bien bang fminsearch mot bien. Dung khi ||x_k - x_{k-1}|| <= tol.
  x = x0(:);  n = numel(x);
  op = optimset('TolX', tol * 1e-3, 'TolFun', 1e-14, 'MaxFunEvals', 2000);
  fprintf('%4s  %-36s %14s %12s\n', 'k', 'x_k', 'f(x_k)', '||dx||');
  for k = 1:1000
    xo = x;
    for i = 1:n
      x(i) = fminsearch(@(t) f([x(1:i-1); t; x(i+1:end)]), x(i), op);
    end
    fprintf('%4d  %-36s %14.8g %12.3e\n', k, mat2str(x', 6), f(x), norm(x - xo));
    if norm(x - xo) <= tol, break; end
  end
end
