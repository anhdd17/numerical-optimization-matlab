function x = gs_quad(H, q, x0, tol)
% GS_QUAD  Gauss-Seidel cho f = 1/2 x'Hx + q'x (H doi xung xac dinh duong).
%   Moi vong: x_i = -(q_i + sum_{j~=i} H_ij x_j) / H_ii, i = 1..n.
%   Dung khi ||x_k - x_{k-1}|| <= tol.
  x = x0(:);  n = numel(x);
  fprintf('%4s  %-36s %12s\n', 'k', 'x_k', '||dx||');
  fprintf('%4d  %-36s\n', 0, mat2str(x', 6));
  for k = 1:1000
    xo = x;
    for i = 1:n
      x(i) = -(q(i) + H(i, :) * x - H(i, i) * x(i)) / H(i, i);
    end
    fprintf('%4d  %-36s %12.3e\n', k, mat2str(x', 6), norm(x - xo));
    if norm(x - xo) <= tol, break; end
  end
end
