function x = heavy_ball(g, x0, alpha, beta, tol)
% HEAVY_BALL  Gradient quan tinh: x_{k+1} = x_k - alpha*grad f(x_k) + beta*(x_k - x_{k-1}).
%   Ham bac hai voi tri rieng mu..L: alpha = 4/(sqrt(L)+sqrt(mu))^2,
%   beta = ((sqrt(L)-sqrt(mu))/(sqrt(L)+sqrt(mu)))^2. Dung khi ||grad f|| <= tol.
  x = x0(:);  xp = x;
  fprintf('%5s  %-36s %12s\n', 'k', 'x_k', '||grad||');
  for k = 0:100000
    d = g(x);
    if k <= 5 || norm(d) <= tol, fprintf('%5d  %-36s %12.3e\n', k, mat2str(x', 6), norm(d)); end
    if norm(d) <= tol, break; end
    xn = x - alpha * d + beta * (x - xp);
    xp = x;  x = xn;
  end
end
