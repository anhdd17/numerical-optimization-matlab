function x = grad_desc(f, g, x0, tol, H)
% GRAD_DESC  Gradient: x_{k+1} = x_k - alpha_k * grad f(x_k), dung khi ||grad f|| <= tol.
%   Co H (ham bac hai 1/2 x'Hx + q'x): buoc chinh xac alpha = g'g / g'Hg.
%   Khong co H: tim kiem chinh xac bang fminbnd tren [0, amax] (amax nhan doi
%   den khi f tang).
  x = x0(:);
  fprintf('%5s  %-34s %14s %12s %10s\n', 'k', 'x_k', 'f(x_k)', '||grad||', 'alpha');
  for k = 0:100000
    d = g(x);
    if k <= 5 || norm(d) <= tol
      if k == 0, fprintf('%5d  %-34s %14.8g %12.3e\n', k, mat2str(x', 6), f(x), norm(d));
      else, fprintf('%5d  %-34s %14.8g %12.3e %10.4g\n', k, mat2str(x', 6), f(x), norm(d), a); end
    end
    if norm(d) <= tol, break; end
    if nargin >= 5
      a = (d' * d) / (d' * H * d);
    else
      phi = @(s) f(x - s * d);  amax = 1e-3;
      while phi(2 * amax) < phi(amax), amax = 2 * amax; end
      a = fminbnd(phi, 0, 2 * amax, optimset('TolX', 1e-12));
    end
    x = x - a * d;
  end
end
