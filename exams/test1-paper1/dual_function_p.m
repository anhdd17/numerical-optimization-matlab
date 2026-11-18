function [g, x] = dual_function_p(nu, lambda, a)
% DUAL_FUNCTION_P  Ham doi ngau Lagrange cua bai toan (P) (phan B):
%
%   min sum_i x_i ln(a + x_i)   s.t.  sum_i x_i = 1,  x_i >= 0
%
%   L(x, nu, lambda) = sum_i [x_i ln(a + x_i) - (nu + lambda_i) x_i] + nu.
%   Voi mu_i = nu + lambda_i, inf theo x_i > -a dat tai x_i = a/w_i - a,
%   w_i = W(a e^(1 - mu_i)) (W: ham Lambert), gia tri -a (w_i - 1)^2 / w_i, nen
%
%       g(nu, lambda) = nu - a * sum_i (w_i - 1)^2 / w_i.
%
%   [g, x] = dual_function_p(nu, lambda, a)
%   nu     - nhan tu cua rang buoc dang thuc (vo huong)
%   lambda - nhan tu cua x_i >= 0 (n x 1, lambda >= 0)
%   g      - gia tri ham doi ngau; x - diem dat inf cua L (n x 1)
  lambda = lambda(:);
  w = lambert_w0(a * exp(1 - nu - lambda));
  g = nu - a * sum((w - 1).^2 ./ w);
  x = a ./ w - a;
end
