function [x, nu, lambda, level] = water_filling(alpha, P)
% WATER_FILLING  Giai bai toan phan bo cong suat (Boyd & Vandenberghe, muc 5.5.3)
%
%   min  -sum_i log(alpha_i + x_i)   s.t.  x >= 0,  sum_i x_i = P
%
%   bang cach giai TRUC TIEP he dieu kien KKT. Bai toan loi, rang buoc affine,
%   nen KKT la dieu kien can va du. Tu KKT:
%       x_i = max(0, L - alpha_i),  L = 1/nu  ("muc nuoc"),  sum_i x_i = P.
%   Ham L -> sum max(0, L - alpha_i) tang chat tren [min alpha, +inf), nen L duy
%   nhat; tim L chinh xac bang cach sap xep alpha (khong lap).
%
%   [x, nu, lambda, level] = water_filling(alpha, P)
%
% Dau vao:  alpha - vector n x 1, alpha_i > 0;  P > 0 - tong cong suat
% Dau ra:   x      - nghiem toi uu
%           nu     - nhan tu cua rang buoc dang thuc (nu = 1/L)
%           lambda - nhan tu cua x_i >= 0: lambda_i = nu - 1/(alpha_i + x_i)
%           level  - muc nuoc L

  alpha = alpha(:);
  n = numel(alpha);
  if any(alpha <= 0) || P <= 0
    error('water_filling:input', 'Can alpha_i > 0 va P > 0.');
  end
  a = sort(alpha);
  % Thu k = n, n-1, ..., 1 kenh duoc "do nuoc": L = (P + sum a(1:k))/k,
  % hop le khi L > a(k) (va L <= a(k+1) neu k < n).
  for k = n:-1:1
    L = (P + sum(a(1:k))) / k;
    if L > a(k)
      break;
    end
  end
  level = L;
  x = max(0, L - alpha);
  nu = 1 / L;
  lambda = nu - 1 ./ (alpha + x);
  lambda(abs(lambda) < 1e-14) = 0;
end
