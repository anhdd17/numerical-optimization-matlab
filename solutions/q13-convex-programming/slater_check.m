function [ok, xhat, sstar] = slater_check(prob, x0, isAffine)
% SLATER_CHECK  Kiem tra dieu kien Slater cho bai toan quy hoach loi
%
%   min f(x)   s.t.  g_i(x) <= 0 (g_i loi),   h_j(x) = 0 (h_j affine)
%
%   bang cach giai bai toan phu (bien (x, s)):
%
%       min s   s.t.  g_i(x) <= s  (g_i khong affine),   g_i(x) <= 0  (g_i affine),
%                     h(x) = 0,   s >= -1.
%
%   Slater (dang tinh: rang buoc affine chi can thoa, khong can thoa chat) dung
%   khi va chi khi gia tri toi uu s* < 0; khi do xhat la mot diem Slater.
%   Chan duoi s >= -1 giu cho bai toan phu bi chan.
%
%   [ok, xhat, sstar] = slater_check(prob, x0)
%   [ok, xhat, sstar] = slater_check(prob, x0, isAffine)
%
% Dau vao:
%   prob     - struct: .g (x -> m x 1), .h (x -> p x 1) (handle, co the bo trong)
%   x0       - diem xuat phat cho solver (n x 1)
%   isAffine - (tuy chon) vector logic m x 1: rang buoc g_i nao la affine.
%              Mac dinh: tat ca la false (Slater dang chat cho moi g_i).
%
% Dau ra:
%   ok    - true neu Slater dung (s* < -1e-6)
%   xhat  - diem tim duoc (diem Slater neu ok)
%   sstar - gia tri toi uu cua bai toan phu
%
% Dung solve_nlp.m cua cau 4 (sqp tren Octave, fmincon tren MATLAB).

  x0 = x0(:);
  n = numel(x0);
  if ~isfield(prob, 'g') || isempty(prob.g)
    ok = true;  xhat = x0;  sstar = -1;          % khong co rang buoc bat dang thuc
    return;
  end
  m = numel(prob.g(x0));
  if nargin < 3 || isempty(isAffine), isAffine = false(m, 1); end
  isAffine = logical(isAffine(:));
  w = double(~isAffine);                         % chi rang buoc phi tuyen duoc noi s
  aux.f = @(z) z(end);
  aux.g = @(z) [prob.g(z(1:n)) - w * z(end); -1 - z(end)];
  if isfield(prob, 'h') && ~isempty(prob.h)
    aux.h = @(z) prob.h(z(1:n));
  end
  s0 = max(prob.g(x0)) + 1;                      % diem xuat phat chap nhan duoc cho s
  [z, sstar] = solve_nlp(aux, [x0; s0]);
  xhat = z(1:n);
  ok = sstar < -1e-6;
end
