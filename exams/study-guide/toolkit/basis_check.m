function basis_check(A, b, c, B, sense)
% BASIS_CHECK  Kiem tra cot B (chi so) cua QHTT dang chinh tac Ax = b, x >= 0:
%   co so? -> chap nhan duoc? -> toi uu?  sense = 'max' hoac 'min'.
  c = c(:);  b = b(:);  n = size(A, 2);
  AB = A(:, B);
  if size(AB, 1) ~= size(AB, 2) || abs(det(AB)) < 1e-12
    fprintf('det A_B = 0 hoac khong vuong => KHONG PHAI CO SO\n');  return;
  end
  xB = AB \ b;  x = zeros(n, 1);  x(B) = xB;
  y = AB' \ c(B);  D = c' - y' * A;
  fprintf('det A_B = %g\nx_B = %s, x = %s\ny = %s\nDelta_j = c_j - y''A_j = %s\n', ...
    det(AB), mat2str(xB', 6), mat2str(x', 6), mat2str(y', 6), mat2str(D + 0, 6));
  if any(xB < -1e-12)
    fprintf('=> LA CO SO nhung KHONG chap nhan duoc (x_B co phan tu am)\n');
  elseif (strcmp(sense, 'max') && all(D <= 1e-12)) || (strcmp(sense, 'min') && all(D >= -1e-12))
    fprintf('=> CO SO TOI UU, F = %g\n', c' * x);
  else
    fprintf('=> CO SO CHAP NHAN DUOC nhung chua toi uu, F = %g\n', c' * x);
  end
end
