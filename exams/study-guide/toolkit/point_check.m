function point_check(A, b, c, x, sense)
% POINT_CHECK  Kiem tra vector x voi QHTT dang chinh tac Ax = b, x >= 0:
%   phuong an? -> cuc bien (cot ung x_j > 0 doc lap tuyen tinh)? -> toi uu?
  c = c(:);  b = b(:);  x = x(:);
  r = A * x - b;
  fprintf('Ax - b = %s, min x = %g\n', mat2str(r', 6), min(x));
  if norm(r) > 1e-9 || any(x < -1e-12)
    fprintf('=> KHONG PHAI PHUONG AN\n');  return;
  end
  J = find(x > 1e-12)';
  fprintf('J = {j : x_j > 0} = %s, rank A_J = %d\n', mat2str(J), rank(A(:, J)));
  if rank(A(:, J)) < numel(J)
    fprintf('=> LA PHUONG AN nhung KHONG cuc bien, F = %g\n', c' * x);  return;
  end
  % Bo sung cot de co co so (neu x suy bien), roi kiem tra uoc luong
  B = J;
  for j = 1:size(A, 2)
    if numel(B) < size(A, 1) && ~any(B == j) && rank(A(:, [B j])) > rank(A(:, B)), B = [B j]; end
  end
  y = A(:, B)' \ c(B);  D = c' - y' * A;
  fprintf('Co so {%s}: y = %s, Delta_j = %s\n', num2str(B), mat2str(y', 6), mat2str(D + 0, 6));
  if (strcmp(sense, 'max') && all(D <= 1e-12)) || (strcmp(sense, 'min') && all(D >= -1e-12))
    fprintf('=> PHUONG AN CUC BIEN TOI UU, F = %g\n', c' * x);
  else
    fprintf('=> PHUONG AN CUC BIEN (uoc luong chua thoa: chua ket luan toi uu tu co so nay), F = %g\n', c' * x);
  end
end
