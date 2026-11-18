function x = golden(f, a, b, tol)
% GOLDEN  Lat cat vang tim cuc tieu cua f don mot tren [a, b].
%   tau = 0.618; giu [a, x2] neu f(x1) < f(x2), nguoc lai giu [x1, b].
%   Dung khi b - a <= tol; tra ve trung diem.
  t = (sqrt(5) - 1) / 2;
  x1 = b - t * (b - a);  x2 = a + t * (b - a);  f1 = f(x1);  f2 = f(x2);  k = 0;
  while b - a > tol
    k = k + 1;
    if f1 < f2
      b = x2;  x2 = x1;  f2 = f1;  x1 = b - t * (b - a);  f1 = f(x1);
    else
      a = x1;  x1 = x2;  f1 = f2;  x2 = a + t * (b - a);  f2 = f(x2);
    end
  end
  x = (a + b) / 2;
  fprintf('Lat cat vang: %d buoc, [a, b] = [%.6f, %.6f], x* = %.6f, f = %.6f\n', k, a, b, x, f(x));
end
