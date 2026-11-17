function check_solution(c, A, b, f, isOctave)
% CHECK_SOLUTION  Doi chieu gia tri toi uu voi glpk (Octave) hoac linprog (MATLAB).
  n = numel(c);
  if isOctave
    [~, fref] = glpk(c, A, b, zeros(n, 1), [], repmat('S', 1, size(A, 1)), ...
      repmat('C', 1, n), 1);
    name = 'glpk';
  elseif exist('linprog', 'file')
    [~, fref] = linprog(c, [], [], A, b, zeros(n, 1), []);
    name = 'linprog';
  else
    return;
  end
  fprintf('Doi chieu voi %s: f* = %g (sai lech %.2e)\n\n', name, fref, abs(fref - f));
end
