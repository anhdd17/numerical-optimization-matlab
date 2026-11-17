function h = draw_arrow(p, v, color, hl)
% DRAW_ARROW  Ve mui ten tu diem p theo vector v (dau mui ten co kich thuoc
%             co dinh, hien thi giong nhau tren MATLAB va GNU Octave).
%             hl: (tuy chon) chieu dai dau mui ten theo don vi cua truc.
  p = p(:)'; v = v(:)';
  q = p + v;
  L = norm(v);
  u = v / L;
  w = [-u(2), u(1)];
  if nargin < 4, hl = min(0.12, 0.35 * L); end
  head = [q; q - hl * u + 0.45 * hl * w; q - hl * u - 0.45 * hl * w];
  h = plot([p(1) q(1) - 0.5 * hl * u(1)], [p(2) q(2) - 0.5 * hl * u(2)], '-', ...
    'Color', color, 'LineWidth', 2);
  fill(head(:, 1), head(:, 2), color, 'EdgeColor', color);
end
