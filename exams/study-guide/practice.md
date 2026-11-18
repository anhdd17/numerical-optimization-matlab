# Đề luyện và ngân hàng tự luận

[← Tài liệu ôn thi](README.md)

Đề luyện số 1 soạn theo khuôn [Đề số 1](../test1-paper1/) (8 câu, số liệu từ [`run_practice.m`](run_practice.m)). [Đề luyện số 2](#đề-luyện-số-2) đúng format thi 12 câu + tự luận, phủ các chủ đề còn lại (số liệu từ [`run_practice2.m`](run_practice2.m)). Câu hỏi mẫu cho từng chủ đề trong 18 câu nằm trong các file `topics/`; ngân hàng tự luận ở [cuối trang](#ngân-hàng-tự-luận).

## Đề luyện số 1

*Soạn theo khuôn Đề số 1, thời gian 90 phút.*

**Câu 1.** Dùng phương pháp Gauss–Seidel tìm cực trị địa phương của $`F(x) = 3x_1^2 + 2x_2^2 + 4x_3^2 + 2x_1x_2 - 2x_2x_3 + x_1 - 4x_3`$ với $`x^0 = (1, 1, 1)`$, $`\varepsilon = 10^{-4}`$.

**Câu 2.** Cho bài toán $`\max F = 2x_1 + x_2 + x_3`$ với $`x_1 + x_2 + x_3 = 2`$, $`x_1 - x_2 + x_4 = 4`$, $`x_j \ge 0`$. Ma trận $`B = (A_1, A_2)`$:
A. Là cơ sở tối ưu. B. Là cơ sở chấp nhận được. C. Là cơ sở. D. Không phải là cơ sở.

**Câu 3.** Dùng phương pháp Newton tìm cực tiểu của $`F(x) = x_1^4 + x_1^2 + 2x_2^2 + x_3^2 + x_1x_2 - x_2x_3 - 4x_1 + 2x_3`$ với $`x^0 = (1, 1, 1)`$, $`\varepsilon = 10^{-4}`$.

**Câu 4.** Bài toán $`f(x) = e^{-x} + x^2 - 2\sin 3x \to \min`$ có phương án xấp xỉ ($`\varepsilon = 10^{-4}`$) tối ưu địa phương trong đoạn $`[0, 2]`$ là $`x^* = \ ?`$

**Câu 5.** Cho bài toán $`\max F = x_1 + 2x_2 + x_3 + x_4`$ với $`x_1 + x_2 + x_3 + x_4 = 6`$, $`2x_1 + x_2 + 3x_3 + x_5 = 9`$, $`x_j \ge 0`$. Vector $`x = (1, 2, 1, 2, 2)`$:
A. Là phương án cực biên tối ưu. B. Là phương án cực biên. C. Là phương án. D. Không phải là phương án.

**Câu 6.** Dùng phương pháp Gradient tìm cực tiểu của $`f(x) = (x_1 - 1)^2 + 10(x_2 - x_1^2)^2`$ với $`x^0 = (0, 0)`$, $`\varepsilon = 10^{-3}`$.

**Câu 7.** Dùng phương pháp Gradient quán tính tìm cực tiểu của $`f(x) = 3x_1^2 + 2x_2^2 + x_3^2 + x_1x_2 - x_2x_3 + 2x_1 - x_2 + 3x_3`$ với $`x^0 = (1, 1, 1)`$, $`\varepsilon = 10^{-4}`$.

**Câu 8.** Dùng phương pháp tựa Newton tìm cực tiểu của $`F(x, y, z) = x^4 + y^4 + 2x^2 + 3y^2 + z^2 - xy + yz - 2x + y - 3z`$ với $`x^0 = (1, 1, 1)`$, $`\varepsilon = 10^{-4}`$.

**Phần B (tự luận).** Cho $`a`$ là tháng sinh của học viên. Xét bài toán (P):

```math
\min f(x) = \Big(x_1 - \frac a4\Big)^2 + \sum_{i=2}^n x_i^2\quad \text{s.t.}\quad x_1 + \dots + x_n = 1,\quad x_i \ge 0 .
```

a) (P) có phải bài toán tối ưu lồi không? Viết bài toán đối ngẫu và điều kiện KKT.
b) Giải (P) với $`n = 2`$ (theo $`a`$).

---

## Đáp án và lời giải đề luyện số 1

*Làm xong đề rồi mới đọc phần này.* Số liệu lấy từ output của `run_practice.m`.

| Câu | Đáp án |
|---|---|
| 1 | $`x^* \approx (-0.3235,\ 0.4706,\ 0.6176) = (-11/34,\ 8/17,\ 21/34)`$, $`F^* = -95/68 \approx -1.39706`$ |
| 2 | **C** (là cơ sở, nhưng $`x_B = (3, -1)`$ có phần tử âm) |
| 3 | $`x^* \approx (0.8846,\ -0.5385,\ -1.2692)`$ |
| 4 | $`x^* \approx 0.5015`$, $`f(x^*) \approx -1.13848`$ |
| 5 | **C** (là phương án, không cực biên) |
| 6 | $`x^* = (1, 1)`$ (phương pháp cho $`(0.9990,\ 0.9980)`$ sau 441 bước) |
| 7 | $`x^* = (-0.325,\ -0.05,\ -1.525) = (-13/40,\ -1/20,\ -61/40)`$ |
| 8 | $`x^* \approx (0.3632,\ -0.3558,\ 1.6779)`$, $`F^* \approx -3.091285`$ |
| B | $`a \le 4`$: $`x^* = \big(\frac{4 + a}{8}, \frac{4 - a}{8}\big)`$; $`a \ge 4`$: $`x^* = (1, 0)`$ |

**Câu 1.** Lập ma trận: $`H = [6\ 2\ 0;\ 2\ 4\ -2;\ 0\ -2\ 8]`$, $`q = (1, 0, -4)`$. Các định thức con chính $`6, 20, 136`$ đều dương. Công thức Gauss–Seidel:

```math
x_1 = -\frac{1 + 2x_2}{6},\qquad x_2 = \frac{x_3 - x_1}{2},\qquad x_3 = \frac{x_2 + 2}{4}.
```

Vòng 1: $`x_1 = -3/6 = -0.5`$, $`x_2 = (1 + 0.5)/2 = 0.75`$, $`x_3 = 2.75/4 = 0.6875`$. Sau 9 vòng, $`\Vert x_9 - x_8 \Vert = 5.5 \cdot 10^{-5}`$ và $`x_9 = (-0.323546, 0.470603, 0.617651)`$. Nghiệm đúng $`-H^{-1}q = (-11/34, 8/17, 21/34)`$.

**Câu 2.** $`A_B = [1\ 1;\ 1\ -1]`$ có $`\det = -2 \ne 0`$, nên là cơ sở. Giải $`A_B x_B = b`$: $`x_1 + x_2 = 2`$, $`x_1 - x_2 = 4`$ cho $`x_B = (3, -1)`$, có phần tử âm, nên không chấp nhận được. **Đáp án C.** Thêm: cơ sở $`(A_1, A_4)`$ cho $`x = (2, 0, 0, 2)`$, $`\Delta = (0, -1, -1, 0) \le 0`$, là cơ sở tối ưu với $`F = 4`$.

**Câu 3.** Phần hằng của Hessian là $`[2\ 1\ 0;\ 1\ 4\ -1;\ 0\ -1\ 2]`$ (định thức con chính $`2, 7, 12`$), cộng thêm $`12x_1^2 \ge 0`$, nên $`F`$ lồi chặt. Khử biến: từ $`\partial F/\partial x_3 = 2x_3 - x_2 + 2 = 0`$ và $`\partial F/\partial x_2 = x_1 + 4x_2 - x_3 = 0`$ được $`x_3 = (x_2 - 2)/2`$ và $`x_2 = -(x_1 + 1)/3.5`$. Thế vào $`\partial F/\partial x_1 = 0`$ được phương trình một ẩn $`4x_1^3 + \tfrac{12}{7}x_1 - \tfrac{30}{7} = 0`$, có nghiệm $`x_1 = 0.884636`$. Newton hội tụ trong 3 bước ($`\Vert g \Vert`$: $`5.83 \to 0.126 \to 1.3\cdot 10^{-3} \to 1.5 \cdot 10^{-7}`$).

**Câu 4.** $`f'(x) = -e^{-x} + 2x - 6\cos 3x`$, với $`f'(0) = -7`$ và $`f'(2) = -1.90`$. Các điểm dừng: $`0.501505`$ (cực tiểu, $`f = -1.138483`$) và $`1.769501`$ (cực đại, $`f = 4.9566`$). Vì $`f'(2) \lt 0`$ nên $`x = 2`$ là cực tiểu địa phương tại biên, nhưng $`f(2) = 4.694`$ lớn hơn nhiều. Lát cắt vàng trên $`[0, 2]`$ cần 21 bước ($`k \ge \ln(5\cdot 10^{-5})/\ln\tau = 20.6`$), cho $`x^* = 0.501502`$; `fminbnd` cho $`0.501504`$.

![Câu 4 đề luyện](figures/practice_q4.png)

**Câu 5.** $`1 + 2 + 1 + 2 = 6`$ và $`2 + 2 + 3 + 2 = 9`$, $`x \ge 0`$, nên là phương án. Nó có 5 thành phần dương, nhiều hơn $`m = 2`$, nên các cột phụ thuộc tuyến tính và $`x`$ không cực biên. **Đáp án C** ($`F = 8`$). Thêm: phương án cực biên tối ưu là $`(0, 6, 0, 0, 3)`$ với $`F = 12`$.

**Câu 6.** $`f \ge 0`$, và $`f = 0`$ khi và chỉ khi $`x_1 = 1`$, $`x_2 = x_1^2 = 1`$. Gradient với tìm kiếm chính xác zic-zac theo hai trục (các bước đầu $`\alpha = 0.162,\ 0.05,\ 0.070,\ 0.05, \dots`$), dừng sau 441 bước tại $`(0.999012, 0.997975)`$ với $`\Vert g \Vert = 9.9\cdot 10^{-4}`$. Với $`\varepsilon = 10^{-3}`$ trên gradient, sai số theo $`x`$ vẫn khoảng $`2 \cdot 10^{-3}`$, nên ghi đáp án là $`(1, 1)`$ kèm giá trị xấp xỉ.

**Câu 7.** $`H = [6\ 1\ 0;\ 1\ 4\ -1;\ 0\ -1\ 2]`$, $`q = (2, -1, 3)`$, với trị riêng $`1.5505,\ 4,\ 6.4495`$. Tham số $`\alpha = 4/(\sqrt L + \sqrt\mu)^2 = 0.279241`$, $`\beta = 0.116963`$. Phương pháp quán tính dừng sau 14 bước tại $`(-0.324993, -0.049993, -1.52499)`$. Nghiệm đúng $`-H^{-1} q = (-13/40, -1/20, -61/40)`$.

**Câu 8.** Phần bậc hai có ma trận $`[4\ -1\ 0;\ -1\ 6\ 1;\ 0\ 1\ 2] \succ 0`$ (định thức con chính $`4, 23, 42`$), nên $`F`$ lồi chặt. Rút $`z = (3 - y)/2`$ từ $`\partial F/\partial z = 0`$. BFGS (Armijo, $`B_0 = I`$) dừng sau 7 bước tại $`(0.363156, -0.355765, 1.67788)`$; `fminsearch` cho $`(0.363162, -0.355767, 1.67788)`$, $`F^* = -3.0912851`$.

### Phần B

**a)** $`\nabla^2 f = 2I \succ 0`$, nên $`f`$ lồi chặt; ràng buộc affine; vậy (P) là bài toán lồi. Điểm $`x = (1/n, \dots, 1/n) \gt 0`$ thỏa Slater, nên có đối ngẫu mạnh và KKT là điều kiện cần và đủ. Đặt $`c = (a/4, 0, \dots, 0)`$ và $`L = \sum_i (x_i - c_i)^2 + \nu(1 - \sum_i x_i) - \sum_i \lambda_i x_i`$. Cực tiểu theo từng $`x_i \in \mathbb{R}`$ của $`(x_i - c_i)^2 - \mu_i x_i`$ đạt tại $`x_i = c_i + \mu_i/2`$, giá trị $`-\mu_i c_i - \mu_i^2/4`$. Do đó

```math
g(\nu, \lambda) = \nu - \frac a4 (\nu + \lambda_1) - \frac14 \sum_{i=1}^n (\nu + \lambda_i)^2,\qquad (D):\ \max_{\nu \in \mathbb{R},\ \lambda \ge 0} g(\nu, \lambda).
```

KKT: $`2(x_1 - a/4) - \nu - \lambda_1 = 0`$; $`2x_i - \nu - \lambda_i = 0`$ ($`i \ge 2`$); $`\sum x_i = 1`$; $`x \ge 0`$; $`\lambda \ge 0`$; $`\lambda_i x_i = 0`$.

**b)** Với $`n = 2`$, xét ba trường hợp:

| Trường hợp | Từ KKT | Điều kiện | Kết quả |
|---|---|---|---|
| $`x_1, x_2 \gt 0`$ ($`\lambda = 0`$) | $`x_1 = \frac a4 + \frac\nu2`$, $`x_2 = \frac\nu2`$, tổng bằng 1 nên $`\nu = 1 - \frac a4`$ | $`x_2 = \frac{4 - a}{8} \gt 0 \iff a \lt 4`$ | $`x^* = \big(\frac{4+a}{8}, \frac{4-a}{8}\big)`$, $`f^* = \frac{(4-a)^2}{32}`$ |
| $`x_2 = 0`$, $`x_1 = 1`$ | $`\lambda_1 = 0`$, $`\nu = 2 - \frac a2`$, $`\lambda_2 = -\nu = \frac a2 - 2`$ | $`\lambda_2 \ge 0 \iff a \ge 4`$ | $`x^* = (1, 0)`$, $`f^* = \big(1 - \frac a4\big)^2`$ |
| $`x_1 = 0`$, $`x_2 = 1`$ | $`\nu = 2`$, $`\lambda_1 = -\frac a2 - 2`$ | $`\lambda_1 \lt 0`$ | loại |

Như vậy với tháng sinh 1–3 nghiệm nằm trong, với tháng 4 thì $`x^* = (1, 0)`$ và $`f^* = 0`$, với tháng 5–12 nghiệm nằm trên biên với $`\lambda_2 = a/2 - 2 \gt 0`$.
Kiểm chứng bằng chương trình: $`a = 2`$ cho $`x^* = (0.75, 0.25)`$, $`f^* = 0.125 = g(\nu^*, \lambda^*)`$; $`a = 8`$ cho $`x^* = (1, 0)`$, $`\nu^* = -2`$, $`\lambda^* = (0, 2)`$, $`f^* = 1 = g`$. `sqp` cho cùng nghiệm, kể cả với $`n = 4`$ (khi $`a = 2`$: $`x^* = (0.625, 0.125, 0.125, 0.125)`$).

---

## Đề luyện số 2

*Đúng format thi: 12 câu trắc nghiệm/điền đáp án (6 điểm) và tự luận (4 điểm), 90 phút. Đề này dùng các chủ đề mà Đề số 1 và Đề luyện số 1 chưa ra: câu 1–6, 9, 13–18 của ngân hàng.*

**Câu 1.** Dùng thuật toán đơn hình giải $`\max F = 4x_1 + 3x_2`$ với $`2x_1 + x_2 \le 10`$, $`x_1 + x_2 \le 7`$, $`x_2 \le 5`$, $`x \ge 0`$ (cơ sở xuất phát gồm các biến bù, quy tắc Dantzig). Phương án tối ưu $`x^* = \ ?`$, $`F^* = \ ?`$, số bước xoay $`= \ ?`$

**Câu 2.** Áp dụng phương pháp hai pha cho $`\min F = x_1 - 2x_2`$ với $`x_1 + x_2 \ge 2`$, $`-x_1 + x_2 \le 1`$, $`x \ge 0`$. Kết luận:
A. Pha 1 cho $`w^* \gt 0`$ nên bài toán không có phương án. B. $`w^* = 0`$; pha 2 kết luận hàm mục tiêu không bị chặn dưới. C. $`w^* = 0`$; nghiệm tối ưu $`x^* = (0, 2)`$. D. $`w^* = 0`$; nghiệm tối ưu $`x^* = (1/2, 3/2)`$.

**Câu 3.** Giải bằng bài toán M: $`\min F = 2x_1 + 3x_2 + x_3`$ với $`x_1 + x_2 + x_3 = 5`$, $`2x_1 + x_2 \ge 4`$, $`x \ge 0`$. Số biến giả cần thêm $`= \ ?`$, $`x^* = \ ?`$, $`F^* = \ ?`$

**Câu 4.** Bài toán $`\min (x_1 - 2)^2 + (x_2 - 1)^2`$ với $`x_1^2 - x_2 \le 0`$, $`x_1 + x_2 \le 2`$. Điểm $`\bar x = (1, 1)`$:
A. Là điểm KKT với $`\lambda = (2/3, 2/3)`$ và là nghiệm tối ưu toàn cục. B. Là điểm KKT nhưng không phải cực tiểu. C. Không phải điểm KKT vì có nhân tử âm. D. Không thỏa LICQ nên không kết luận được.

**Câu 5.** Bài toán $`\min x_1^2 + 2x_2^2`$ với $`x_1 + x_2 \ge 3`$. Hàm đối ngẫu Lagrange $`g(\lambda) = \ ?`$, nhân tử tối ưu $`\lambda^* = \ ?`$, giá trị đối ngẫu $`d^* = \ ?`$

**Câu 6.** Dùng chia đôi theo đạo hàm tìm cực tiểu của $`f(x) = x^4 - 3x^2 + x`$ trên $`[0, 2]`$. Khoảng chứa cực tiểu sau 3 bước là $`?`$, và nghiệm với $`\varepsilon = 10^{-4}`$ là $`x^* = \ ?`$

**Câu 7.** Dùng phương pháp Hooke–Jeeves ($`\Delta = 0.5`$, giảm bước một nửa) tìm cực tiểu của $`F(x) = 2x_1^2 + x_1x_2 + x_2^2 - 3x_1 - x_2`$ từ $`x^0 = (0, 0)`$, $`\varepsilon = 10^{-6}`$. Điểm sau bước dò đầu tiên $`= \ ?`$, nghiệm $`x^* = \ ?`$

**Câu 8.** Bài toán $`\min -\sum_{i=1}^3 \ln(\alpha_i + x_i)`$ với $`\sum x_i = 1`$, $`x \ge 0`$, $`\alpha = (0.5, 1, 2)`$. Nghiệm $`x^* = \ ?`$, mức nước $`1/\nu^* = \ ?`$

**Câu 9.** Dùng Frank–Wolfe (tìm kiếm chính xác) cho $`\min (x_1 - 2)^2 + (x_2 - 2)^2`$ với $`x_1 + 2x_2 \le 4`$, $`x \ge 0`$, $`x_0 = (0, 0)`$. Điểm $`x_1, x_2 = \ ?`$ và nghiệm $`x^* = \ ?`$

**Câu 10.** Dùng phương pháp chiếu gradient Rosen cho $`\min (x_1 - 2)^2 + (x_2 - 2)^2`$ với $`x_1 + x_2 \le 2`$ (1), $`-x_1 \le 0`$ (2), $`-x_2 \le 0`$ (3), $`x_0 = (2, 0)`$. Ước lượng nhân tử $`u`$ tại $`x_0`$, ràng buộc bị bỏ, hướng $`d_0`$ và điểm $`x_1`$ là gì?

**Câu 11.** Dùng phương pháp Newton cho $`\min e^{x_1} + e^{2x_2} + x_3^2`$ với $`x_1 + x_2 + x_3 = 1`$, $`x_0 = (1, 0, 0)`$. Giá trị $`\lambda(x_0)^2/2 = \ ?`$, nghiệm $`x^* = \ ?`$

**Câu 12.** Dùng đơn hình Beale cho $`\min f = x_1^2 + 2x_2^2 - 2x_1 - 6x_2`$ với $`x_1 + x_2 + x_3 = 2`$, $`x \ge 0`$, cơ sở xuất phát $`\{x_3\}`$. Ở bước 1, biến vào, $`\theta_1`$, $`\theta_2`$ là gì? Nghiệm $`x^* = \ ?`$, $`f^* = \ ?`$

**Phần B (tự luận).** Cho $`a`$ là tháng sinh của học viên. Xét bài toán (P):

```math
\min f_0(x) = x_1^2 + x_2^2\quad \text{s.t.}\quad x_1 + x_2 \ge a .
```

a) Chứng minh (P) là bài toán lồi và thỏa Slater. Giải (P) bằng KKT.
b) Viết hàm đối ngẫu, bài toán đối ngẫu và kiểm tra đối ngẫu mạnh.
c) Phương pháp hàm chắn: viết bài toán định tâm với tham số $`t \gt 0`$, tìm đường trung tâm $`x^*(t)`$ dạng đóng, nhân tử ước lượng $`\lambda(t)`$, và chứng minh khe đối ngẫu bằng $`1/t`$.
d) Với $`t^{(0)} = 1`$, $`\mu = 10`$, cần bao nhiêu lần định tâm để khe $`\le 10^{-6}`$?

---

## Đáp án và lời giải đề luyện số 2

Mọi số liệu lấy từ output của [`run_practice2.m`](run_practice2.m). Script chạy phương pháp đề yêu cầu bằng các hàm của bài tập lớn, giải tay song song và đối chiếu với solver.

| Câu | Đáp án |
|---|---|
| 1 | $`x^* = (3, 4)`$, $`F^* = 24`$, 2 bước xoay |
| 2 | **B** |
| 3 | 1 biến giả; $`x^* = (2, 0, 3)`$, $`F^* = 7`$ |
| 4 | **A** |
| 5 | $`g(\lambda) = 3\lambda - \tfrac38\lambda^2`$, $`\lambda^* = 4`$, $`d^* = 6`$ |
| 6 | $`[1, 1.25]`$; $`x^* = 1.1309`$ |
| 7 | $`(0.5, 0)`$; $`x^* = (5/7, 1/7) = (0.7143, 0.1429)`$ |
| 8 | $`x^* = (0.75, 0.25, 0)`$, mức nước $`1.25`$ ($`\nu^* = 0.8`$) |
| 9 | $`x_1 = (2, 0)`$, $`x_2 = (1, 1)`$; $`x^* = (1.6, 1.2)`$, $`f^* = 0.8`$ |
| 10 | $`u = (0, -4)`$ cho $`W = \{1, 3\}`$, bỏ (3); $`d_0 = (-2, 2)`$; $`x_1 = (1, 1)`$ là nghiệm, $`\lambda_1 = 2`$ |
| 11 | $`\lambda(x_0)^2/2 = 0.8528`$; $`x^* = (0.4003, -0.1464, 0.7461)`$ |
| 12 | $`x_2`$ vào, $`\theta_1 = 2`$, $`\theta_2 = 1.5`$ (tạo biến tự do); $`x^* = (2/3, 4/3, 0)`$, $`f^* = -16/3`$ |

**Câu 1.** Dạng chính tắc thêm $`x_3, x_4, x_5`$. Bước 1: $`x_1`$ vào ($`\Delta = 4`$ lớn nhất); tỉ số $`10/2 = 5 \lt 7/1`$, nên $`x_3`$ ra. Bước 2: $`x_2`$ vào, $`x_4`$ ra. Bảng cuối có mọi $`\Delta_j \le 0`$: $`x^* = (3, 4, 0, 0, 1)`$, tức $`x_5 = 1`$ (ràng buộc $`x_2 \le 5`$ không chặt). Kiểm tra bằng đỉnh: $`(5, 0) \to 20`$, $`(3, 4) \to 24`$, $`(2, 5) \to 23`$.

**Câu 2.** Dòng 1 có biến bù hệ số $`-1`$ nên cần 1 biến giả (dòng 2 có sẵn cột đơn vị $`x_4`$). Pha 1 cho $`w^* = 0`$. Pha 2 có cột với $`\Delta_j \gt 0`$ mà mọi phần tử trên cột $`\le 0`$, nên `exitflag = -3`. Lý do trực tiếp: trên tia $`x = (s, s + 1)`$ chấp nhận được, $`F = -s - 2 \to -\infty`$. Đáp án C và D là các đỉnh, nhưng bài toán không có nghiệm tối ưu.

**Câu 3.** Bẫy: cột $`x_3 = (1, 0)^T`$ đã là cột đơn vị của dòng 1, nên chỉ dòng 2 (có $`-x_4`$) cần biến giả. Như vậy là **1** biến giả, không phải 2. Thay $`x_3 = 5 - x_1 - x_2`$: $`F = x_1 + 2x_2 + 5`$ với $`2x_1 + x_2 \ge 4`$, $`x_1 + x_2 \le 5`$. Cực tiểu tại $`(2, 0)`$, nên $`x^* = (2, 0, 3)`$, $`F^* = 7`$. Kết quả trùng với `glpk`.

**Câu 4.** Cả hai ràng buộc tích cực tại $`(1, 1)`$; $`\nabla g_1 = (2, -1)`$ và $`\nabla g_2 = (1, 1)`$ độc lập tuyến tính (LICQ). Giải $`\lambda_1 (2, -1) + \lambda_2 (1, 1) = -\nabla f = (2, 0)`$ được $`\lambda_1 = \lambda_2 = 2/3 \ge 0`$. Bài toán lồi ($`f`$ lồi, $`g_1, g_2`$ lồi), nên điểm KKT là nghiệm toàn cục: $`f^* = 1`$. `sqp` cho cùng kết quả.

**Câu 5.** $`L = x_1^2 + 2x_2^2 + \lambda(3 - x_1 - x_2)`$, cực tiểu tại $`x_1 = \lambda/2`$, $`x_2 = \lambda/4`$, nên $`g(\lambda) = 3\lambda - \tfrac38\lambda^2`$. Cực đại tại $`\lambda^* = 4`$, $`d^* = 6`$. Nghiệm gốc $`x^* = (2, 1)`$, $`f^* = 6 = d^*`$ (Slater thỏa).

**Câu 6.** $`f'(x) = 4x^3 - 6x + 1`$. Bước 1: $`f'(1) = -1 \lt 0 \Rightarrow [1, 2]`$. Bước 2: $`f'(1.5) = 5.5 \gt 0 \Rightarrow [1, 1.5]`$. Bước 3: $`f'(1.25) = 1.3125 \gt 0 \Rightarrow [1, 1.25]`$. Chạy đến $`10^{-4}`$ được $`x^* = 1.1309`$ (15 bước).
**Bẫy:** $`f'`$ đổi dấu **hai lần** trên $`[0, 2]`$: có cực đại địa phương tại $`0.1699`$ ($`f = 0.0841`$) và cực tiểu tại $`1.1309`$ ($`f = -1.0702`$). Chia đôi chỉ đúng vì $`f'(1) \lt 0`$ đã loại $`[0, 1]`$. Phải so với $`f(0) = 0`$ và $`f(2) = 6`$ để khẳng định đây là cực tiểu toàn cục trên đoạn.

**Câu 7.** Bước dò quanh $`(0, 0)`$ với $`F = 0`$: thử $`x_1 + 0.5`$ được $`F(0.5, 0) = -1`$, nhận. Thử $`x_2 \pm 0.5`$ được $`F(0.5, 0.5) = -1`$ (không giảm thật sự) và $`F(0.5, -0.5) = -0.5`$, nên giữ $`(0.5, 0)`$. $`F`$ lồi chặt ($`H = [4\ 1;\ 1\ 2]`$), nên nghiệm là $`-H^{-1}q = (5/7, 1/7)`$, $`F^* = -8/7`$. `hooke_jeeves` cho $`(0.714287, 0.142857)`$ sau 181 lần tính $`F`$. Nghiệm có mẫu số 7, không nằm trên lưới $`0.5 \cdot 2^{-j}`$, nên phương pháp không thể "may mắn" rơi trúng nghiệm.

**Câu 8.** Water-filling: $`x_i = \max(0, L - \alpha_i)`$ với mức nước $`L = 1/\nu`$. Thử hai kênh đầu: $`(L - 0.5) + (L - 1) = 1 \Rightarrow L = 1.25 \lt \alpha_3 = 2`$, hợp lệ. Vậy $`x^* = (0.75, 0.25, 0)`$, $`\nu^* = 0.8`$, $`\lambda_3 = \nu^* - 1/\alpha_3 = 0.3`$, $`f^* = -1.139434`$. `water_filling.m` cho cùng kết quả.

**Câu 9.** Các đỉnh của miền: $`(0, 0), (4, 0), (0, 2)`$.
- Bước 1: $`\nabla f(0) = (-4, -4)`$ chọn $`y = (4, 0)`$; trên đoạn $`f(4\gamma, 0)`$ đạt cực tiểu tại $`\gamma = 1/2`$, nên $`x_1 = (2, 0)`$.
- Bước 2: $`\nabla f = (0, -4)`$ chọn $`y = (0, 2)`$, $`\gamma = 1/2`$, nên $`x_2 = (1, 1)`$ với $`f = 2`$.
- Nghiệm là hình chiếu của $`(2, 2)`$ lên $`x_1 + 2x_2 = 4`$: $`x^* = (2, 2) - 0.4\,(1, 2) = (1.6, 1.2)`$, $`f^* = 0.8`$.

**Nhận xét quan trọng:** Frank–Wolfe chuẩn zic-zac giữa hai đỉnh và hội tụ chỉ với tốc độ $`O(1/k)`$. Đến khe $`\le 10^{-2}, 10^{-3}, 10^{-4}`$ cần lần lượt 753, 7657, 76769 bước. Biến thể away-step chỉ cần 5 bước. Trong phòng thi, hãy ghi nghiệm đúng (tính qua KKT hoặc hình chiếu) kèm số bước mà chương trình đã chạy.

**Câu 10.** Tại $`x_0 = (2, 0)`$ có $`W = \{1, 3\}`$, $`N = [1\ 1;\ 0\ -1]`$, $`\nabla f(x_0) = (0, -4)`$. Vì $`N`$ vuông khả nghịch nên $`P = 0`$, $`d = 0`$. Ước lượng nhân tử $`u = -(NN^T)^{-1}N\nabla f = (0, -4)`$, trong đó $`u_3 = -4 \lt 0`$, nên **bỏ (3)**.
Với $`W = \{1\}`$: $`P = I - \tfrac12 [1\ 1;\ 1\ 1] = \tfrac12[1\ -1;\ -1\ 1]`$, $`d_0 = -P\nabla f = (-2, 2)`$. Tìm kiếm chính xác trên $`(2 - 2\alpha, 2\alpha)`$ cho $`\alpha = 0.5`$, nên $`x_1 = (1, 1)`$.
Tại $`x_1`$: $`d = 0`$, $`u_1 = 2 \ge 0`$, dừng. Nghiệm $`x^* = (1, 1)`$, $`f^* = 2`$, $`\lambda_1 = 2`$.

**Câu 11.** Tại $`x_0`$: $`g = (e, 2, 0)`$, $`H = \mathrm{diag}(e, 4, 2)`$. Giải hệ KKT $`[H\ A^T;\ A\ 0][\Delta x; w] = [-g; 0]`$ được $`\Delta x = (-0.5064, -0.1645, 0.6709)`$, $`\lambda^2/2 = \Delta x^T H \Delta x/2 = 0.8528`$.
Newton chấp nhận được hội tụ sau 4 bước: $`x^* = (0.400294, -0.146426, 0.746132)`$, $`\nu^* = -1.492264`$, $`f^* = 2.795109`$.
Kiểm tra tay: điều kiện KKT $`e^{x_1} = 2e^{2x_2} = 2x_3 = s`$ với $`s = -\nu`$, nên $`x_1 = \ln s`$, $`x_2 = \tfrac12\ln(s/2)`$, $`x_3 = s/2`$. Từ $`x_1 + x_2 + x_3 = 1`$ được phương trình một biến của $`s`$, cho $`s = 1.4923`$.

**Câu 12.** Tại cơ sở $`\{x_3\}`$: $`\partial f/\partial x_1 = -2`$, $`\partial f/\partial x_2 = -6`$, nên $`x_2`$ vào. $`\theta_1 = 2`$ ($`x_3 = 2 - x_2`$ về 0) và $`\theta_2 = 1.5`$ ($`4x_2 - 6 = 0`$). Vì $`\theta_2 \lt \theta_1`$, tạo biến tự do $`u_1`$; bảng in của chương trình:

| k | $`x_k`$ | $`f`$ | vào | ra | $`\theta_1`$ | $`\theta_2`$ | bước |
|---|---|---|---|---|---|---|---|
| 0 | (0, 0, 2) | 0 | $`x_2`$ | – | 2 | 1.5 | tạo $`u_1`$ |
| 1 | (0, 1.5, 0.5) | −4.5 | $`x_1`$ | $`x_3`$ | 0.5 | 1 | xoay |
| 2 | (0.5, 1.5, 0) | −5.25 | $`u_1`$ | – | 3 | 0.3333 | tạo $`u_2`$, loại $`u_1`$ |
| 3 | (2/3, 4/3, 0) | −16/3 | | | | | tối ưu |

Kiểm tra bằng KKT: $`2x_1 - 2 + \lambda = 0`$, $`4x_2 - 6 + \lambda = 0`$, $`x_1 + x_2 = 2`$, cho $`\lambda = 2/3`$, $`x^* = (2/3, 4/3)`$. `sqp` cho cùng nghiệm.

### Phần B, đề luyện số 2

**a)** $`f_0`$ có $`\nabla^2 f_0 = 2I \succ 0`$ nên lồi chặt; ràng buộc $`a - x_1 - x_2 \le 0`$ là affine. Điểm $`(a, a)`$ có $`a - 2a = -a \lt 0`$, nên Slater thỏa; KKT là điều kiện cần và đủ, nghiệm duy nhất. Hệ KKT là $`2x_i - \lambda = 0`$, $`\lambda \ge 0`$, $`\lambda(a - x_1 - x_2) = 0`$. Nếu $`\lambda = 0`$ thì $`x = 0`$, vi phạm ràng buộc (vì $`a \ge 1`$). Vậy ràng buộc tích cực: $`x_1 = x_2 = a/2`$, $`\lambda^* = a`$, $`f^* = a^2/2`$.

**b)** $`L = x_1^2 + x_2^2 + \lambda(a - x_1 - x_2)`$ đạt cực tiểu tại $`x_i = \lambda/2`$:

```math
g(\lambda) = a\lambda - \frac{\lambda^2}{2},\qquad (D):\ \max_{\lambda \ge 0} g(\lambda) \;\Rightarrow\; \lambda^* = a,\quad d^* = \frac{a^2}{2} = f^* .
```

**c)** Bài toán định tâm (một ràng buộc, $`m = 1`$):

```math
\min_x\ t\,(x_1^2 + x_2^2) - \ln(x_1 + x_2 - a) .
```

Điều kiện dừng là $`2t x_i = 1/(x_1 + x_2 - a)`$. Hai phương trình như nhau, nên $`x_1 = x_2 = s`$, và $`2ts(2s - a) = 1`$, tức $`4ts^2 - 2tas - 1 = 0`$. Lấy nghiệm có $`2s \gt a`$:

```math
x^*(t) = s(t)\,(1, 1),\qquad s(t) = \frac{ta + \sqrt{t^2a^2 + 4t}}{4t} \;\xrightarrow{t \to \infty}\; \frac a2 .
```

Nhân tử ước lượng là $`\lambda(t) = -1/(t f_1(x^*(t))) = 1/(t(2s - a)) = 2s(t)`$. Từ đó khe đối ngẫu:

```math
f_0(x^*(t)) - g(\lambda(t)) = 2s^2 - \big(2as - 2s^2\big) = 2s(2s - a) = \frac1t = \frac mt .
```

Kiểm chứng số với $`a = 2`$: $`t = 1`$ cho $`x^*(t) = (1.207107, 1.207107)`$ và $`\lambda = 2.414214`$; $`t = 100`$ cho $`(1.002494, 1.002494)`$. Khe bằng $`1/t`$ đến sai số $`10^{-15}`$ (cả với $`a = 7`$).

**d)** Cần $`m/t \le 10^{-6}`$ với $`t = 10^k`$, nên $`k = \lceil \log_{10}(1/10^{-6}) \rceil = 6`$ lần tăng $`t`$. Cộng lần định tâm đầu tại $`t^{(0)} = 1`$ là **7 lần định tâm**. `barrier_method` cũng báo 7 bước ngoài và cho $`x = (a/2, a/2)`$ với $`a = 2, 7`$.

---

## Ngân hàng tự luận

Phần tự luận chiếm 4 điểm. Dưới đây là các bài có lời giải đầy đủ (dạng 1) và các dàn ý "trình bày thuật toán" (dạng 2); khung trả lời chung ở [mục 5 của tài liệu](README.md#5-phần-tự-luận-4-điểm).

**Dạng 1: quy hoạch lồi, đối ngẫu, KKT (lời giải đầy đủ)**

| Bài | Nội dung | Ở đâu |
|---|---|---|
| Đề số 1, Phần B | $`\min \sum x_i \ln(a + x_i)`$ trên đơn hình: lồi, đối ngẫu (Lambert W), KKT, $`n = 2`$ | [test1-paper1](../test1-paper1/#phần-b) |
| Đề luyện số 1, Phần B | $`\min (x_1 - a/4)^2 + \sum_{i \ge 2} x_i^2`$ trên đơn hình: chia trường hợp theo $`a`$ | [trên trang này](#phần-b) |
| Đề luyện số 2, Phần B | $`\min x_1^2 + x_2^2`$, $`x_1 + x_2 \ge a`$: KKT, đối ngẫu, đường trung tâm của hàm chắn, khe $`1/t`$ | [trên trang này](#phần-b-đề-luyện-số-2) |
| Tự luận 4.1 | Giải bài toán bằng KKT theo các tập tích cực | [kkt-duality-convex](topics/kkt-duality-convex.md) |
| Tự luận 4.2 | Các điều kiện chính quy; ví dụ cực tiểu không phải điểm KKT | [kkt-duality-convex](topics/kkt-duality-convex.md) |
| Tự luận 5.1 | Đối ngẫu của một QP, kiểm tra đối ngẫu mạnh | [kkt-duality-convex](topics/kkt-duality-convex.md) |
| Tự luận 5.2 | Đối ngẫu LP, khôi phục nghiệm gốc bằng điều kiện bù | [kkt-duality-convex](topics/kkt-duality-convex.md) |
| Tự luận 13.1 | Water-filling: lồi, Slater, KKT, đối ngẫu, giải | [kkt-duality-convex](topics/kkt-duality-convex.md) |
| Tự luận 13.2 | Hàm bậc hai tách biến trên đơn hình (xét trường hợp) | [kkt-duality-convex](topics/kkt-duality-convex.md) |
| Tự luận 13.3 | Chứng minh KKT đủ cho bài toán lồi; khi nào KKT cần | [kkt-duality-convex](topics/kkt-duality-convex.md) |

**Dạng 2: trình bày thuật toán và ví dụ (dàn ý ở mục "Nếu ra tự luận" của từng chủ đề)**

| Câu | Chủ đề | File |
|---|---|---|
| 1, 2, 3, 18 | Đơn hình, hai pha, bài toán M, Beale | [lp-and-qp](topics/lp-and-qp.md) |
| 6, 7 | Chia đôi, lát cắt vàng | [line-search](topics/line-search.md) |
| 8–12 | Gauss–Seidel, Hooke–Jeeves, gradient/quán tính, Newton, tựa Newton | [unconstrained](topics/unconstrained.md) |
| 14–17 | Frank–Wolfe, Rosen/Murtagh–Sargent, Newton có ràng buộc, hàm chắn | [constrained-methods](topics/constrained-methods.md) |
