# Các phương pháp cho quy hoạch lồi có ràng buộc (câu 14, 15, 16, 17)

[← Tài liệu ôn thi](../README.md)

Bốn câu này cùng giải một loại bài toán, là quy hoạch lồi có ràng buộc. Vì vậy với câu điền đáp án dạng "dùng phương pháp X giải … cho nghiệm $`x^* =`$", **đáp án chính là nghiệm tối ưu của bài toán lồi**, không phụ thuộc vào phương pháp. Cách làm nhanh gồm hai bước:

1. Tìm nghiệm bằng KKT giải tay, hoặc bằng `solve_qp` / `solve_nlp` (các hàm của câu 4, 5 trong repo; trên Octave chúng gọi `qp` / `sqp`).
2. Chạy hàm của đúng phương pháp đề yêu cầu (`frank_wolfe`, `gradient_projection`, `newton_eq`, `barrier_method`) với `display = true` để chép bảng lặp vào bài.

Mọi số liệu dưới đây lấy từ output của [`check_constrained.m`](../check_constrained.m). Chạy bằng lệnh `cd exams/study-guide; check_constrained` (lệnh này tự gọi `exam_setup` để nạp đường dẫn).

| Câu | Bài toán | Mỗi bước giải | Tốc độ |
|---|---|---|---|
| 14 Frank–Wolfe | $`\min f`$, ràng buộc tuyến tính, miền bị chặn | 1 bài QHTT | $`O(1/k)`$ |
| 15 Rosen / Murtagh–Sargent | $`\min f`$, $`g_i \le 0`$, $`h_j = 0`$ | chiếu gradient lên tập tích cực | tuyến tính / siêu tuyến tính |
| 16 Newton có ràng buộc | $`\min f`$, $`Ax = b`$ | 1 hệ KKT tuyến tính | bậc hai |
| 17 Hàm chắn | $`\min f_0`$, $`f_i \le 0`$, $`Ax = b`$ | định tâm bằng Newton câu 16 | khe $`m/t`$ giảm theo cấp số nhân |

---

## Câu 14. Frank–Wolfe

### Lý thuyết cốt lõi

Bài toán $`\min f(x)`$ với $`x \in D = \{Ax \le b,\ A_{eq}x = b_{eq},\ x \ge lb\}`$, trong đó $`f`$ lồi khả vi và $`D`$ là đa diện **bị chặn**.

```math
y_k = \arg\min_{y \in D} \nabla f(x_k)^T y \ \ \text{(QHTT)},\qquad d_k = y_k - x_k,\qquad x_{k+1} = x_k + \gamma_k d_k,\ \ \gamma_k \in [0, 1].
```

- **Khe Frank–Wolfe** $`G_k = \nabla f(x_k)^T (x_k - y_k) \ge f(x_k) - f^*`$. Nó là chứng chỉ sai số mà không cần biết $`f^*`$, và $`G_k = 0`$ khi và chỉ khi $`x_k`$ tối ưu.
  *Chứng minh:* $`f^* \ge f(x_k) + \nabla f(x_k)^T(x^* - x_k) \ge f(x_k) + \nabla f(x_k)^T(y_k - x_k) = f(x_k) - G_k`$.
- **Chọn bước:** tìm kiếm chính xác $`\gamma_k = \arg\min_{[0,1]} f(x_k + \gamma d_k)`$, hoặc bước mở $`\gamma_k = 2/(k+2)`$ (cho $`k = 0, 1, 2, 3`$ được $`1, 0.667, 0.5, 0.4`$).
- **Hội tụ:** $`f(x_k) - f^* \le \dfrac{2C}{k + 2}`$ với $`C = L\,\mathrm{diam}(D)^2`$, tức tốc độ $`O(1/k)`$. Khi $`x^*`$ nằm giữa một cạnh, dãy **zic-zac** giữa hai đầu mút của cạnh đó.
- **Bước lùi (away-step):** viết $`x_k = \sum_{v \in S_k} w_v v`$, xét thêm hướng $`d^A = x_k - v_k`$ với $`v_k = \arg\max_{v \in S_k} \nabla f(x_k)^T v`$ (ra xa đỉnh xấu nhất), $`\gamma \le w_{v_k}/(1 - w_{v_k})`$. Biến thể này loại được đỉnh xấu và cho hội tụ tuyến tính khi $`f`$ lồi mạnh.
- **Tính thưa:** $`x_k`$ là tổ hợp lồi của tối đa $`k + 1`$ đỉnh.
- **Với hàm bậc hai** $`f = \tfrac12 x^THx + c^Tx`$: tìm kiếm chính xác có công thức $`\gamma = \min\big(1,\ -\nabla f^T d / (d^T H d)\big)`$.

### Thuật toán

1. Chọn $`x_0 \in D`$ và $`\varepsilon`$.
2. Giải QHTT để được $`y_k`$. Nếu $`D`$ có ít đỉnh, chỉ cần tính $`\nabla f(x_k)^T v`$ tại mọi đỉnh $`v`$ rồi lấy đỉnh nhỏ nhất.
3. Tính $`G_k`$; nếu $`G_k \le \varepsilon`$ thì dừng.
4. $`d_k = y_k - x_k`$; chọn $`\gamma_k`$; $`x_{k+1} = x_k + \gamma_k d_k`$; quay lại bước 2.

### Dạng câu hỏi có thể gặp

| Dạng | Cách làm nhanh | Lệnh Matlab |
|---|---|---|
| "Dùng Frank–Wolfe giải … cho $`x^* =`$" | Đáp án là nghiệm QP/NLP: KKT tay hoặc solver | `solve_qp(H, c, A, b, [], [], lb)`; bảng lặp: `frank_wolfe(f, g, A, b, [], [], lb, x0, struct('display', true))` |
| Đỉnh $`y_0`$, khe $`G_0`$, bước $`\gamma_0`$, điểm $`x_1`$ | Tính $`\nabla f(x_0)^T v`$ tại các đỉnh; $`\varphi'(\gamma) = 0`$ | như trên với `'maxIter', 2` |
| Giá trị $`\gamma_k = 2/(k+2)`$; cận sai số từ khe | Thế số; $`f(x_k) - f^* \le G_k`$ | — |
| Lý thuyết: tốc độ, vì sao zic-zac, bước lùi | Thuộc mục lý thuyết | — |

### Câu hỏi mẫu

**14.1.** Bài toán $`\min (x_1 - 2.5)^2 + 4(x_2 - 1.25)^2`$ với $`x_1 + 2x_2 \le 4`$, $`x \ge 0`$. Frank–Wolfe xuất phát từ $`x_0 = (1, 0)`$. Đỉnh $`y_0`$ của bài toán con tuyến tính là:
A. $`(0, 0)`$ B. $`(4, 0)`$ C. $`(0, 2)`$ D. $`(2, 1)`$

<details><summary>Đáp án</summary>

**C.** $`\nabla f(x_0) = (-3, -10)`$. Các đỉnh của $`D`$ là $`(0,0), (4,0), (0,2)`$, cho $`\nabla f^T v = 0, -12, -20`$, nên $`y_0 = (0, 2)`$. Khe $`G_0 = (-3,-10)^T(1,-2) = 17`$. Với $`d_0 = (-1, 2)`$: $`\varphi(\gamma) = (1.5 + \gamma)^2 + 4(2\gamma - 1.25)^2`$, $`\varphi' = 34\gamma - 17 = 0`$, nên $`\gamma_0 = 0.5`$ và $`x_1 = (0.5, 1)`$. Chương trình in đúng $`y_0 = (0, 2)`$, $`G_0 = 17`$, $`\gamma_0 = 0.5`$, $`x_1 = (0.5, 1)`$.

</details>

**14.2.** Dùng Frank–Wolfe giải bài toán ở câu 14.1. Nghiệm $`x^* = \ ?`$, $`f^* = \ ?`$

<details><summary>Đáp án</summary>

$`x^* = (2, 1)`$, $`f^* = 0.5`$. Tâm $`(2.5, 1.25)`$ vi phạm ràng buộc, nên ràng buộc tích cực. KKT cho $`2(x_1 - 2.5) + \lambda = 0`$, $`8(x_2 - 1.25) + 2\lambda = 0`$, $`x_1 + 2x_2 = 4`$, suy ra $`\lambda = 1`$. `solve_qp` cho $`(2, 1)`$, $`f^* = 0.5`$.
Lưu ý: FW thường (tìm kiếm chính xác) zic-zac giữa $`(0,2)`$ và $`(4,0)`$, cần khoảng 16000 bước mới đạt $`G \le 10^{-3}`$. FW có bước lùi đạt nghiệm sau 5 bước. Khi làm bài, ghi bảng vài bước đầu rồi nêu nghiệm đúng.

</details>

**14.3.** Bài toán $`\min (x_1 - 1)^2 + (x_2 - 2)^2`$ với $`x_1 + x_2 \le 2`$, $`x \ge 0`$, $`x_0 = (0, 0)`$, tìm kiếm chính xác. Điểm $`x_2`$ (sau 2 bước) là:
A. $`(0, 2)`$ B. $`(0.5, 1.5)`$ C. $`(1, 1)`$ D. $`(2, 0)`$

<details><summary>Đáp án</summary>

**B.**
- **Bước 0:** $`\nabla f(0,0) = (-2, -4)`$; tại các đỉnh $`(0,0), (2,0), (0,2)`$ được $`0, -4, -8`$, nên $`y_0 = (0,2)`$, $`G_0 = 8`$. $`\varphi(\gamma) = 1 + (2\gamma - 2)^2`$ đạt cực tiểu tại $`\gamma = 1`$, nên $`x_1 = (0, 2)`$.
- **Bước 1:** $`\nabla f(0,2) = (-2, 0)`$, nên $`y_1 = (2, 0)`$, $`G_1 = 4`$. $`d = (2, -2)`$, $`\varphi(\gamma) = (2\gamma - 1)^2 + 4\gamma^2`$, $`\varphi' = 16\gamma - 4 = 0`$, nên $`\gamma = 0.25`$ và $`x_2 = (0.5, 1.5)`$.
- **Bước 2:** $`G_2 = 0`$, nên $`x_2`$ tối ưu (chính là hình chiếu của $`(1, 2)`$ lên $`x_1 + x_2 = 2`$, với $`\lambda = 1`$, $`f^* = 0.5`$).

Chương trình in đúng: khe $`(8, 4, 0)`$, các bước $`\gamma = (1, 0.25)`$.

</details>

### Nếu ra tự luận

*"Trình bày giải thuật Frank–Wolfe và thực hiện 2 bước cho ví dụ."* Dàn ý:
1. Phát biểu bài toán (ràng buộc tuyến tính, $`D`$ bị chặn).
2. Ý tưởng: tuyến tính hóa $`f`$ tại $`x_k`$; cực tiểu của hàm tuyến tính trên đa diện đạt tại một đỉnh.
3. Các bước của thuật toán (mục trên) và hai cách chọn bước.
4. Chứng minh khe là chứng chỉ sai số (3 dòng ở trên); nêu tốc độ $`O(1/k)`$ và hiện tượng zic-zac, cách khắc phục bằng bước lùi.
5. Ví dụ: lấy câu 14.3 ở trên (2 bước ra nghiệm đúng, dễ trình bày) hoặc câu 14.1 (minh họa zic-zac).

Báo cáo đầy đủ: [câu 14](../../../solutions/q14-frank-wolfe/).

### Bẫy

- Phải xét $`\nabla f^T v`$ tại **mọi** đỉnh của $`D`$, kể cả $`(0, 0)`$.
- Với tìm kiếm chính xác, $`\gamma`$ phải nằm trong $`[0, 1]`$: nếu cực tiểu của $`\varphi`$ lớn hơn 1 thì lấy $`\gamma = 1`$ (bước 0 của câu 14.3).
- Dừng theo khe $`G_k`$, không theo $`\Vert x_{k+1} - x_k \Vert`$: bước có thể rất nhỏ trong khi còn cách xa nghiệm.
- Nếu $`D`$ không bị chặn thì bài QHTT con có thể không có nghiệm hữu hạn, và Frank–Wolfe không áp dụng được.

---

## Câu 15. Chiếu gradient Rosen và Murtagh–Sargent

### Lý thuyết cốt lõi

Tại điểm chấp nhận được $`x`$, **tập làm việc** $`W`$ gồm các đẳng thức và các bất đẳng thức tích cực. $`N`$ là ma trận mà mỗi dòng là một pháp vector $`\nabla g_i^T`$ với $`i \in W`$ (cần hạng đủ, tức LICQ).

```math
P = I - N^T (N N^T)^{-1} N,\qquad d = -P\nabla f(x),\qquad u = -(N N^T)^{-1} N \nabla f(x)\ \ \Rightarrow\ \ \nabla f + N^T u = -d .
```

- $`P = P^T = P^2`$ và $`NP = 0`$: $`P`$ là phép chiếu trực giao lên $`\{d : Nd = 0\}`$. Nếu $`d \ne 0`$ thì $`\nabla f^T d = -\Vert P\nabla f \Vert^2 \lt 0`$, tức $`d`$ là hướng giảm.
- Khi **$`d = 0`$**: nếu $`u_i \ge 0`$ với mọi bất đẳng thức trong $`W`$ thì $`x`$ là điểm KKT, tức tối ưu (bài toán lồi). Nếu có $`u_j \lt 0`$ thì **bỏ ràng buộc có $`u_j`$ âm nhất** khỏi $`W`$, tính lại $`P`$; khi đó $`\hat d \ne 0`$ và $`a_j^T \hat d \lt 0`$ (rời khỏi biên $`j`$).
- **Bước lớn nhất** giữ tính chấp nhận được:

```math
\alpha_{\max} = \min\Big\{ \frac{b_i - a_i^T x}{a_i^T d} : i \notin W,\ a_i^T d \gt 0 \Big\},\qquad \alpha = \arg\min_{0 \le \alpha \le \alpha_{\max}} f(x + \alpha d).
```

  Nếu $`\alpha = \alpha_{\max}`$ thì thêm ràng buộc chặn vào $`W`$.
- **Ràng buộc phi tuyến:** $`N`$ là ma trận Jacobi; sau mỗi bước phải phục hồi về mặt ràng buộc bằng Newton, $`z \leftarrow z - N^T(NN^T)^{-1} g_W(z)`$.
- **Murtagh–Sargent** (metric biến thiên): thay $`P`$ bằng

```math
P_H = H - H N^T (N H N^T)^{-1} N H,\qquad d = -P_H \nabla f,\qquad u = -(N H N^T)^{-1} N H \nabla f,
```

  với $`H \approx`$ nghịch đảo Hessian của hàm Lagrange, cập nhật bằng BFGS. Với $`H = I`$ thì trở về Rosen; với $`W = \varnothing`$ thì là tựa Newton của câu 12. Tốc độ: Rosen hội tụ tuyến tính (là gradient trên không gian con), M–S hội tụ siêu tuyến tính.

### Thuật toán

1. Chọn $`x_0`$ **chấp nhận được**; $`W`$ là tập tích cực tại $`x_0`$.
2. Tính $`N`$, $`P`$ (hoặc $`P_H`$), $`d`$, $`u`$.
3. Nếu $`d = 0`$: dừng nếu $`u \ge 0`$ trên các bất đẳng thức; nếu không, bỏ ràng buộc có $`u_j`$ âm nhất và quay lại bước 2.
4. Tính $`\alpha_{\max}`$ và $`\alpha`$; $`x \leftarrow x + \alpha d`$; thêm ràng buộc mới trở thành tích cực vào $`W`$; (M–S: cập nhật $`H`$); quay lại bước 2.

### Dạng câu hỏi có thể gặp

| Dạng | Cách làm nhanh | Lệnh Matlab |
|---|---|---|
| "Dùng Rosen / Murtagh–Sargent giải … cho $`x^* =`$" | Nghiệm QP/NLP (KKT tay hoặc solver) | `solve_qp`, `solve_nlp`; bảng lặp: `gradient_projection(prob, x0, struct('lineSearch', 'exact', 'display', true))`, thêm `'metric', 'bfgs'` cho M–S |
| Tính $`P`$ cho $`N`$ cho trước | $`N`$ một dòng $`a^T`$: $`P = I - aa^T/(a^Ta)`$ | `P = eye(n) - N'*((N*N')\N)` |
| Hướng $`d_0`$, nhân tử $`u`$, ràng buộc bị bỏ | Công thức trên; bỏ $`u_j`$ âm nhất | `u = -(N*N')\(N*g)` |
| $`\alpha_{\max}`$, điểm $`x_1`$ | Tỉ số $`(b_i - a_i^Tx)/(a_i^Td)`$ với $`a_i^Td \gt 0`$ | `fminbnd(@(a) f(x + a*d), 0, amax)` |
| Lý thuyết: tính chất của $`P`$, điều kiện dừng, khác nhau giữa Rosen và M–S | Mục lý thuyết | — |

### Câu hỏi mẫu

Bài toán chung cho 15.1–15.3 (Bazaraa):

```math
\min f = 2x_1^2 + 2x_2^2 - 2x_1x_2 - 4x_1 - 6x_2\quad \text{với}\quad \underbrace{x_1 + x_2 \le 2}_{(1)},\ \underbrace{x_1 + 5x_2 \le 5}_{(2)},\ \underbrace{-x_1 \le 0}_{(3)},\ \underbrace{-x_2 \le 0}_{(4)},\qquad x_0 = (0, 0).
```

**15.1.** Tại $`x_0 = (0,0)`$ với $`W = \{3, 4\}`$, phương pháp Rosen bỏ ràng buộc nào và điểm $`x_1`$ là:
A. Bỏ (3), $`x_1 = (1, 0)`$ B. Bỏ (4), $`x_1 = (0, 1)`$ C. Bỏ (4), $`x_1 = (0, 1.5)`$ D. Không bỏ, $`x_0`$ tối ưu

<details><summary>Đáp án</summary>

**B.** $`N = -I`$ nên $`P = 0`$, $`d = 0`$. $`\nabla f(x_0) = (-4, -6)`$ cho $`u = (-4, -6)`$, cả hai âm, nên bỏ (4) (âm nhất). Với $`W = \{3\}`$: $`P = \mathrm{diag}(0, 1)`$, $`d = (0, 6)`$. Bước chặn: (1) cho $`\alpha \le 2/6`$, (2) cho $`\alpha \le 5/30`$, nên $`\alpha_{\max} = 1/6`$. Cực tiểu của $`f(0, 6\alpha)`$ nằm tại $`\alpha = 0.25 \gt 1/6`$, nên $`\alpha = 1/6`$ và $`x_1 = (0, 1)`$; ràng buộc (2) trở thành tích cực.

</details>

**15.2.** Tại $`x_1 = (0, 1)`$, sau khi bỏ (3) còn $`W = \{2\}`$ với $`N = (1, 5)`$. Ma trận chiếu $`P`$ và hướng $`d`$ là:

<details><summary>Đáp án</summary>

Với $`W = \{2, 3\}`$ trước đó, $`u = (0.4, -5.6)`$, nên bỏ (3). Với $`N = (1, 5)`$ thì $`NN^T = 26`$ và $`N^TN`$ có các dòng $`(1, 5)`$, $`(5, 25)`$, nên $`P = I - N^TN/26`$ có các dòng $`(25, -5)/26`$ và $`(-5, 1)/26`$. Với $`\nabla f(x_1) = (-6, -2)`$: $`d = -P\nabla f(x_1) = (140, -28)/26 \approx (5.385, -1.077)`$, cùng hướng với $`(5, -1)`$.

Chương trình in đúng $`26P = [25\ {-5};\ {-5}\ 1]`$ và $`26d = (140, -28)`$.

</details>

**15.3.** Dùng phương pháp Rosen giải bài toán trên. Nghiệm $`x^* = \ ?`$, nhân tử của ràng buộc (2) là $`?`$

<details><summary>Đáp án</summary>

$`x^* = (35/31, 24/31) \approx (1.1290, 0.7742)`$, $`f^* = -222/31 \approx -7.16129`$, $`\lambda_2 = 32/31 \approx 1.0323`$ (các nhân tử khác bằng 0). Giải tay KKT với chỉ (2) tích cực: $`4x_1 - 2x_2 - 4 + \lambda = 0`$, $`4x_2 - 2x_1 - 6 + 5\lambda = 0`$, $`x_1 + 5x_2 = 5`$. Rosen với tìm kiếm chính xác cần 2 bước (bước 1: $`\alpha = 0.2097 \lt \alpha_{\max} = 0.2321`$), Murtagh–Sargent cũng 2 bước, và hàm chắn (câu 17) cho cùng kết quả.

</details>

### Nếu ra tự luận

*"Trình bày giải thuật Rosen (Murtagh–Sargent) và áp dụng cho ví dụ."* Dàn ý:
1. Bài toán và tập làm việc.
2. Ma trận chiếu cùng 3 tính chất ($`P^2 = P`$, $`NP = 0`$, $`d`$ là hướng giảm, kèm chứng minh 1 dòng).
3. Định lý dừng và quy tắc bỏ ràng buộc.
4. Công thức $`\alpha_{\max}`$.
5. Các bước thuật toán.
6. Phần M–S: công thức $`P_H`$, cập nhật BFGS cho $`H`$, lý do hội tụ nhanh hơn.
7. Ví dụ Bazaraa: trình bày bước 0 (câu 15.1) và bước 1 (câu 15.2), kết luận ở câu 15.3.

Báo cáo đầy đủ: [câu 15](../../../solutions/q15-gradient-projection/).

### Bẫy

- $`x_0`$ phải **chấp nhận được**; nếu chưa có thì phải tìm trước (pha 1).
- $`d = 0`$ chưa chắc là tối ưu: phải xét **dấu của $`u`$**.
- Chỉ bỏ **một** ràng buộc mỗi lần (ràng buộc có $`u`$ âm nhất).
- $`\alpha_{\max}`$ chỉ xét các ràng buộc **ngoài** $`W`$ có $`a_i^T d \gt 0`$.
- Với ràng buộc cong, đi theo $`d`$ sẽ ra khỏi miền, nên cần bước phục hồi.

---

## Câu 16. Newton cho quy hoạch lồi với ràng buộc $`Ax = b`$

### Lý thuyết cốt lõi

Bài toán $`\min f(x)`$ với $`Ax = b`$, trong đó $`f`$ lồi, khả vi hai lần, $`A \in \mathbb{R}^{p \times n}`$ có hạng $`p`$. Điều kiện tối ưu (không cần điều kiện chính quy vì ràng buộc affine):

```math
\nabla f(x^*) + A^T\nu^* = 0,\qquad Ax^* = b .
```

**Hàm bậc hai** $`f = \tfrac12 x^TPx + q^Tx`$: nghiệm thu được ngay từ hệ KKT tuyến tính

```math
\begin{pmatrix} P & A^T \\ A & 0 \end{pmatrix}\begin{pmatrix} x^* \\ \nu^* \end{pmatrix} = \begin{pmatrix} -q \\ b \end{pmatrix}.
```

**Xuất phát chấp nhận được** ($`Ax_0 = b`$). Bước Newton $`\Delta x`$ là nghiệm của

```math
\begin{pmatrix} \nabla^2 f(x) & A^T \\ A & 0 \end{pmatrix}\begin{pmatrix} \Delta x \\ w \end{pmatrix} = \begin{pmatrix} -\nabla f(x) \\ 0 \end{pmatrix},\qquad \lambda(x)^2 = \Delta x^T \nabla^2 f(x)\,\Delta x = -\nabla f(x)^T\Delta x .
```

- $`A\Delta x = 0`$ nên mọi điểm lặp đều chấp nhận được.
- $`\lambda^2/2 = f(x) - \min \hat f`$ (với $`\hat f`$ là xấp xỉ bậc hai); dừng khi $`\lambda^2/2 \le \varepsilon`$.
- Tìm kiếm quay lui: $`f(x + t\Delta x) \le f(x) - \alpha t\lambda^2`$.
- Phương pháp tương đương với Newton không ràng buộc sau khi khử biến, nên hội tụ bậc hai gần nghiệm; với hàm bậc hai thì đạt nghiệm sau **1 bước**.

**Xuất phát bất kỳ.** Phần dư

```math
r(x, \nu) = \begin{pmatrix} \nabla f(x) + A^T\nu \\ Ax - b \end{pmatrix},\qquad
\begin{pmatrix} \nabla^2 f & A^T \\ A & 0 \end{pmatrix}\begin{pmatrix} \Delta x \\ \Delta\nu \end{pmatrix} = -\,r(x, \nu).
```

- Quay lui theo $`\Vert r \Vert`$: $`\Vert r(y + t\Delta y) \Vert \le (1 - \alpha t)\Vert r(y) \Vert`$.
- $`Ax^+ - b = (1 - t)(Ax - b)`$, nên chỉ cần một bước $`t = 1`$ là mọi điểm sau đều chấp nhận được.
- $`f`$ có thể tăng trong quá trình lặp, vì tìm kiếm nhằm giảm $`\Vert r \Vert`$ chứ không phải $`f`$.

### Thuật toán

1. **Chấp nhận được:** giải hệ KKT để được $`\Delta x`$ và $`\lambda^2`$ → nếu $`\lambda^2/2 \le \varepsilon`$ thì dừng → quay lui trên $`f`$ để được $`t`$ → $`x \leftarrow x + t\Delta x`$.
2. **Bất kỳ:** giải hệ với vế phải $`-r`$ → quay lui trên $`\Vert r \Vert`$ → $`(x, \nu) \leftarrow (x, \nu) + t(\Delta x, \Delta\nu)`$ → dừng khi $`\Vert Ax - b \Vert`$ và $`\Vert r \Vert`$ đều $`\le \varepsilon`$.

### Dạng câu hỏi có thể gặp

| Dạng | Cách làm nhanh | Lệnh Matlab |
|---|---|---|
| "Dùng Newton giải $`\min f`$, $`Ax = b`$, cho $`x^* =`$" | Bậc hai: giải hệ KKT tuyến tính. Không bậc hai: khử biến bằng $`Ax = b`$ rồi giải $`\nabla = 0`$ | `newton_eq(f, g, h, A, b, x0, struct('display', true))` |
| Bước $`\Delta x`$, $`w`$, $`\lambda^2`$ tại điểm cho trước | Lập ma trận KKT cỡ $`n + p`$, giải | `K = [H A'; A zeros(p)]; s = K \ [-g; 0]` |
| Phần dư sau một bước $`t`$ | $`\Vert Ax^+ - b \Vert = (1 - t)\Vert Ax - b \Vert`$ | — |
| Lý thuyết: vì sao điểm lặp chấp nhận được, ý nghĩa của $`\lambda^2/2`$ | Mục lý thuyết | — |

### Câu hỏi mẫu

**16.1.** Bài toán $`\min \tfrac12\Vert x \Vert^2`$ với $`x_1 + x_2 + x_3 = 3`$, $`x_1 - x_3 = 1`$, xuất phát từ $`x_0 = (2, 0, 1)`$. Giá trị $`\lambda(x_0)^2/2`$ và nghiệm $`x^*`$ là:
A. $`0.75`$; $`(1.5, 1, 0.5)`$ B. $`1.5`$; $`(1.5, 1, 0.5)`$ C. $`0.75`$; $`(1, 1, 1)`$ D. $`2.5`$; $`(2, 0, 1)`$

<details><summary>Đáp án</summary>

**A.** Điều kiện tối ưu $`x + A^T\nu = 0`$ cho $`x = -A^T\nu`$, và $`AA^T = \mathrm{diag}(3, 2)`$, nên $`\nu = -(AA^T)^{-1}b = (-1, -0.5)`$, $`x^* = (1.5, 1, 0.5)`$, $`f^* = 1.75`$. Vì $`f`$ bậc hai nên $`\lambda^2/2 = f(x_0) - f^* = 2.5 - 1.75 = 0.75`$, và một bước $`t = 1`$ cho nghiệm. Xuất phát bất kỳ từ $`x_0 = 0`$ ($`\Vert Ax_0 - b \Vert = 3.1623`$) cũng chỉ cần 1 bước.

</details>

**16.2.** Bài toán $`\min x_1^4 + x_2^2`$ với $`x_1 + x_2 = 2`$, tại $`x_0 = (1, 1)`$. Tính bước Newton $`\Delta x`$, $`\lambda^2/2`$ và nghiệm.

<details><summary>Đáp án</summary>

$`\nabla f = (4, 2)`$, $`\nabla^2 f = \mathrm{diag}(12, 2)`$. Hệ KKT gồm ba phương trình $`12\Delta x_1 + w = -4`$, $`2\Delta x_2 + w = -2`$, $`\Delta x_1 + \Delta x_2 = 0`$.
Đặt $`\Delta x = (s, -s)`$: $`12s + w = -4`$ và $`-2s + w = -2`$ cho $`s = -1/7`$. Vậy $`\Delta x = (-1/7, 1/7)`$, $`w = -16/7`$, $`\lambda^2 = (12 + 2)/49 = 2/7`$, $`\lambda^2/2 = 1/7 \approx 0.142857`$.
Nghiệm: khử biến $`x_2 = 2 - x_1`$, điều kiện $`4x_1^3 = 2(2 - x_1)`$ tức $`2x_1^3 + x_1 - 2 = 0`$, cho $`x_1 = 0.835122`$. Do đó $`x^* \approx (0.8351, 1.1649)`$, $`f^* \approx 1.843348`$, $`\nu^* = -2x_2 \approx -2.32976`$. `newton_eq` cần 3 bước, với $`\lambda^2/2`$: $`0.143 \to 2.5\cdot10^{-3} \to 1.1\cdot10^{-6} \to 2\cdot10^{-13}`$ (hội tụ bậc hai).

</details>

**16.3.** Newton xuất phát không chấp nhận được, $`\Vert Ax_0 - b \Vert = 3.1623`$. Sau một bước với $`t = 0.5`$, $`\Vert Ax_1 - b \Vert`$ bằng:
A. 0 B. 1.5811 C. 3.1623 D. Không xác định được

<details><summary>Đáp án</summary>

**B.** $`Ax^+ - b = (1 - t)(Ax - b)`$, nên $`0.5 \times 3.1623 = 1.5811`$. Đáp án A chỉ đúng khi $`t = 1`$.

</details>

### Nếu ra tự luận

*"Trình bày thuật toán Newton cho quy hoạch lồi ràng buộc tuyến tính, với điểm xuất phát chấp nhận được và bất kỳ."* Dàn ý:
1. Điều kiện tối ưu.
2. Bước Newton là nghiệm của bài toán con bậc hai, dẫn đến hệ KKT.
3. Mệnh đề: $`A\Delta x = 0`$, $`\nabla f^T\Delta x = -\lambda^2`$, ý nghĩa của $`\lambda^2/2`$.
4. Giải thuật 1 (chấp nhận được).
5. Phần dư $`r`$ và giải thuật 2 (bất kỳ), tính chất $`(1 - t)`$.
6. Tốc độ bậc hai (qua khử biến).
7. Ví dụ: câu 16.2, tính trọn bước đầu như trên.

Báo cáo đầy đủ: [câu 16](../../../solutions/q16-newton-equality-constrained/).

### Bẫy

- Vế phải dòng dưới bằng 0 khi xuất phát chấp nhận được, và bằng $`-(Ax - b)`$ khi xuất phát bất kỳ.
- Trong hệ KKT phải dùng $`\nabla^2 f`$, không phải $`I`$.
- $`\lambda^2/2`$ là lượng giảm dự báo, khác $`\Vert \nabla f \Vert`$ (vì $`\nabla f`$ không bằng 0 tại nghiệm có ràng buộc).
- Ngoài miền xác định (ví dụ $`\log x`$ với $`x \le 0`$) phải giảm $`t`$.

---

## Câu 17. Hàm chắn Newton

### Lý thuyết cốt lõi

Bài toán $`\min f_0(x)`$ với $`f_i(x) \le 0`$ ($`i = 1..m`$) và $`Ax = b`$, lồi, thỏa **Slater**. Hàm chắn logarit và bài toán định tâm:

```math
\phi(x) = -\sum_{i=1}^m \log(-f_i(x)),\qquad x^*(t) = \arg\min\ t f_0(x) + \phi(x)\ \ \text{với}\ \ Ax = b .
```

- $`\phi`$ lồi và tiến ra $`+\infty`$ ở biên, nên các điểm lặp luôn ở trong miền.
- Đạo hàm: $`\nabla\phi = \sum_i \dfrac{\nabla f_i}{-f_i}`$, $`\nabla^2\phi = \sum_i \dfrac{\nabla f_i \nabla f_i^T}{f_i^2} + \sum_i \dfrac{\nabla^2 f_i}{-f_i}`$.
- **Đường trung tâm** $`\{x^*(t) : t \gt 0\}`$. Điểm đối ngẫu và khe:

```math
\lambda_i^*(t) = -\frac{1}{t\,f_i(x^*(t))} \gt 0,\qquad f_0(x^*(t)) - p^* \le f_0(x^*(t)) - q(\lambda^*(t), \nu^*(t)) = \frac{m}{t}.
```

  Tức là KKT đúng hết, trừ điều kiện bù được thay bằng $`-\lambda_i f_i = 1/t`$.
- **Số bước ngoài:** từ $`t^{(0)}`$, nhân $`t \leftarrow \mu t`$ cho đến khi $`m/t \le \varepsilon`$, cần $`\Big\lceil \dfrac{\log\big(m/(t^{(0)}\varepsilon)\big)}{\log\mu} \Big\rceil`$ lần nhân, cộng thêm bước đầu.
- **Pha 1:** giải $`\min s`$ với $`f_i(x) \le s`$, $`Ax = b`$ (và $`s \ge -1`$). Kết luận: $`s^* \lt 0`$ thì có điểm chấp nhận chặt; $`s^* \gt 0`$ thì bài toán vô nghiệm; $`s^* = 0`$ thì không có điểm chấp nhận chặt (Slater sai), không dùng hàm chắn được.
- **Chọn $`\mu`$:** $`\mu`$ lớn thì ít bước ngoài nhưng mỗi lần định tâm cần nhiều bước Newton hơn. Trong thực tế $`\mu = 10`$ đến $`100`$ cho tổng số bước Newton gần như nhau.

### Thuật toán

1. Có $`x`$ chấp nhận chặt (nếu chưa có thì chạy pha 1); chọn $`t = t^{(0)}`$, $`\mu \gt 1`$, $`\varepsilon`$.
2. **Định tâm:** dùng Newton có ràng buộc đẳng thức (câu 16) giải $`\min t f_0 + \phi`$, $`Ax = b`$, xuất phát từ $`x`$ hiện tại.
3. $`x \leftarrow x^*(t)`$. Nếu $`m/t \le \varepsilon`$ thì dừng.
4. $`t \leftarrow \mu t`$, quay lại bước 2.

### Dạng câu hỏi có thể gặp

| Dạng | Cách làm nhanh | Lệnh Matlab |
|---|---|---|
| "Dùng hàm chắn giải … cho $`x^* =`$" | Nghiệm của bài toán lồi (KKT, solver) | `solve_qp` / `solve_nlp`; bảng lặp: `barrier_method(prob, x0, struct('display', true))` |
| Khe sau $`k`$ bước, số bước ngoài | $`m/t`$; công thức $`\lceil\log(m/(t^{(0)}\varepsilon))/\log\mu\rceil`$ | `ceil(log(m/(t0*eps))/log(mu))` |
| $`x^*(t)`$ cho bài nhỏ (1–2 biến) | Giải $`t\nabla f_0 + \nabla\phi = 0`$ | `fminbnd(@(x) t*f0(x) - log(-f1(x)), lo, hi)` |
| Ước lượng nhân tử $`\lambda_i = -1/(t f_i)`$ | Thế số | `out.lambda` |
| Lý thuyết: đường trung tâm, pha 1, ảnh hưởng của $`\mu`$ | Mục lý thuyết | — |

### Câu hỏi mẫu

**17.1.** Giải bài toán Bazaraa (câu 15) bằng hàm chắn với $`m = 4`$, $`t^{(0)} = 1`$, $`\mu = 10`$, $`\varepsilon = 10^{-8}`$. Số bước ngoài và khe sau bước ngoài thứ 3 là:
A. 9 và 0.04 B. 10 và 0.04 C. 10 và 0.4 D. 8 và 0.004

<details><summary>Đáp án</summary>

**B.** $`\lceil \log_{10}(4/10^{-8}) \rceil = \lceil 8.6 \rceil = 9`$ lần nhân, cộng bước đầu thành 10. Các giá trị $`t`$ là $`1, 10, 100`$, nên khe sau bước 3 là $`m/t = 4/100 = 0.04`$. Chương trình: 10 bước ngoài, tổng 61 bước Newton, khe $`(4, 0.4, 0.04, 0.004, \dots)`$. Kết quả $`x = (1.1290323, 0.7741935)`$ và $`\lambda \to (0, 1.0323, 0, 0)`$, đúng nghiệm của câu 15.

</details>

**17.2.** Bài toán $`\min x`$ với $`1 - x \le 0`$. Điểm trung tâm $`x^*(t)`$ và nhân tử ước lượng $`\lambda(t)`$ là:

<details><summary>Đáp án</summary>

Bài toán định tâm $`\min\ tx - \log(x - 1)`$ có điều kiện $`t - 1/(x - 1) = 0`$, nên $`x^*(t) = 1 + 1/t`$. Khe $`f_0(x^*(t)) - p^* = 1/t = m/t`$ (đúng bằng cận, với $`m = 1`$). Nhân tử $`\lambda = -1/(t(1 - x^*(t))) = 1`$ với mọi $`t`$, đúng bằng nhân tử KKT. Kiểm tra bằng `fminbnd`: $`t = 1, 10, 100`$ cho $`x^*(t) = 2, 1.1, 1.01`$.

</details>

**17.3.** Bài toán $`\min (x_1 - 3)^2 + (x_2 - 2)^2`$ với $`x_1^2 + x_2^2 \le 5`$, $`x_1 + 2x_2 \le 4`$, xuất phát $`x_0 = (3, 3)`$ (không chấp nhận được). Nghiệm của phương pháp hàm chắn là $`x^* = \ ?`$ Pha 1 kết luận gì?

<details><summary>Đáp án</summary>

Pha 1 cho $`s = -0.8948 \lt 0`$, nên tìm được điểm chấp nhận chặt và chuyển sang pha 2. Kết quả $`x^* = (2, 1)`$, $`f^* = 2`$, $`\lambda = (1/3, 2/3)`$, khớp với KKT giải tay (cả hai ràng buộc tích cực tại $`(2, 1)`$: $`2(x - c) + \lambda_1(4, 2) + \lambda_2(1, 2) = 0`$).

</details>

### Nếu ra tự luận

*"Trình bày thuật toán hàm chắn Newton giải quy hoạch lồi tổng quát."* Dàn ý:
1. Bài toán và giả thiết Slater.
2. Thay hàm chỉ $`I_-`$ bằng $`-\tfrac1t\log(-u)`$, dẫn đến bài toán định tâm.
3. Tính lồi của $`\phi`$; $`\nabla\phi`$, $`\nabla^2\phi`$.
4. Đường trung tâm, định lý khe $`m/t`$ (chứng minh: $`x^*(t)`$ cực tiểu hàm Lagrange với $`\lambda^*(t)`$, nên $`q = f_0 - m/t`$).
5. Thuật toán: vòng ngoài, định tâm bằng Newton câu 16, số bước ngoài.
6. Pha 1.
7. Ví dụ: câu 17.2 (giải tay trọn vẹn) và câu 17.1 (bảng $`t`$, khe).

Báo cáo đầy đủ: [câu 17](../../../solutions/q17-barrier-method/).

### Bẫy

- Hàm chắn cần điểm **chấp nhận chặt** ($`f_i \lt 0`$), không phải chỉ chấp nhận được.
- Khe là $`m/t`$ với $`m`$ là **số ràng buộc bất đẳng thức**, không phải số biến.
- Số bước ngoài phải cộng thêm bước đầu tiên ($`t = t^{(0)}`$).
- $`m/t`$ là cận của $`f_0(x^*(t)) - p^*`$; sai số thực có thể nhỏ hơn (ví dụ 17.1 thực tế khoảng $`1/t`$).
- Nếu pha 1 cho $`s^* = 0`$ thì Slater sai, nên không kết luận "vô nghiệm" mà kết luận "không dùng được hàm chắn".
