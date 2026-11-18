# Tối ưu không ràng buộc nhiều biến (câu 8, 9, 10, 11, 12)

[← Tài liệu ôn thi](../README.md)

Ký hiệu: $`g_k = \nabla f(x_k)`$, $`H_k = \nabla^2 f(x_k)`$. Với hàm bậc hai $`f = \tfrac12 x^T H x + q^T x`$: $`\mu, L`$ là trị riêng nhỏ nhất và lớn nhất của $`H`$, $`\kappa = L/\mu`$. Số liệu trong các câu hỏi mẫu lấy từ output của [`check_unconstrained.m`](../check_unconstrained.m).

**Nguyên tắc chung cho câu điền đáp án:** nếu $`f`$ lồi chặt thì mọi phương pháp hội tụ đều cho cùng một $`x^*`$. Vì vậy hãy tìm đáp án bằng cách nhanh nhất (`-H\q` hoặc `fminunc`/`fminsearch`), rồi chạy phương pháp đề yêu cầu để có bảng lặp ghi vào bài ([quy trình](../README.md#3-quy-trình-làm-một-câu-trắc-nghiệm--điền-đáp-án)).

| Phương pháp | Công thức lặp | Tốc độ | Cần |
|---|---|---|---|
| Gauss–Seidel (câu 8) | lần lượt $`x_i \leftarrow \mathrm{argmin}_t f(x_1, .., t, .., x_n)`$ | tuyến tính | $`f`$ |
| Hooke–Jeeves (câu 9) | dò $`\pm\Delta e_i`$ + bước mẫu $`x_n + (x_n - x_b)`$, giảm $`\Delta`$ | tuyến tính, chậm | chỉ $`f`$ |
| Gradient (câu 10) | $`x_{k+1} = x_k - \alpha_k g_k`$ | tuyến tính, tỉ số $`\frac{\kappa - 1}{\kappa + 1}`$ | $`\nabla f`$ |
| Quán tính (câu 10) | $`x_{k+1} = x_k - \alpha g_k + \beta (x_k - x_{k-1})`$ | tuyến tính, tỉ số $`\frac{\sqrt\kappa - 1}{\sqrt\kappa + 1}`$ | $`\nabla f`$ |
| Newton (câu 11) | $`H_k d_k = -g_k`$, $`x_{k+1} = x_k + d_k`$ | bậc hai | $`\nabla f`$, $`\nabla^2 f`$ |
| Tựa Newton (câu 12) | $`d_k = -B_k g_k`$, cập nhật $`B_k \approx H^{-1}`$ | siêu tuyến tính | $`\nabla f`$ |

Các ví dụ chung đã chạy trong bài tập lớn: Bazaraa $`F = (x_1 - 2)^4 + (x_1 - 2x_2)^2`$ từ $`x^0 = (0, 3)`$; Rosenbrock từ $`(-1.2, 1)`$; hàm bậc hai $`n = 10`$ với $`A = \mathrm{tridiag}(-1, 2, -1)`$ ($`\kappa = 48.37`$). Kết quả: Gauss–Seidel 207 vòng; Hooke–Jeeves 3409 lần tính $`f`$; gradient 472 bước; quán tính 78 bước; Newton 1 bước (hàm bậc hai); BFGS 5–10 bước.

---

## Câu 8. Gauss–Seidel (tọa độ tuần hoàn)

### Lý thuyết cốt lõi

- Mỗi vòng cực tiểu $`f`$ lần lượt theo $`x_1, \dots, x_n`$, các biến khác giữ nguyên **giá trị mới nhất**.
- Với hàm bậc hai, $`\partial f/\partial x_i = 0`$ cho công thức đóng (đây chính là phương pháp Gauss–Seidel giải hệ $`Hx = -q`$):

```math
x_i \leftarrow -\frac{q_i + \sum_{j \ne i} H_{ij} x_j}{H_{ii}} .
```

- **Hội tụ:** với $`H \succ 0`$, phương pháp luôn hội tụ; mỗi bước tọa độ làm $`f`$ không tăng. Phương pháp chậm khi các biến liên kết mạnh (thung lũng hẹp nằm chéo trục), vì đường đi zic-zac song song với các trục.
- Với hàm không bậc hai, mỗi bước là một bài toán một chiều, giải bằng tìm khoảng chứa cực tiểu kết hợp lát cắt vàng.
- Hàm không trơn có thể làm phương pháp kẹt ở điểm không phải cực tiểu.

### Thuật toán

1. Cho $`x^0`$ và $`\varepsilon`$.
2. Vòng $`k`$: với $`i = 1..n`$, $`x_i \leftarrow \mathrm{argmin}_t f(x_1, \dots, x_{i-1}, t, x_{i+1}, \dots, x_n)`$.
3. Dừng khi $`\Vert x_k - x_{k-1}\Vert \le \varepsilon`$ (hoặc khi độ giảm của $`f`$ đủ nhỏ).

### Dạng câu hỏi có thể gặp

| Dạng | Cách làm nhanh | Matlab |
|---|---|---|
| Điền $`x^*`$ (hàm bậc hai) | Lập $`H`$ ($`H_{ii} = 2 a_{ii}`$), giải $`x^* = -H^{-1}q`$ | `gs_quad(H,q,x0,1e-4)`, `-H\q` |
| Điền $`x^*`$ (hàm bất kỳ) | `fminsearch`, rồi chạy Gauss–Seidel | `gs_min(f,x0,1e-4)`, `gauss_seidel_min` |
| Điểm sau vòng 1 | Thay lần lượt vào công thức đóng | `gs_quad` (dòng $`k = 1`$) |

### Câu hỏi mẫu

**8.1.** Áp dụng Gauss–Seidel cho $`f = x_1^2 + x_1x_2 + x_2^2 - 3x_1`$ từ $`x^0 = (0, 0)`$. Điểm sau vòng 1 và cực tiểu là:
A. $`(1.5, -0.75)`$; $`(2, -1)`$ B. $`(1.5, 0)`$; $`(2, -1)`$ C. $`(1.5, -0.75)`$; $`(1, -0.5)`$ D. $`(3, -1.5)`$; $`(2, -1)`$

<details><summary>Đáp án</summary>

**A.** $`\partial f/\partial x_1 = 2x_1 + x_2 - 3 = 0`$ cho $`x_1 = (3 - x_2)/2 = 1.5`$; $`\partial f/\partial x_2 = x_1 + 2x_2 = 0`$ cho $`x_2 = -x_1/2 = -0.75`$. Nghiệm của hệ là $`(2, -1)`$. Các vòng tiếp theo: $`(1.875, -0.9375)`$, $`(1.96875, -0.984375)`$, …

</details>

**8.2.** Với $`f = \tfrac12 x^T H x + q^T x`$, phương pháp Gauss–Seidel:
A. Hội tụ với mọi $`H \succ 0`$ B. Chỉ hội tụ khi $`H`$ chéo trội C. Hội tụ sau $`n`$ vòng D. Hội tụ bậc hai

<details><summary>Đáp án</summary>

**A.** Mỗi bước tọa độ là cực tiểu chính xác của một hàm lồi chặt, nên $`f`$ giảm, và dãy hội tụ về nghiệm duy nhất. Chéo trội là điều kiện đủ, không phải điều kiện cần. Tốc độ là tuyến tính, không phải bậc hai; hội tụ sau hữu hạn vòng chỉ xảy ra trong trường hợp đặc biệt ($`H`$ chéo: 1 vòng).

</details>

### Nếu ra tự luận

"Trình bày giải thuật Gauss–Seidel tìm cực tiểu hàm $`n`$ biến, áp dụng cho ví dụ." Dàn ý: ý tưởng; các bước; công thức đóng cho hàm bậc hai; nhận xét hội tụ (vì sao zic-zac khi biến liên kết mạnh); bảng 3 vòng cho ví dụ, ví dụ câu 1 Đề số 1 ([lời giải](../../test1-paper1/)). Báo cáo đầy đủ: [câu 8](../../../solutions/q08-gauss-seidel/).

### Bẫy

- Dùng giá trị cũ thay vì giá trị vừa cập nhật (đó là phương pháp Jacobi, không phải Gauss–Seidel).
- Sai hệ số khi lập $`H`$ từ đề.
- Tiêu chuẩn dừng theo độ giảm của $`f`$ dừng quá sớm khi hội tụ chậm; dùng tiêu chuẩn theo $`x`$.

---

## Câu 9. Tìm kiếm trực tiếp (Hooke–Jeeves)

### Lý thuyết cốt lõi

- **Chỉ dùng giá trị hàm**, không cần đạo hàm, áp dụng được cho hàm không khả vi.
- **Bước dò** $`\mathrm{Dò}(y, \Delta)`$: với $`i = 1..n`$, nếu $`f(y + \Delta e_i) \lt f(y)`$ thì nhận điểm đó; nếu không thì thử $`y - \Delta e_i`$. Mỗi lần dò tốn từ $`n`$ đến $`2n`$ lần tính $`f`$.
- **Bước theo mẫu:** sau khi dò thành công từ $`x_b`$ đến $`x_n`$, nhảy tới $`x_p = x_n + (x_n - x_b) = 2x_n - x_b`$ rồi dò quanh $`x_p`$.
- Khi dò quanh điểm cơ sở thất bại thì giảm bước $`\Delta \leftarrow \beta \Delta`$ (thường $`\beta = 1/2`$). Dừng khi $`\Delta \lt \varepsilon`$.
- **Định lý:** nếu $`\nabla f`$ liên tục Lipschitz (hằng số $`L`$) và dò thất bại tại $`x`$ với bước $`\Delta`$, thì $`\Vert \nabla f(x) \Vert \le \sqrt n L \Delta / 2`$. Do đó khi $`\Delta \to 0`$, các điểm cơ sở tiến tới điểm dừng.
- **May mắn của lưới:** các điểm thử nằm trên lưới $`x^0 + \Delta \mathbb{Z}^n`$. Nếu $`x^*`$ nằm trên lưới thì thuật toán có thể "rơi trúng" nghiệm chính xác.

### Thuật toán

1. $`x_b = x^0`$, $`\Delta = \Delta_0`$.
2. $`x_n = \mathrm{Dò}(x_b, \Delta)`$. Nếu $`f(x_n) \lt f(x_b)`$ thì sang bước 3. Nếu không: dừng khi $`\Delta \lt \varepsilon`$, còn lại thì $`\Delta \leftarrow \beta\Delta`$ và lặp bước 2.
3. Bước mẫu: $`x_p = 2x_n - x_b`$, $`x_b \leftarrow x_n`$, $`x_n = \mathrm{Dò}(x_p, \Delta)`$. Nếu $`f(x_n) \lt f(x_b)`$ thì lặp bước 3; nếu không thì quay về bước 2.

### Dạng câu hỏi có thể gặp

| Dạng | Cách làm nhanh | Matlab |
|---|---|---|
| Điểm sau bước dò / bước mẫu đầu tiên | Tính tay 4–6 giá trị $`f`$ | `hooke_jeeves(f,x0,struct('delta',0.5))` |
| Điền $`x^*`$ | Tìm cực tiểu bằng `fminsearch` | `hooke_jeeves` |
| Lý thuyết | Không dùng đạo hàm; định lý chặn $`\Vert\nabla f\Vert`$ | — |

### Câu hỏi mẫu

**9.1.** Áp dụng Hooke–Jeeves cho $`F = (x_1 - 2)^4 + (x_1 - 2x_2)^2`$ từ $`x^0 = (0, 3)`$ với $`\Delta = 0.5`$. Điểm sau bước dò đầu tiên và điểm theo mẫu là:
A. $`(0.5, 2.5)`$ và $`(1, 2)`$ B. $`(0.5, 3)`$ và $`(1, 3)`$ C. $`(0.5, 2.5)`$ và $`(1.5, 2)`$ D. $`(0, 2.5)`$ và $`(0, 2)`$

<details><summary>Đáp án</summary>

**A.** $`F(0, 3) = 52`$. Thử $`x_1 + 0.5`$: $`F(0.5, 3) = 35.3125 \lt 52`$, nhận. Thử $`x_2 + 0.5`$: $`F(0.5, 3.5) = 47.3125`$, không nhận; thử $`x_2 - 0.5`$: $`F(0.5, 2.5) = 25.3125`$, nhận. Bước mẫu $`x_p = 2(0.5, 2.5) - (0, 3) = (1, 2)`$ với $`F = 10`$. Chạy tiếp, `hooke_jeeves` rơi trúng $`x^* = (2, 1)`$ (nằm trên lưới bước 0.5) sau 92 lần tính $`F`$.

</details>

**9.2.** Phát biểu nào **sai** về Hooke–Jeeves?
A. Không cần đạo hàm B. Mỗi bước dò tốn từ $`n`$ đến $`2n`$ lần tính $`f`$ C. Hội tụ bậc hai gần nghiệm D. Bước mẫu giúp đi nhanh dọc thung lũng

<details><summary>Đáp án</summary>

**C.** Phương pháp chỉ hội tụ tuyến tính và chậm (hàm bậc hai $`n = 10`$ cần 3409 lần tính $`f`$).

</details>

### Nếu ra tự luận

"Trình bày giải thuật tìm kiếm trực tiếp, áp dụng 1–2 vòng." Dàn ý: bước dò, bước mẫu, giảm bước (có sơ đồ); định lý chặn gradient kèm chứng minh Taylor ngắn; bảng tính tay như câu 9.1. Báo cáo đầy đủ: [câu 9](../../../solutions/q09-hooke-jeeves/).

### Bẫy

- Bước mẫu là $`2x_n - x_b`$, không phải $`x_n + \Delta`$.
- Quên quay về dò quanh $`x_b`$ khi dò quanh $`x_p`$ thất bại.
- Kết quả "chính xác tuyệt đối" có thể do may mắn của lưới, không chứng tỏ phương pháp nhanh.

---

## Câu 10. Gradient và gradient quán tính

### Lý thuyết cốt lõi

- $`-g`$ là **hướng giảm nhanh nhất**: theo Cauchy–Schwarz, $`\min_{\Vert d\Vert = 1} g^T d = -\Vert g\Vert`$, đạt tại $`d = -g/\Vert g\Vert`$.
- **Bước cố định** (hàm bậc hai): hội tụ khi và chỉ khi $`0 \lt \alpha \lt 2/L`$. Bước tốt nhất $`\alpha = 2/(L + \mu)`$ cho tỉ số hội tụ $`(\kappa - 1)/(\kappa + 1)`$.
- **Tìm kiếm chính xác** $`\alpha_k = \mathrm{argmin}_\alpha f(x_k - \alpha g_k)`$, với hàm bậc hai có công thức đóng:

```math
\alpha_k = \frac{g_k^T g_k}{g_k^T H g_k},\qquad g_{k+1}^T g_k = 0\ \ (\text{zic-zac}),\qquad
f(x_k) - f^* \le \Big(\frac{\kappa - 1}{\kappa + 1}\Big)^{2k} \big(f(x_0) - f^*\big).
```

- **Armijo:** thử $`\alpha = \alpha_0, \rho\alpha_0, \rho^2\alpha_0, \dots`$ đến khi $`f(x - \alpha g) \le f(x) - c\alpha\Vert g\Vert^2`$ (thường $`c = 10^{-4}`$, $`\rho = 1/2`$).
- **Quán tính (heavy-ball, Polyak):** $`x_{k+1} = x_k - \alpha g_k + \beta (x_k - x_{k-1})`$, với tham số tối ưu

```math
\alpha^* = \frac{4}{(\sqrt L + \sqrt \mu)^2},\qquad \beta^* = \Big(\frac{\sqrt L - \sqrt\mu}{\sqrt L + \sqrt\mu}\Big)^2,\qquad \text{tỉ số } \frac{\sqrt\kappa - 1}{\sqrt\kappa + 1}.
```

  Số bước tỉ lệ với $`\sqrt\kappa`$ thay vì $`\kappa`$. Với $`\beta = 0`$ thì đây chính là gradient bước cố định.

### Thuật toán

1. Cho $`x^0`$, $`\varepsilon`$ (quán tính thêm $`x^{-1} = x^0`$, $`\alpha`$, $`\beta`$).
2. Lặp: tính $`g_k`$; dừng nếu $`\Vert g_k\Vert \le \varepsilon`$; chọn $`\alpha_k`$; cập nhật $`x_{k+1}`$.

### Dạng câu hỏi có thể gặp

| Dạng | Cách làm nhanh | Matlab |
|---|---|---|
| Điền $`x^*`$ | Hàm bậc hai: `-H\q`; hàm khác: `fminsearch` | `grad_desc`, `heavy_ball` (toolkit), `gradient_descent`, `momentum_gradient` |
| Bước $`\alpha_0`$ tìm kiếm chính xác, điểm $`x_1`$ | $`\alpha = g^Tg / g^THg`$ | tính tay |
| Miền $`\alpha`$ hội tụ, $`\alpha`$ tối ưu, tỉ số | Từ trị riêng của $`H`$ | `eig(H)` |
| $`\alpha, \beta`$ tối ưu cho quán tính | Công thức trên | `eig(H)` |

### Câu hỏi mẫu

**10.1.** Với $`f = x_1^2 + 4x_2^2`$ và $`x^0 = (1, 1)`$, bước tìm kiếm chính xác đầu tiên của phương pháp gradient và điểm $`x^1`$ là?

<details><summary>Đáp án</summary>

$`g_0 = (2, 8)`$, $`H = \mathrm{diag}(2, 8)`$, nên $`\alpha_0 = \frac{4 + 64}{2\cdot 4 + 8\cdot 64} = \frac{68}{520} = \frac{17}{130} = 0.130769`$ và $`x^1 = (1, 1) - 0.130769\,(2, 8) = (0.738462,\ -0.046154)`$.

</details>

**10.2.** Với $`f = x_1^2 + 4x_2^2`$, phương pháp gradient bước cố định hội tụ khi và chỉ khi:
A. $`0 \lt \alpha \lt 0.5`$ B. $`0 \lt \alpha \lt 0.25`$ C. $`0 \lt \alpha \lt 0.2`$ D. $`0 \lt \alpha \lt 1`$

<details><summary>Đáp án</summary>

**B.** $`H = \mathrm{diag}(2, 8)`$ nên $`L = 8`$, và điều kiện là $`\alpha \lt 2/L = 0.25`$. Bước tốt nhất là $`2/(L + \mu) = 0.2`$, với tỉ số $`(\kappa - 1)/(\kappa + 1) = 3/5`$.

</details>

**10.3.** Cũng hàm trên, tham số tối ưu của gradient quán tính là:
A. $`\alpha = 2/9, \beta = 1/9`$ B. $`\alpha = 0.2, \beta = 0.6`$ C. $`\alpha = 1/4, \beta = 1/3`$ D. $`\alpha = 2/9, \beta = 1/3`$

<details><summary>Đáp án</summary>

**A.** $`\alpha = 4/(\sqrt8 + \sqrt2)^2 = 4/18 = 2/9`$ và $`\beta = \big(\frac{\sqrt8 - \sqrt2}{\sqrt8 + \sqrt2}\big)^2 = (1/3)^2 = 1/9`$. Tỉ số hội tụ $`(\sqrt\kappa - 1)/(\sqrt\kappa + 1) = 1/3`$. Từ $`(1, 1)`$ đến $`\Vert g\Vert \le 10^{-6}`$: gradient ($`\alpha = 0.2`$) cần 32 bước, quán tính cần 18 bước.

</details>

### Nếu ra tự luận

"Trình bày thuật toán gradient và gradient quán tính, so sánh, áp dụng." Dàn ý:
1. Chứng minh $`-g`$ là hướng giảm nhanh nhất.
2. Ba cách chọn bước; với hàm bậc hai: công thức $`\alpha_k`$, tính vuông góc $`g_{k+1} \perp g_k`$, tỉ số hội tụ.
3. Quán tính: công thức, tham số tối ưu, vì sao giảm zic-zac.
4. Ví dụ: 2–3 bước tính tay như câu 10.1, và bảng so sánh số bước.

Báo cáo đầy đủ: [câu 10](../../../solutions/q10-gradient-momentum/).

### Bẫy

- $`\alpha \ge 2/L`$ làm phương pháp phân kỳ.
- Tìm kiếm chính xác nhanh hay chậm tùy điểm xuất phát (câu 6 Đề số 1 chỉ 5 bước, câu 6 đề luyện 441 bước).
- Với $`\varepsilon`$ áp lên $`\Vert g \Vert`$, sai số theo $`x`$ có thể lớn hơn $`\varepsilon`$ nhiều khi $`\kappa`$ lớn.

---

## Câu 11. Newton

### Lý thuyết cốt lõi

- Cực tiểu xấp xỉ Taylor bậc hai: $`m(d) = f + g^T d + \tfrac12 d^T H d`$, nên $`H_k d_k = -g_k`$.
- **Hàm bậc hai** với $`H \succ 0`$: hội tụ trong **đúng 1 bước** từ mọi $`x^0`$.
- **Hội tụ bậc hai** gần $`x^*`$ khi $`\nabla^2 f(x^*) \succ 0`$ và Hessian liên tục Lipschitz: $`\Vert x_{k+1} - x^*\Vert \le C\Vert x_k - x^*\Vert^2`$. Số chữ số đúng tăng gấp đôi mỗi bước.
- **Hessian suy biến tại $`x^*`$:** chỉ hội tụ tuyến tính. Ví dụ $`f = x^4`$: $`x_{k+1} = x_k - \frac{4x_k^3}{12x_k^2} = \tfrac23 x_k`$. Hàm Bazaraa cũng có tỉ số $`2/3`$.
- **$`H_k`$ không xác định dương:** $`d_k`$ có thể không phải hướng giảm. Hai cách khắc phục: Newton giảm bước (Armijo theo hướng $`d_k`$, chỉ dùng được khi $`g^T d \lt 0`$) và Newton cải biên ($`H_k + \tau I \succ 0`$, tăng $`\tau`$ đến khi phân tích Cholesky thành công).
- Bất biến với phép đổi biến affine, nên không bị ảnh hưởng bởi $`\kappa`$.

### Thuật toán

1. Cho $`x^0`$, $`\varepsilon`$.
2. Lặp: nếu $`\Vert g_k\Vert \le \varepsilon`$ thì dừng; giải $`H_k d_k = -g_k`$ (cải biên: thay $`H_k`$ bằng $`H_k + \tau I`$ nếu cần); chọn $`\alpha_k`$ (thuần túy: $`\alpha_k = 1`$; giảm bước: Armijo); $`x_{k+1} = x_k + \alpha_k d_k`$.

### Dạng câu hỏi có thể gặp

| Dạng | Cách làm nhanh | Matlab |
|---|---|---|
| Điền $`x^*`$ | Kiểm tra lồi chặt, khử biến tuyến tính, `fminunc` | `newton(g,h,x0,1e-4)` (toolkit), `newton_method` |
| Một bước Newton từ $`x^0`$ | Giải hệ $`2\times 2`$ hoặc $`3\times 3`$ | `x0 - h(x0)\g(x0)` |
| Tốc độ hội tụ, số bước với hàm bậc hai | Lý thuyết ở trên | — |
| $`H`$ không xác định dương | `eig(H)`, kiểm tra $`g^Td`$ | `newton_method(..., struct('variant','modified'))` |

### Câu hỏi mẫu

**11.1.** Áp dụng Newton cho $`f(x) = x^4`$ từ $`x^0 = 1`$. Dãy lặp và tốc độ hội tụ là:
A. $`x_k = (2/3)^k`$, tuyến tính B. $`x_k = (1/2)^k`$, tuyến tính C. $`x_1 = 0`$ D. hội tụ bậc hai

<details><summary>Đáp án</summary>

**A.** $`x_{k+1} = x_k - 4x_k^3/(12x_k^2) = \tfrac23 x_k`$, cho $`0.6667, 0.4444, 0.2963, \dots`$ Vì $`f''(0) = 0`$ (Hessian suy biến tại nghiệm) nên chỉ hội tụ tuyến tính.

</details>

**11.2.** Với $`f = x_1^4 + x_1x_2 + (1 + x_2)^2`$ và $`x^0 = (0, 0)`$, bước Newton thuần túy cho kết quả gì?

<details><summary>Đáp án</summary>

$`g(0) = (0, 2)`$ và $`H(0) = [0\ 1;\ 1\ 2]`$, với trị riêng $`-0.414`$ và $`2.414`$, nên $`H`$ **không xác định dương**. Giải $`H d = -g`$ được $`d = (-2, 0)`$, $`g^T d = 0`$ (không phải hướng giảm), và $`f(x^0 + d) = 17 \gt f(x^0) = 1`$. Phải dùng Newton cải biên: 6 bước đến $`x^* = (0.695884, -1.347942)`$, $`f^* = -0.582445`$ (khớp với `fminsearch`).

</details>

**11.3.** Với $`f = 2x_1^2 + x_1x_2 + x_2^2 - x_1 + 2x_2`$, một bước Newton từ $`(5, -7)`$ cho:
A. $`(4/7, -9/7)`$ B. $`(0, 0)`$ C. $`(1, -1)`$ D. Cần nhiều bước

<details><summary>Đáp án</summary>

**A.** Hàm bậc hai với $`H = [4\ 1;\ 1\ 2] \succ 0`$, nên một bước cho đúng $`x^* = -H^{-1}q = (0.571429, -1.285714) = (4/7, -9/7)`$.

</details>

### Nếu ra tự luận

"Trình bày thuật toán Newton, các cải tiến, áp dụng." Dàn ý:
1. Dẫn xuất từ Taylor bậc hai.
2. Các bước thuật toán.
3. Định lý hội tụ bậc hai (phát biểu) và trường hợp hàm bậc hai.
4. Ba nhược điểm và hai cải tiến (giảm bước, cải biên).
5. Ví dụ: bảng $`\Vert g\Vert`$ cho thấy hội tụ bậc hai (câu 3 Đề số 1: $`0.43 \to 1.4\cdot10^{-3} \to 4.5\cdot 10^{-11}`$).

Báo cáo đầy đủ: [câu 11](../../../solutions/q11-newton-method/).

### Bẫy

- Quên kiểm tra $`H \succ 0`$ trước khi kết luận bước Newton là hướng giảm.
- Lập Hessian sai (thiếu hệ số 2 ở số hạng bình phương, hoặc thiếu $`12x_1^2`$).
- Viết `inv(H)*g` thay vì `H\g`: kết quả vẫn đúng, nhưng chậm hơn và kém ổn định.

---

## Câu 12. Tựa Newton

### Lý thuyết cốt lõi

- $`B_k \approx H^{-1}`$, hướng $`d_k = -B_k g_k`$, tìm kiếm theo tia (Wolfe hoặc chính xác). Đặt $`s_k = x_{k+1} - x_k`$ và $`y_k = g_{k+1} - g_k`$.
- **Phương trình cát tuyến:** $`B_{k+1} y_k = s_k`$ (vì $`y_k \approx H s_k`$).
- Công thức cập nhật, với $`\rho_k = 1/(y_k^T s_k)`$:

```math
\text{BFGS:}\ \ B_{k+1} = (I - \rho_k s_k y_k^T)\, B_k\, (I - \rho_k y_k s_k^T) + \rho_k s_k s_k^T,\qquad
\text{DFP:}\ \ B_{k+1} = B_k - \frac{B_k y_k y_k^T B_k}{y_k^T B_k y_k} + \frac{s_k s_k^T}{y_k^T s_k}.
```

- Cả hai giữ $`B_{k+1} \succ 0`$ nếu $`B_k \succ 0`$ và $`y_k^T s_k \gt 0`$ (điều kiện độ cong). Điều kiện Wolfe mạnh bảo đảm $`y_k^Ts_k \gt 0`$.
- **Hàm bậc hai và tìm kiếm chính xác:** các hướng liên hợp với nhau theo $`H`$, hội tụ sau tối đa $`n`$ bước và $`B_n = H^{-1}`$.
- Hội tụ **siêu tuyến tính**. BFGS tự điều chỉnh tốt hơn DFP (Rosenbrock: BFGS 39 bước, DFP 840 bước).
- Thường lấy $`B_0 = I`$, nên bước đầu là một bước gradient. Có thể chuẩn hóa $`B_0 = \frac{s^Ty}{y^Ty} I`$.

### Thuật toán

1. $`x^0`$, $`B_0 = I`$.
2. Lặp: dừng nếu $`\Vert g_k\Vert \le \varepsilon`$; $`d_k = -B_k g_k`$; chọn $`\alpha_k`$ (Wolfe/Armijo/chính xác); $`x_{k+1} = x_k + \alpha_k d_k`$; tính $`s_k, y_k`$; nếu $`y_k^Ts_k \gt 0`$ thì cập nhật $`B_{k+1}`$, ngược lại giữ nguyên $`B`$.

### Dạng câu hỏi có thể gặp

| Dạng | Cách làm nhanh | Matlab |
|---|---|---|
| Điền $`x^*`$ | Lồi chặt (bậc bốn + bậc hai xác định dương), `fminunc`; khử biến tuyến tính | `bfgs(f,g,x0,1e-4)` (toolkit), `quasi_newton` (`'bfgs'`/`'dfp'`) |
| Tính $`B_1`$ từ $`s, y`$ | Thay vào công thức, kiểm tra $`B_1 y = s`$ | tính tay / `mat2str` |
| Số bước với hàm bậc hai | $`\le n`$ (với tìm kiếm chính xác) | — |

### Câu hỏi mẫu

**12.1.** Với $`B_0 = I`$, $`s_0 = (1, 0)`$, $`y_0 = (2, 1)`$, ma trận $`B_1`$ theo BFGS là:
A. $`[0.75\ -0.5;\ -0.5\ 1]`$ B. $`[0.7\ -0.4;\ -0.4\ 0.8]`$ C. $`[0.5\ 0;\ 0\ 1]`$ D. $`[1\ 0;\ 0\ 1]`$

<details><summary>Đáp án</summary>

**A.** $`\rho = 1/(y^Ts) = 1/2`$, $`B_1 = (I - \tfrac12 s y^T)(I - \tfrac12 y s^T) + \tfrac12 ss^T = [0.75\ -0.5;\ -0.5\ 1]`$. Kiểm tra cát tuyến: $`B_1 y = (1.5 - 0.5,\ -1 + 1) = (1, 0) = s`$. Phương án B là kết quả của DFP (cũng thỏa $`B_1 y = s`$).

</details>

**12.2.** BFGS với tìm kiếm chính xác cho $`f = \tfrac12 x^TAx - b^Tx`$, $`A \succ 0`$ cỡ $`3 \times 3`$, hội tụ sau:
A. 1 bước B. tối đa 3 bước C. vô hạn bước D. tối đa 6 bước

<details><summary>Đáp án</summary>

**B.** Với $`A = [2\ 0\ 0;\ 0\ 4\ 1;\ 0\ 1\ 3]`$, $`b = (1, 1, 1)`$: `quasi_newton` hội tụ sau 3 bước và $`\Vert B_3 - A^{-1}\Vert = 1.5\cdot 10^{-15}`$.

</details>

### Nếu ra tự luận

"Trình bày thuật toán tựa Newton (BFGS, DFP), áp dụng." Dàn ý:
1. Động cơ: tránh tính và giải với Hessian.
2. Phương trình cát tuyến.
3. Công thức BFGS/DFP; chứng minh $`B_{k+1} y_k = s_k`$ (thay trực tiếp); điều kiện giữ xác định dương.
4. Tính chất với hàm bậc hai.
5. Ví dụ: 2 bước tính tay ($`d_0 = -g_0`$, $`\alpha_0`$, $`s_0`$, $`y_0`$, $`B_1`$, $`d_1`$).

Báo cáo đầy đủ: [câu 12](../../../solutions/q12-quasi-newton/).

### Bẫy

- Nhầm $`B`$ (xấp xỉ $`H^{-1}`$) với xấp xỉ $`H`$: hai dạng BFGS có công thức khác nhau.
- Cập nhật khi $`y^Ts \le 0`$ làm mất tính xác định dương.
- Kết quả "xấu bất thường" (DFP 840 bước) là hiện tượng thật của DFP, không phải lỗi chương trình.
