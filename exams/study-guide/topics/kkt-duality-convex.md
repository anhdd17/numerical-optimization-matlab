# Điều kiện KKT, đối ngẫu Lagrange, quy hoạch lồi (câu 4, 5, 13)

[← Tài liệu ôn thi](../README.md)

Ba câu này là phần lý thuyết của môn và là nguồn chính của **tự luận** (Phần B của [Đề số 1](../../test1-paper1/README.md#phần-b) thuộc đúng nhóm này). Báo cáo đầy đủ: [câu 4](../../../solutions/q04-kkt-conditions/), [câu 5](../../../solutions/q05-lagrange-duality/), [câu 13](../../../solutions/q13-convex-programming/). Mọi số liệu dưới đây lấy từ output của [`check_kkt_duality.m`](../check_kkt_duality.m) (chạy trên GNU Octave).

Bài toán chung trong cả file:

```math
(P):\quad \min f(x)\quad \text{s.t.}\quad g_i(x) \le 0\ (i = 1..m),\qquad h_j(x) = 0\ (j = 1..p),
```

với $`I(x) = \{i : g_i(x) = 0\}`$ là tập ràng buộc tích cực và hàm Lagrange $`\mathcal{L}(x, \lambda, \mu) = f(x) + \sum_i \lambda_i g_i(x) + \sum_j \mu_j h_j(x)`$.

---

## Câu 4. Điều kiện chính quy và KKT

### Lý thuyết cốt lõi

**Nón.** Nón tiếp tuyến $`T(x^*)`$ gồm các hướng đi được trong $`D`$ (theo nghĩa giới hạn). Nón tuyến tính hóa là

```math
L(x^*) = \{ d : \nabla g_i(x^*)^T d \le 0\ (i \in I(x^*)),\ \ \nabla h_j(x^*)^T d = 0 \}.
```

Luôn có $`T \subseteq L`$. Nếu $`x^*`$ là cực tiểu địa phương thì $`\nabla f(x^*)^T d \ge 0`$ với mọi $`d \in T(x^*)`$.

**Điều kiện chính quy** (bảo đảm thay được $`T`$ bằng $`L`$):

| Tên | Phát biểu tại $`x^*`$ |
|---|---|
| LICQ | Các $`\nabla g_i(x^*)`$ ($`i \in I`$) và $`\nabla h_j(x^*)`$ độc lập tuyến tính |
| MFCQ | Các $`\nabla h_j`$ độc lập tuyến tính và $`\exists d`$: $`\nabla h_j^T d = 0`$, $`\nabla g_i^T d \lt 0`$ ($`i \in I`$) |
| Slater (bài toán lồi) | $`\exists \hat x`$: $`g_i(\hat x) \lt 0`$ mọi $`i`$, $`h_j(\hat x) = 0`$ |
| Ràng buộc tuyến tính | Mọi $`g_i, h_j`$ affine |
| Abadie | $`T(x^*) = L(x^*)`$ |

```math
\text{LICQ} \Rightarrow \text{MFCQ} \Rightarrow \text{Abadie},\qquad \text{Slater} \Rightarrow \text{MFCQ tại mọi } x \in D,\qquad \text{affine} \Rightarrow \text{Abadie}.
```

**Điều kiện KKT** tại $`x^*`$: tồn tại $`\lambda, \mu`$ với

```math
\begin{aligned}
&\nabla f(x^*) + \textstyle\sum_i \lambda_i \nabla g_i(x^*) + \sum_j \mu_j \nabla h_j(x^*) = 0 && \text{(dừng)}\\
&g_i(x^*) \le 0,\quad h_j(x^*) = 0 && \text{(chấp nhận được)}\\
&\lambda_i \ge 0\ \ (\mu_j \text{ tùy ý dấu}) && \text{(dấu nhân tử)}\\
&\lambda_i\, g_i(x^*) = 0 && \text{(bù)}
\end{aligned}
```

**Định lý KKT (điều kiện cần).** Nếu $`x^*`$ là cực tiểu địa phương và tại đó có một điều kiện chính quy thì $`x^*`$ là điểm KKT. *Ý chứng minh:* hệ $`\nabla g_i^T d \le 0\ (i \in I)`$, $`\nabla h_j^T d = 0`$, $`\nabla f^T d \lt 0`$ vô nghiệm (do $`T = L`$), nên theo bổ đề Farkas, $`-\nabla f`$ thuộc nón sinh bởi $`\nabla g_i`$ ($`i \in I`$, hệ số $`\ge 0`$) cộng không gian sinh bởi $`\nabla h_j`$.

**Tính chất nhân tử.** LICQ suy ra nhân tử duy nhất. Tại cực tiểu, MFCQ khi và chỉ khi tập nhân tử khác rỗng và bị chặn (Gauvin).
**Fritz John:** không cần chính quy, luôn có $`(\lambda_0, \lambda, \mu) \ne 0`$, $`\lambda_0 \ge 0`$; chính quy bảo đảm $`\lambda_0 \gt 0`$.

**Điều kiện đủ.** (i) Bài toán lồi thì mọi điểm KKT là cực tiểu toàn cục (chứng minh ở [câu 13](#câu-13-quy-hoạch-lồi-chính-quy-và-kkt)). (ii) Cấp 2: KKT và $`d^T \nabla^2_{xx}\mathcal{L}\, d \gt 0`$ trên nón tới hạn suy ra cực tiểu địa phương chặt.
**Ý nghĩa độ nhạy:** nới $`g_i \le 0`$ thành $`g_i \le \varepsilon`$ thì $`f^*`$ giảm khoảng $`\lambda_i \varepsilon`$.

### Thuật toán

**Kiểm tra chính quy tại một điểm $`\bar x`$:**
1. Kiểm tra $`\bar x`$ chấp nhận được và tìm $`I(\bar x)`$.
2. Lập ma trận các gradient tích cực (cùng $`\nabla h_j`$). Nếu hạng bằng số vector thì LICQ thỏa.
3. Nếu LICQ không thỏa: tìm $`d`$ với $`\nabla h_j^T d = 0`$, $`\nabla g_i^T d \lt 0`$ (thử $`d`$ đơn giản như $`(1,1)`$). Có $`d`$ thì MFCQ thỏa.

**Giải KKT theo trường hợp tập tích cực:**
1. Với mỗi tập con $`A \subseteq \{1..m\}`$ (ưu tiên trường hợp dễ: $`A = \emptyset`$ trước): đặt $`g_i = 0`$ cho $`i \in A`$, $`\lambda_i = 0`$ cho $`i \notin A`$.
2. Giải hệ dừng cùng các ràng buộc tích cực và $`h = 0`$.
3. **Loại** nghiệm vi phạm $`g_i \le 0`$ (với $`i \notin A`$) hoặc có $`\lambda_i \lt 0`$.
4. Điểm còn lại là điểm KKT. Kết luận cực tiểu nhờ tính lồi, điều kiện cấp 2, hoặc so sánh giá trị $`f`$.

Mẹo hình học: với $`f = \Vert x - p \Vert^2`$ (chiếu điểm $`p`$), thử chiếu lên từng ràng buộc. Hình chiếu nào thỏa các ràng buộc còn lại thì đó là nghiệm.

### Dạng câu hỏi có thể gặp

| Dạng câu | Cách làm nhanh | Lệnh Matlab |
|---|---|---|
| "$`\bar x`$ có phải điểm KKT?" | Tìm $`I(\bar x)`$, giải hệ dừng, xét dấu $`\lambda`$ | `r = kkt_check(prob, x)` (`.KKT`, `.lambda`) |
| "LICQ / MFCQ thỏa tại $`\bar x`$?" | Hạng các gradient tích cực; tìm $`d`$ | `r.LICQ`, `r.MFCQ` |
| Điền nhân tử $`\lambda, \mu`$ | Giải hệ tuyến tính $`\sum \lambda_i \nabla g_i = -\nabla f`$ | `[grad g...] \ (-gradf)` |
| Điền $`x^*`$ của bài có ràng buộc | Giải KKT theo trường hợp, đối chiếu solver | `solve_nlp(prob, x0)` |
| "Cực tiểu nhưng không là điểm KKT" | Tìm chỗ hỏng chính quy (gradient tích cực cùng phương) | `r.stationarity` |

`kkt_check` dùng `lsqnonneg`: tại điểm **không** KKT nó in nhân tử không âm *gần nhất*, không phải nghiệm của hệ dừng. Muốn chỉ ra nhân tử âm phải giải hệ dừng bằng tay hoặc bằng `\`.

### Câu hỏi mẫu

**Câu 4.1.** Bài toán $`\min (x_1 - 3)^2 + (x_2 - 3)^2`$ với $`g_1 = x_1 + x_2 - 4 \le 0`$, $`g_2 = -x_1 \le 0`$. Tại $`\bar x = (0, 4)`$:
A. Là điểm KKT với $`\lambda = (2, 0)`$. B. Là điểm KKT với nhân tử không duy nhất. C. LICQ thỏa nhưng không là điểm KKT. D. Không chấp nhận được.

<details><summary>Đáp án</summary>

**C.** $`g_1(0, 4) = 0`$ và $`g_2 = 0`$ nên cả hai tích cực. $`\nabla g_1 = (1, 1)`$ và $`\nabla g_2 = (-1, 0)`$ độc lập tuyến tính, nên LICQ thỏa. $`\nabla f(0, 4) = (-6, 2)`$. Hệ dừng $`-6 + \lambda_1 - \lambda_2 = 0`$, $`2 + \lambda_1 = 0`$ cho $`\lambda = (-2, -8)`$, có nhân tử âm nên không phải điểm KKT. Điểm KKT thật là $`(2, 2)`$ với $`\lambda = (2, 0)`$, $`f^* = 2`$ (hình chiếu của $`(3,3)`$ lên $`x_1 + x_2 = 4`$). `kkt_check` cho `KKT = 1, lambda = [2 0]` tại $`(2,2)`$; `solve_nlp` cho $`x = (2, 2)`$, $`f = 2`$.

</details>

**Câu 4.2.** Ràng buộc $`-x_1 \le 0`$, $`-x_2 \le 0`$, $`-x_1 - 2x_2 \le 0`$, xét tại $`(0, 0)`$, hàm mục tiêu $`f = x_1 + x_2`$:
A. LICQ thỏa. B. MFCQ thỏa, LICQ không thỏa, nhân tử không duy nhất. C. Không thỏa điều kiện chính quy nào, không là điểm KKT. D. LICQ thỏa, nhân tử duy nhất.

<details><summary>Đáp án</summary>

**B.** Ba ràng buộc đều tích cực, có 3 gradient $`(-1, 0), (0, -1), (-1, -2)`$ trong $`\mathbb{R}^2`$ nên phụ thuộc tuyến tính, LICQ không thỏa. Với $`d = (1, 1)`$, các tích vô hướng là $`-1, -1, -3`$, đều âm, nên MFCQ thỏa. Hệ dừng $`(1, 1) = \lambda_1(1, 0) + \lambda_2(0, 1) + \lambda_3(1, 2)`$ cho $`\lambda = (1 - t,\ 1 - 2t,\ t)`$ với $`t \in [0, 1/2]`$: tập nhân tử bị chặn nhưng không duy nhất (đúng định lý Gauvin). `kkt_check`: `LICQ = 0, MFCQ = 1, KKT = 1`, nhân tử duy nhất `= 0`.

</details>

**Câu 4.3.** Bài toán $`\min x_1^2 + x_2^2 + x_3^2`$ với $`x_1 + 2x_2 + 3x_3 = 14`$. Điền $`x^*`$, $`\mu^*`$, $`f^*`$.

<details><summary>Đáp án</summary>

Điều kiện dừng $`2x + \mu(1, 2, 3) = 0`$ cho $`x = -\tfrac{\mu}{2}(1, 2, 3)`$. Thế vào ràng buộc: $`-\tfrac{\mu}{2}(1 + 4 + 9) = 14`$, nên $`\mu = -2`$. Vậy $`x^* = (1, 2, 3)`$, $`f^* = 14`$. Ràng buộc đẳng thức nên $`\mu \lt 0`$ là hợp lệ. Bài toán lồi nên đây là cực tiểu toàn cục. `kkt_check`: `KKT = 1, mu = -2`.

</details>

### Nếu ra tự luận

**Tự luận 4.1 — Giải bài toán bằng điều kiện KKT.**
*Đề:* $`\min f = (x_1 - 2)^2 + (x_2 - 3)^2`$ với $`g_1 = x_1 + x_2 - 3 \le 0`$, $`g_2 = x_1 - 2x_2 \le 0`$. Kiểm tra điều kiện chính quy, viết KKT và giải.

*Lời giải.*
1. **Tính lồi, chính quy.** $`\nabla^2 f = 2I \succ 0`$ nên $`f`$ lồi chặt; $`g_1, g_2`$ affine. Mọi ràng buộc affine nên điều kiện chính quy luôn thỏa (Abadie), KKT là điều kiện cần. Bài toán lồi nên KKT cũng là điều kiện đủ. Nghiệm duy nhất vì $`f`$ lồi chặt. (Slater cũng thỏa: $`\hat x = (0, 1)`$ cho $`g = (-2, -2)`$.)
2. **KKT:** $`2(x_1 - 2) + \lambda_1 + \lambda_2 = 0`$, $`2(x_2 - 3) + \lambda_1 - 2\lambda_2 = 0`$, $`g_i \le 0`$, $`\lambda_i \ge 0`$, $`\lambda_i g_i = 0`$.
3. **Xét trường hợp:**

| Tích cực | Giải | Kiểm tra | Kết luận |
|---|---|---|---|
| Không có | $`x = (2, 3)`$ | $`g_1 = 2 \gt 0`$ | loại |
| Chỉ $`g_2`$ | hình chiếu lên $`x_1 = 2x_2`$: $`x = (2.8, 1.4)`$ | $`g_1 = 1.2 \gt 0`$ | loại |
| Cả hai | $`x = (2, 1)`$, $`\lambda = (4/3, -4/3)`$ | $`\lambda_2 \lt 0`$ | loại |
| Chỉ $`g_1`$ | $`x_1 - 2 = x_2 - 3`$, $`x_1 + x_2 = 3`$ cho $`x = (1, 2)`$, $`\lambda_1 = 2`$ | $`g_2 = -3 \le 0`$, $`\lambda_1 \ge 0`$ | **điểm KKT** |

4. **Kết luận:** $`x^* = (1, 2)`$, $`\lambda^* = (2, 0)`$, $`f^* = 2`$. Tại $`x^*`$ chỉ $`g_1`$ tích cực với $`\nabla g_1 = (1,1) \ne 0`$ nên LICQ thỏa và nhân tử duy nhất. Ý nghĩa: nới thành $`x_1 + x_2 \le 3 + \varepsilon`$ thì $`f^*`$ giảm khoảng $`2\varepsilon`$.
*Kiểm chứng:* `kkt_check` tại $`(1,2)`$ cho `KKT = 1, lambda = [2 0]`; trường hợp cả hai tích cực cho $`\lambda = (1.333, -1.333)`$; `solve_nlp` cho $`x = (1, 2)`$, $`f = 2`$.

**Tự luận 4.2 — Phát biểu điều kiện chính quy; ví dụ KKT không cần khi thiếu chính quy.**
*Lời giải mẫu.*
1. Phát biểu bài toán (P), tập tích cực, các điều kiện LICQ, MFCQ, Slater, Abadie và quan hệ giữa chúng (bảng và sơ đồ suy ra ở trên). Nêu chứng minh LICQ ⇒ MFCQ: hệ $`\nabla g_i^T d = -1`$ ($`i \in I`$), $`\nabla h_j^T d = 0`$ có ma trận hệ số với các hàng độc lập tuyến tính nên có nghiệm $`d`$. Và Slater ⇒ MFCQ: với $`d = \hat x - x^*`$, tính lồi cho $`\nabla g_i(x^*)^T d \le g_i(\hat x) - g_i(x^*) = g_i(\hat x) \lt 0`$ với $`i \in I`$.
2. Phát biểu định lý KKT (điều kiện cần với chính quy) và ý chứng minh qua Farkas.
3. **Ví dụ** (câu 4, ví dụ 5): $`\min x_1`$ với $`g_1 = x_2 - x_1^3 \le 0`$, $`g_2 = -x_2 \le 0`$. Từ $`0 \le x_2 \le x_1^3`$ suy ra $`x_1 \ge 0`$, nên $`x^* = (0, 0)`$ là cực tiểu toàn cục. Tại đó $`\nabla g_1 = (0, 1)`$ và $`\nabla g_2 = (0, -1)`$ phụ thuộc tuyến tính, nên LICQ hỏng. Không có $`d`$ với $`d_2 \lt 0`$ và $`-d_2 \lt 0`$, nên MFCQ hỏng. $`T = \{d_1 \ge 0, d_2 = 0\} \ne L = \{d_2 = 0\}`$, nên Abadie hỏng. Hệ dừng $`(1, 0) + \lambda_1 (0, 1) + \lambda_2 (0, -1) = 0`$ có thành phần thứ nhất $`1 = 0`$, vô nghiệm: **cực tiểu nhưng không phải điểm KKT**. Chỉ có nhân tử Fritz John với $`\lambda_0 = 0`$.
4. Ví dụ ngược lại: điểm KKT không phải cực tiểu (câu 4, ví dụ 1: $`\min x_1 + x_2`$ trên đường tròn $`x_1^2 + x_2^2 = 2`$ có điểm KKT $`(1,1)`$ với $`\mu = -\tfrac12`$ là cực đại). Vậy KKT chỉ đủ khi có tính lồi hoặc điều kiện cấp 2.

### Bẫy

- Nhân tử của ràng buộc **đẳng thức** có thể âm; chỉ $`\lambda_i`$ của bất đẳng thức cần $`\ge 0`$.
- Phải đưa ràng buộc về dạng $`g \le 0`$ trước khi xét dấu ($`x_1 + x_2 \ge 1`$ thành $`1 - x_1 - x_2 \le 0`$).
- "Hệ dừng có nghiệm" chưa đủ: còn phải kiểm tra chấp nhận được, dấu $`\lambda`$ và điều kiện bù.
- Điểm KKT của bài toán **không lồi** có thể là cực đại hoặc điểm yên ngựa.
- `kkt_check` tại điểm không KKT in nhân tử không âm gần nhất. Không được chép số đó làm "nhân tử".

---

## Câu 5. Đối ngẫu Lagrange

### Lý thuyết cốt lõi

**Hàm đối ngẫu** $`q(\lambda, \mu) = \inf_{x \in X} \mathcal{L}(x, \lambda, \mu)`$ (có thể bằng $`-\infty`$). **Bài toán đối ngẫu** $`(D)`$: $`\max q(\lambda, \mu)`$ với $`\lambda \ge 0`$, giá trị tối ưu ký hiệu $`d^*`$.
- $`q`$ luôn **lõm** (inf của các hàm affine theo $`(\lambda, \mu)`$), nên $`(D)`$ luôn là bài toán lồi, kể cả khi $`(P)`$ không lồi.
- **Đối ngẫu yếu:** $`d^* \le p^*`$. *Chứng minh:* với $`\bar x`$ chấp nhận được và $`\lambda \ge 0`$, ta có $`q(\lambda, \mu) \le \mathcal{L}(\bar x, \lambda, \mu) = f(\bar x) + \sum \lambda_i g_i(\bar x) + \sum \mu_j h_j(\bar x) \le f(\bar x)`$. Lấy inf theo $`\bar x`$ rồi sup theo $`(\lambda, \mu)`$.
- **Hệ quả:** nếu $`f(\bar x) = q(\bar\lambda, \bar\mu)`$ thì cả hai tối ưu. Nếu $`p^* = -\infty`$ thì $`(D)`$ vô nghiệm.
- **Đối ngẫu mạnh (Slater):** nếu $`(P)`$ lồi, thỏa Slater và $`p^*`$ hữu hạn thì $`d^* = p^*`$ và $`(D)`$ đạt nghiệm. Chứng minh dùng định lý tách tập lồi $`\{(u, v, t) : g(x) \le u, h(x) = v, f(x) \le t\}`$; Slater loại siêu phẳng thẳng đứng.
- **Điểm yên ngựa:** $`\mathcal{L}(x^*, \lambda, \mu) \le \mathcal{L}(x^*, \lambda^*, \mu^*) \le \mathcal{L}(x, \lambda^*, \mu^*)`$ đúng khi và chỉ khi $`x^*`$ tối ưu, $`(\lambda^*, \mu^*)`$ tối ưu đối ngẫu, và $`p^* = d^*`$.
- **Liên hệ KKT:** khi đối ngẫu mạnh, $`\lambda_i^* g_i(x^*) = 0`$ và $`x^*`$ cực tiểu hóa $`\mathcal{L}(\cdot, \lambda^*, \mu^*)`$. Nghiệm đối ngẫu chính là nhân tử KKT.

**Công thức đóng cần thuộc:**

```math
\text{LP chính tắc: } \min c^Tx,\ Ax = b,\ x \ge 0 \ \longleftrightarrow\ \max b^Ty,\ A^Ty \le c\ \ (y \text{ tự do})
```

```math
\text{QP: } \min \tfrac12 x^TQx + c^Tx,\ Ax \le b\ \ (Q \succ 0):\quad x(\lambda) = -Q^{-1}(c + A^T\lambda),\quad q(\lambda) = -\tfrac12 (c + A^T\lambda)^T Q^{-1}(c + A^T\lambda) - b^T\lambda
```

**Bảng quy tắc lập đối ngẫu LP** (đọc từ trái sang phải cho bài min, từ phải sang trái cho bài max):

| Bài toán $`\min c^Tx`$ | Bài toán $`\max b^Ty`$ |
|---|---|
| Ràng buộc $`i`$: $`a_i^Tx \ge b_i`$ | $`y_i \ge 0`$ |
| Ràng buộc $`i`$: $`a_i^Tx \le b_i`$ | $`y_i \le 0`$ |
| Ràng buộc $`i`$: $`a_i^Tx = b_i`$ | $`y_i`$ tự do |
| Biến $`x_j \ge 0`$ | Ràng buộc $`j`$: $`A_j^Ty \le c_j`$ |
| Biến $`x_j \le 0`$ | Ràng buộc $`j`$: $`A_j^Ty \ge c_j`$ |
| Biến $`x_j`$ tự do | Ràng buộc $`j`$: $`A_j^Ty = c_j`$ |

**Định lý đối ngẫu LP:** một bài có nghiệm thì bài kia cũng có, và $`c^Tx^* = b^Ty^*`$ (không cần Slater). Một bài không bị chặn thì bài kia vô nghiệm. **Bù:** $`x_j (c - A^Ty)_j = 0`$ và $`y_i (a_i^Tx - b_i) = 0`$. Trong thuật toán đơn hình, $`y = B^{-T}c_B`$ và dấu hiệu tối ưu chính là $`y`$ chấp nhận được đối ngẫu; $`y_i`$ là giá bóng.

### Thuật toán

**Lập và giải đối ngẫu (bài phi tuyến):**
1. Đưa về dạng $`g_i \le 0`$, $`h_j = 0`$; viết $`\mathcal{L}`$.
2. Cực tiểu $`\mathcal{L}`$ theo $`x`$: nếu $`\mathcal{L}`$ lồi theo $`x`$ thì giải $`\nabla_x \mathcal{L} = 0`$ được $`x(\lambda, \mu)`$. Nếu $`\mathcal{L}`$ tuyến tính theo $`x`$ thì $`q`$ hữu hạn chỉ khi hệ số của $`x`$ bằng 0 (sinh ra ràng buộc của (D)).
3. Thế vào được $`q(\lambda, \mu)`$; ghi rõ miền $`\mathrm{dom}\, q`$.
4. Giải $`\max q`$ với $`\lambda \ge 0`$. Nếu $`(P)`$ lồi và thỏa Slater thì khôi phục $`x^* = x(\lambda^*, \mu^*)`$ và kiểm tra $`q(\lambda^*) = f(x^*)`$.

**Đối ngẫu LP và khôi phục nghiệm bằng điều kiện bù:** từ $`y^*`$ đã biết: (a) mọi $`y_i^* \ne 0`$ thì ràng buộc $`i`$ chặt; (b) mọi ràng buộc đối ngẫu lỏng ($`A_j^Ty^* \ne c_j`$) thì $`x_j = 0`$. Giải hệ tuyến tính còn lại để được $`x^*`$.

### Dạng câu hỏi có thể gặp

| Dạng câu | Cách làm nhanh | Lệnh Matlab |
|---|---|---|
| Tính $`q(\lambda)`$ tại một giá trị | Công thức $`x(\lambda)`$ rồi thế vào $`\mathcal{L}`$ | `res = dual_qp(Q, c, A, b); res.q(lambda, [])` |
| "Bài toán đối ngẫu của LP ... là" | Bảng quy tắc ở trên | `lp_dual(c, A, b)` (dạng chính tắc) |
| Điền $`y^*`$ hoặc $`d^*`$ | Giải (D) hoặc dùng $`d^* = p^*`$ | `r = lp_dual(...)`: `.y`, `.d`, `.compl` |
| Khôi phục $`x^*`$ từ $`y^*`$ | Điều kiện bù | `solve_lp` để đối chiếu |
| "Có khe hở đối ngẫu không?" | Lồi + Slater thì không; không lồi thì có thể có | `res.gap` |
| "$`(D)`$ có đạt nghiệm?" | Slater bảo đảm; thiếu Slater có thể không đạt | |

### Câu hỏi mẫu

**Câu 5.1.** Bài toán $`\min x_1^2 + x_2^2`$ với $`x_1 + x_2 \ge 2`$. Giá trị hàm đối ngẫu tại $`\lambda = 1`$ là:
A. $`1.5`$ B. $`2`$ C. $`1`$ D. $`-\infty`$

<details><summary>Đáp án</summary>

**A.** $`g = 2 - x_1 - x_2 \le 0`$, $`\mathcal{L} = x_1^2 + x_2^2 + \lambda(2 - x_1 - x_2)`$, cực tiểu tại $`x_1 = x_2 = \lambda/2`$, nên $`q(\lambda) = 2\lambda - \lambda^2/2`$ và $`q(1) = 1.5`$. Cực đại tại $`\lambda^* = 2`$ với $`d^* = 2 = p^*`$ ($`x^* = (1, 1)`$). `dual_qp`: `q(1) = 1.5, q(2) = 2, lambda* = 2, d* = p* = 2`.

</details>

**Câu 5.2.** Bài toán $`\min 2x_1 + 3x_2`$ với $`x_1 + x_2 \ge 4`$, $`x_1 + 3x_2 \ge 6`$, $`x \ge 0`$. Bài toán đối ngẫu và giá trị tối ưu là:
A. $`\max 4y_1 + 6y_2`$, $`y_1 + y_2 \le 2`$, $`y_1 + 3y_2 \le 3`$, $`y \ge 0`$; giá trị 9.
B. $`\max 4y_1 + 6y_2`$, $`y_1 + y_2 \ge 2`$, $`y_1 + 3y_2 \ge 3`$, $`y \ge 0`$; giá trị 12.
C. $`\min 4y_1 + 6y_2`$, $`y_1 + y_2 \le 2`$, $`y_1 + 3y_2 \le 3`$, $`y \le 0`$; giá trị 9.
D. $`\max 2y_1 + 3y_2`$, $`y_1 + y_2 \le 4`$, $`y_1 + 3y_2 \le 6`$, $`y \ge 0`$; giá trị 8.

<details><summary>Đáp án</summary>

**A.** Theo bảng: bài min có ràng buộc $`\ge`$ thì $`y_i \ge 0`$; biến $`x_j \ge 0`$ thì ràng buộc đối ngẫu $`\le c_j`$. Gốc: các đỉnh $`(0, 4), (3, 1), (6, 0)`$ cho $`f = 12, 9, 12`$, nên $`x^* = (3, 1)`$, $`p^* = 9`$. Đối ngẫu: $`x_1, x_2 \gt 0`$ nên cả hai ràng buộc đối ngẫu chặt, $`y_1 + y_2 = 2`$, $`y_1 + 3y_2 = 3`$, cho $`y^* = (1.5, 0.5)`$ và $`d^* = 6 + 3 = 9`$. `lp_dual` (thêm biến bù): `x* = [3 1 0 0], p* = 9; y* = [1.5 0.5], d* = 9`, mọi tích bù bằng 0.

</details>

**Câu 5.3.** Với bài toán lồi $`\min x`$ s.t. $`x^2 \le 0`$:
A. $`d^* \lt p^*`$. B. $`d^* = p^*`$ và $`(D)`$ đạt nghiệm. C. $`d^* = p^*`$ nhưng $`(D)`$ không đạt nghiệm. D. $`(D)`$ không có phương án.

<details><summary>Đáp án</summary>

**C.** $`p^* = 0`$ ($`D = \{0\}`$). Với $`\lambda \gt 0`$: $`x + \lambda x^2`$ cực tiểu tại $`x = -1/(2\lambda)`$, nên $`q(\lambda) = -1/(4\lambda)`$; còn $`q(0) = -\infty`$. Khi đó $`\sup_\lambda q = 0 = p^*`$ nhưng không đạt tại $`\lambda`$ hữu hạn nào ($`q(1) = -0.25`$, $`q(10) = -0.025`$, $`q(1000) = -0.00025`$). Nguyên nhân: Slater hỏng (không có $`x`$ với $`x^2 \lt 0`$), và tương ứng không có nhân tử KKT tại $`x^* = 0`$.

</details>

### Nếu ra tự luận

**Tự luận 5.1 — Lập bài toán đối ngẫu của một QP, giải và kiểm tra đối ngẫu mạnh.**
*Đề:* $`\min f = x_1^2 + x_2^2 - 2x_1 - 4x_2`$ với $`x_1 + x_2 \le 1`$.

*Lời giải.*
1. **Lồi, Slater.** $`\nabla^2 f = 2I \succ 0`$, ràng buộc affine nên bài toán lồi. $`\hat x = (0, 0)`$ cho $`g = -1 \lt 0`$, nên Slater thỏa. Do đó có đối ngẫu mạnh và (D) đạt nghiệm.
2. **Hàm Lagrange:** $`\mathcal{L} = x_1^2 + x_2^2 - 2x_1 - 4x_2 + \lambda(x_1 + x_2 - 1)`$, $`\lambda \ge 0`$.
3. **Hàm đối ngẫu.** $`\nabla_x \mathcal{L} = 0`$ cho $`x(\lambda) = (1 - \lambda/2,\ 2 - \lambda/2)`$. Thế vào:

```math
q(\lambda) = \Big(-1 + \frac{\lambda^2}{4}\Big) + \Big(-4 + \frac{\lambda^2}{4}\Big) + \lambda(2 - \lambda) = -5 + 2\lambda - \frac{\lambda^2}{2}.
```

4. **Bài toán đối ngẫu:** $`\max_{\lambda \ge 0} -5 + 2\lambda - \lambda^2/2`$. Có $`q'(\lambda) = 2 - \lambda = 0`$ nên $`\lambda^* = 2 \ge 0`$ và $`d^* = -5 + 4 - 2 = -3`$.
5. **Khôi phục nghiệm gốc:** $`x^* = x(2) = (0, 1)`$, thỏa $`0 + 1 \le 1`$ (chặt, phù hợp $`\lambda^* \gt 0`$). $`f(x^*) = 1 - 4 = -3 = d^*`$: đối ngẫu mạnh.
6. **KKT:** $`2x_1 - 2 + \lambda = 0`$ và $`2x_2 - 4 + \lambda = 0`$ tại $`(0, 1)`$, $`\lambda = 2`$; bù $`\lambda g = 0`$. Thỏa.
*Kiểm chứng:* `dual_qp` cho `lambda* = 2, d* = p* = -3, x(lambda*) = [0 1]`; `q(0) = -5`, `q(1) = -3.5` khớp công thức.

**Tự luận 5.2 — Đối ngẫu LP, khôi phục nghiệm gốc bằng điều kiện bù.**
*Đề:* $`\max z = 3x_1 + 2x_2`$ với $`x_1 + x_2 \le 4`$, $`x_1 + 3x_2 \le 6`$, $`x_1 \le 3.5`$, $`x \ge 0`$. Viết bài toán đối ngẫu, chứng tỏ $`y = (2, 0, 1)`$ là nghiệm đối ngẫu, từ đó tìm nghiệm gốc.

*Lời giải.*
1. **Đối ngẫu** (bài max, ràng buộc $`\le`$ cho $`y_i \ge 0`$; biến $`x_j \ge 0`$ cho ràng buộc đối ngẫu $`\ge c_j`$):

```math
(D):\quad \min w = 4y_1 + 6y_2 + 3.5y_3\quad \text{s.t.}\quad y_1 + y_2 + y_3 \ge 3,\quad y_1 + 3y_2 \ge 2,\quad y \ge 0 .
```

2. **$`y = (2, 0, 1)`$ chấp nhận được:** $`2 + 0 + 1 = 3 \ge 3`$, $`2 + 0 = 2 \ge 2`$, $`y \ge 0`$. Giá trị $`w = 8 + 3.5 = 11.5`$.
3. **Điều kiện bù:** $`y_1 = 2 \gt 0`$ nên $`x_1 + x_2 = 4`$; $`y_3 = 1 \gt 0`$ nên $`x_1 = 3.5`$. Suy ra $`x = (3.5, 0.5)`$. Kiểm tra ràng buộc 2: $`3.5 + 1.5 = 5 \le 6`$, lỏng, phù hợp với $`y_2 = 0`$. Cả hai ràng buộc đối ngẫu đều chặt, phù hợp với $`x_1, x_2 \gt 0`$.
4. $`z = 10.5 + 1 = 11.5 = w`$, nên theo hệ quả đối ngẫu yếu cả hai đều tối ưu: $`x^* = (3.5, 0.5)`$, $`y^* = (2, 0, 1)`$, giá trị tối ưu $`11.5`$.
5. **Giá bóng:** tăng vế phải ràng buộc 1 thêm $`\varepsilon`$ thì $`z^*`$ tăng $`2\varepsilon`$; ràng buộc 2 đang lỏng nên giá bóng bằng 0.
*Kiểm chứng:* `lp_dual` (đưa về $`\min -z`$ với biến bù) cho `x* = [3.5 0.5], z* = 11.5, y = [2 0 1]`; giải trực tiếp (D) bằng `solve_lp` cho $`y = (2, 0, 1)`$, $`w = 11.5`$.

### Bẫy

- Đối ngẫu yếu **luôn** đúng; đối ngẫu mạnh cần lồi + Slater (LP thì chỉ cần có nghiệm). Bài không lồi có thể có khe hở dương (câu 5, ví dụ 3: $`\min -x^2`$ trên $`[0, 2]`$ với $`x \le 1`$, khe hở bằng 1).
- Khi $`\mathcal{L}`$ tuyến tính theo $`x`$, quên điều kiện "hệ số bằng 0" thì sẽ viết sai miền của $`q`$.
- Bài **max** phải đổi chiều: đối ngẫu là min, và dấu $`y`$ theo bảng đọc từ phải sang trái.
- Nghiệm đối ngẫu từ solver dạng $`\min`$ của một bài max có dấu ngược ($`y_{\max} = -y_{\min}`$).
- $`(D)`$ có thể không đạt nghiệm dù không có khe hở (câu 5.3).

---

## Câu 13. Quy hoạch lồi: chính quy và KKT

### Lý thuyết cốt lõi

**Bài toán lồi:** $`f`$ lồi, $`g_i`$ lồi, $`h_j`$ affine ($`h_j = a_j^Tx - b_j`$). Khi đó miền chấp nhận được lồi, cực tiểu địa phương là toàn cục, tập nghiệm lồi; nếu $`f`$ lồi chặt thì nghiệm duy nhất. Với hàm khả vi: $`f`$ lồi khi và chỉ khi $`f(y) \ge f(x) + \nabla f(x)^T(y - x)`$ (bổ đề 1).

**Định lý (KKT là điều kiện đủ, không cần chính quy).** Với bài toán lồi khả vi, mọi điểm KKT là nghiệm tối ưu toàn cục.
*Chứng minh.* Lấy $`x \in D`$. Từ bổ đề 1 và điều kiện dừng:

```math
f(x) - f(x^*) \ge \nabla f(x^*)^T(x - x^*) = -\sum_i \lambda_i \nabla g_i(x^*)^T(x - x^*) - \sum_j \mu_j a_j^T(x - x^*).
```

Với $`h_j`$ affine: $`a_j^T(x - x^*) = h_j(x) - h_j(x^*) = 0`$. Với $`\lambda_i \gt 0`$: theo điều kiện bù $`g_i(x^*) = 0`$, nên $`\nabla g_i(x^*)^T(x - x^*) \le g_i(x) - g_i(x^*) = g_i(x) \le 0`$. Vậy vế phải $`\ge 0`$, tức $`f(x) \ge f(x^*)`$. $`\blacksquare`$

**Định lý (KKT là điều kiện cần, Kuhn–Tucker).** Bài toán lồi thỏa Slater, $`x^*`$ tối ưu, thì tồn tại $`\lambda^* \ge 0, \mu^*`$ với $`\mathcal{L}(x, \lambda^*, \mu^*) \ge \mathcal{L}(x^*, \lambda^*, \mu^*) = f^*`$ cho mọi $`x`$ và $`\lambda_i^* g_i(x^*) = 0`$. Chứng minh: tách hai tập lồi $`\mathcal{A} = \{(u, t) : \exists x, g(x) \le u, f(x) \le t\}`$ và $`\mathcal{B} = \{u \lt 0, t \lt f^*\}`$; Slater cho $`\lambda_0 \gt 0`$.
**Hệ quả:** với bài lồi khả vi thỏa Slater, ba điều sau tương đương: $`x^*`$ tối ưu, $`x^*`$ là điểm KKT, $`(x^*, \lambda^*, \mu^*)`$ là điểm yên ngựa.

**Slater và các dạng:**
- *Chặt:* có $`\hat x`$ với mọi $`g_i(\hat x) \lt 0`$ và $`h(\hat x) = 0`$. Là tính chất của cả bài toán, không cần biết $`x^*`$.
- *Tinh:* ràng buộc affine chỉ cần $`g_i(\hat x) \le 0`$. Nếu mọi ràng buộc affine (LP, QP) thì không cần điều kiện gì.
- Không có ràng buộc đẳng thức: Slater tương đương MFCQ tại một điểm, tương đương MFCQ tại mọi điểm.
- **Chặn nhân tử:** $`\sum_i \lambda_i^* \le \dfrac{f(\hat x) - f^*}{\min_i(-g_i(\hat x))}`$.
- Thiếu Slater: nghiệm có thể không là điểm KKT, hoặc tập nhân tử không bị chặn.

**Không khả vi:** KKT dạng dưới vi phân $`0 \in \partial f(x^*) + \sum \lambda_i \partial g_i(x^*) + \sum \mu_j a_j`$ (ví dụ $`\min \Vert x \Vert_1`$, $`a^Tx = 1`$ cho nghiệm thưa).

### Thuật toán

**Quy trình tự luận "bài toán lồi" (6 bước):**
1. **Lồi:** chỉ ra $`f`$ lồi (Hessian $`\succeq 0`$ hoặc $`\varphi'' \ge 0`$ cho hàm tách biến), $`g_i`$ lồi, $`h_j`$ affine.
2. **Slater:** chỉ ra một $`\hat x`$ cụ thể (thường là tâm của miền, như $`(1/n, \dots, 1/n)`$). Kết luận: KKT là điều kiện cần và đủ, đối ngẫu mạnh.
3. **KKT:** viết đủ 4 nhóm.
4. **Giải:** khử nhân tử theo trường hợp $`x_i \gt 0`$ / $`x_i = 0`$. Thường được công thức dạng $`x_i = \max(0, \cdot)`$ ("water-filling").
5. **Đối ngẫu:** cực tiểu $`\mathcal{L}`$ theo từng biến, được $`q`$; kiểm tra $`q(\lambda^*, \mu^*) = f^*`$.
6. **Kiểm chứng** bằng `solve_nlp` / `water_filling`.

**Kiểm tra Slater bằng máy:** bài toán phụ $`\min s`$ với $`g_i(x) \le s`$ (phi tuyến), $`g_i(x) \le 0`$ (affine, dạng tinh), $`h = 0`$, $`s \ge -1`$. Slater thỏa khi và chỉ khi $`s^* \lt 0`$ (`slater_check`).

### Dạng câu hỏi có thể gặp

| Dạng câu | Cách làm nhanh | Lệnh Matlab |
|---|---|---|
| "Bài toán có lồi không?" | Hessian của $`f`$, $`g_i`$; $`h_j`$ phải affine | `eig(H)` |
| "Slater có thỏa không?" | Tìm $`\hat x`$; nếu miền là một điểm hoặc có $`g_1 = -g_2`$ phi tuyến thì hỏng | `[ok, xhat, s] = slater_check(prob, x0, isAffine)` |
| "Điểm KKT của bài lồi là" | Cực tiểu toàn cục, không cần chính quy | `kkt_check` |
| Water-filling: điền $`x^*`$, $`\nu`$, $`\lambda`$ | Sắp xếp $`\alpha`$, thử số kênh $`k`$: $`L = (P + \sum_{i \le k} \alpha_i)/k`$ | `[x, nu, lam, L] = water_filling(alpha, P)` |
| Chặn trên tổng nhân tử | Công thức chặn với điểm Slater | |

### Câu hỏi mẫu

**Câu 13.1.** Bài toán $`\min x_2`$ với $`g_1 = x_1^2 + x_2^2 - 1 \le 0`$, $`g_2 = 1 - x_1 \le 0`$. Khẳng định đúng:
A. Slater dạng tinh thỏa vì $`g_2`$ affine, nên $`x^* = (1, 0)`$ là điểm KKT.
B. Slater không thỏa (cả dạng tinh); $`x^* = (1, 0)`$ là nghiệm nhưng không là điểm KKT.
C. Bài toán không lồi.
D. Bài toán vô nghiệm.

<details><summary>Đáp án</summary>

**B.** $`x_1 \ge 1`$ và $`x_1^2 + x_2^2 \le 1`$ cho $`D = \{(1, 0)\}`$, nên $`x^* = (1, 0)`$. Dạng tinh cần $`g_1(\hat x) \lt 0`$ (phi tuyến) và $`x_1 \ge 1`$: không có điểm nào. Tại $`(1, 0)`$: $`\nabla g_1 = (2, 0)`$, $`\nabla g_2 = (-1, 0)`$ cùng phương, $`\nabla f = (0, 1)`$; thành phần thứ hai của hệ dừng là $`1 = 0`$, vô nghiệm. `slater_check`: dạng chặt và dạng tinh đều `ok = 0` ($`s^* \approx 0`$); `kkt_check`: `KKT = 0, LICQ = 0`, $`\Vert\nabla_x\mathcal{L}\Vert = 1`$.

</details>

**Câu 13.2.** Water-filling $`\min -\sum_{i=1}^3 \log(\alpha_i + x_i)`$ với $`\sum x_i = 3`$, $`x \ge 0`$, $`\alpha = (1, 2, 4)`$. Điền $`x^*`$, $`\nu`$, $`\lambda_3`$.

<details><summary>Đáp án</summary>

KKT cho $`x_i = \max(0, L - \alpha_i)`$ với $`L = 1/\nu`$. Thử 3 kênh: $`L = (3 + 7)/3 = 3.33 \lt \alpha_3 = 4`$, không hợp lệ. Thử 2 kênh: $`L = (3 + 1 + 2)/2 = 3 \in [2, 4]`$, hợp lệ. Vậy $`x^* = (2, 1, 0)`$, $`\nu = 1/3`$, và $`\lambda_3 = \nu - 1/\alpha_3 = 1/3 - 1/4 = 1/12`$. `water_filling` cho `x = [2 1 0], nu = 0.333333, lambda = [0 0 0.08333]`.

</details>

**Câu 13.3.** Bài toán lồi $`\min x_1^2 + x_2^2`$ với $`x_1 + x_2 \ge 1`$; điểm $`\bar x = (0.5, 0.5)`$ thỏa KKT với $`\lambda = 1`$. Kết luận:
A. $`\bar x`$ là cực tiểu địa phương, chưa chắc toàn cục. B. $`\bar x`$ là cực tiểu toàn cục. C. Cần kiểm tra Slater mới kết luận được. D. Cần LICQ mới kết luận được.

<details><summary>Đáp án</summary>

**B.** Với bài toán lồi, KKT là điều kiện **đủ** cho tối ưu toàn cục mà không cần điều kiện chính quy. Slater/LICQ chỉ cần cho chiều ngược lại (KKT là điều kiện cần). `kkt_check`: `KKT = 1, lambda = 1`.

</details>

### Nếu ra tự luận

Dạng "hàm tách biến trên đơn hình" với $`f = \sum x_i \ln(a + x_i)`$ đã có lời giải đầy đủ ở [Đề số 1, Phần B](../../test1-paper1/README.md#phần-b). Dưới đây là hai biến thể khác cùng khuôn, và một bài chứng minh lý thuyết.

**Tự luận 13.1 — Water-filling: lồi, Slater, KKT, đối ngẫu, giải.**
*Đề:* $`\min f(x) = -\sum_{i=1}^3 \log(\alpha_i + x_i)`$ với $`x_1 + x_2 + x_3 = 1`$, $`x \ge 0`$, $`\alpha = (1, 1.5, 3)`$.

*Lời giải.*
1. **Lồi.** Đặt $`\varphi_i(t) = -\log(\alpha_i + t)`$ trên $`t \gt -\alpha_i`$, có $`\varphi_i'' = 1/(\alpha_i + t)^2 \gt 0`$. Vậy $`f`$ lồi chặt, các ràng buộc affine, (P) là bài toán lồi và có nghiệm duy nhất (miền compact).
2. **Chính quy.** Mọi ràng buộc affine nên KKT là điều kiện cần (Slater dạng tinh). Điểm $`(1/3, 1/3, 1/3) \gt 0`$ còn thỏa Slater chặt. Do đó KKT là điều kiện cần và đủ, và có đối ngẫu mạnh.
3. **Lagrange** (nhân tử $`\lambda_i \ge 0`$ cho $`-x_i \le 0`$, $`\nu`$ cho đẳng thức): $`\mathcal{L} = -\sum_i \log(\alpha_i + x_i) - \sum_i \lambda_i x_i + \nu(\sum_i x_i - 1)`$.
4. **KKT:**

```math
-\frac{1}{\alpha_i + x_i} - \lambda_i + \nu = 0,\qquad \lambda_i \ge 0,\qquad \lambda_i x_i = 0,\qquad x_i \ge 0,\qquad \textstyle\sum_i x_i = 1 .
```

5. **Giải.** Nếu $`x_i \gt 0`$ thì $`\lambda_i = 0`$ và $`x_i = 1/\nu - \alpha_i`$. Nếu $`x_i = 0`$ thì $`\lambda_i = \nu - 1/\alpha_i \ge 0`$, tức $`\alpha_i \ge 1/\nu`$. Gộp lại, với $`L = 1/\nu`$: $`x_i = \max(0, L - \alpha_i)`$ và $`\sum_i \max(0, L - \alpha_i) = 1`$. Hàm vế trái tăng chặt theo $`L`$ nên $`L`$ duy nhất. Thử 3 kênh: $`L = (1 + 5.5)/3 = 2.17 \lt 3`$, không hợp lệ. Thử 2 kênh: $`L = (1 + 1 + 1.5)/2 = 1.75 \in [1.5, 3]`$, hợp lệ. Vậy

```math
x^* = (0.75,\ 0.25,\ 0),\qquad \nu^* = \frac47,\qquad \lambda^* = \Big(0,\ 0,\ \frac47 - \frac13 = \frac{5}{21}\Big),\qquad f^* = -2\ln 1.75 - \ln 3 \approx -2.217844 .
```

6. **Đối ngẫu.** Với $`\mu_i = \nu - \lambda_i`$: $`\inf_{t \gt -\alpha_i} [-\log(\alpha_i + t) + \mu_i t]`$ hữu hạn chỉ khi $`\mu_i \gt 0`$, đạt tại $`t = 1/\mu_i - \alpha_i`$, giá trị $`\log\mu_i + 1 - \alpha_i\mu_i`$. Do đó

```math
q(\nu, \lambda) = -\nu + \sum_{i=1}^3 \big(\log(\nu - \lambda_i) + 1 - \alpha_i(\nu - \lambda_i)\big),\qquad (D):\ \max q\ \ \text{s.t.}\ \lambda \ge 0,\ \nu - \lambda_i \gt 0 .
```

Tại $`(\nu^*, \lambda^*)`$: $`q = -2.217844 = f^*`$, đúng đối ngẫu mạnh.
7. **Ý nghĩa:** kênh 3 có "đáy" $`\alpha_3 = 3`$ cao hơn mức nước $`L = 1.75`$ nên không được cấp; $`\lambda_3 = 5/21`$ là mức tăng của $`f^*`$ trên mỗi đơn vị nếu buộc phải cấp cho kênh 3.
*Kiểm chứng:* `water_filling([1;1.5;3], 1)` cho `x = [0.75 0.25 0], nu = 0.571429, lambda3 = 0.238095`; $`q(\nu^*, \lambda^*) = f^* = -2.217844`$; `solve_nlp` cho cùng $`x`$ và $`f`$.

**Tự luận 13.2 — Hàm bậc hai tách biến trên đơn hình (xét trường hợp).**
*Đề:* $`\min f = \sum_{i=1}^3 \big(\tfrac12 x_i^2 - c_i x_i\big)`$ với $`\sum x_i = 1`$, $`x \ge 0`$, $`c = (2, 1.5, 0)`$.

*Lời giải.*
1. $`\nabla^2 f = I \succ 0`$ nên $`f`$ lồi chặt; ràng buộc affine; Slater: $`(1/3, 1/3, 1/3)`$. KKT là điều kiện cần và đủ, nghiệm duy nhất.
2. $`\mathcal{L} = f + \nu(1 - \sum x_i) - \sum \lambda_i x_i`$. KKT: $`x_i - c_i - \nu - \lambda_i = 0`$, $`\lambda_i \ge 0`$, $`\lambda_i x_i = 0`$, $`\sum x_i = 1`$, $`x \ge 0`$. Từ đó $`x_i = \max(0, c_i + \nu)`$ (phép chiếu $`c`$ lên đơn hình, dịch đi $`\nu`$).
3. Tìm $`\nu`$: giả sử chỉ $`x_1, x_2 \gt 0`$: $`(2 + \nu) + (1.5 + \nu) = 1`$ cho $`\nu = -1.25`$. Kiểm tra $`c_3 + \nu = -1.25 \lt 0`$, nên $`x_3 = 0`$ và $`\lambda_3 = 0 - 0 + 1.25 = 1.25 \ge 0`$. (Giả sử cả ba dương: $`3.5 + 3\nu = 1`$, $`\nu = -5/6`$, $`x_3 = -5/6 \lt 0`$, loại.)
4. **Kết quả:** $`x^* = (0.75, 0.25, 0)`$, $`\nu^* = -1.25`$, $`\lambda^* = (0, 0, 1.25)`$, $`f^* = \tfrac12(0.5625 + 0.0625) - (1.5 + 0.375) = -1.5625`$.
5. **Đối ngẫu:** $`\inf_t (\tfrac12 t^2 - (c_i + \mu_i) t) = -\tfrac12 (c_i + \mu_i)^2`$ với $`\mu_i = \nu + \lambda_i`$, nên $`q(\nu, \lambda) = \nu - \tfrac12 \sum_i (c_i + \nu + \lambda_i)^2`$. Tại $`(\nu^*, \lambda^*)`$: $`-1.25 - \tfrac12(0.5625 + 0.0625 + 0) = -1.5625 = f^*`$.
*Kiểm chứng:* `solve_nlp` cho $`x = (0.75, 0.25, 0)`$, $`f = -1.5625`$.

**Tự luận 13.3 — Phát biểu và chứng minh: với bài toán lồi, KKT là điều kiện đủ; khi nào KKT là điều kiện cần.**
*Lời giải mẫu.*
1. Định nghĩa bài toán lồi, bổ đề 1 (bất đẳng thức gradient của hàm lồi), và chứng minh định lý "KKT là điều kiện đủ" như ở mục lý thuyết (chú ý không dùng điều kiện chính quy).
2. Phát biểu Slater (chặt và tinh), định lý Kuhn–Tucker (KKT là điều kiện cần khi có Slater) và 5 bước chứng minh qua định lý tách. Nhấn mạnh bước dùng Slater: loại $`\lambda_0 = 0`$.
3. **Ví dụ Slater hỏng** (câu 13.1 ở trên, hoặc câu 4 ví dụ 6: $`\min x`$, $`x^2 \le 0`$): nghiệm tồn tại nhưng không phải điểm KKT. **Ví dụ Slater hỏng nhưng KKT vẫn có:** hai hình tròn tiếp xúc $`(x_1 \mp 1)^2 + x_2^2 \le 1`$ với $`f = x_1`$: điểm KKT với tập nhân tử $`\lambda_1 = \lambda_2 + 1/2`$, không bị chặn. Từ đó: Slater là điều kiện đủ để KKT cần, không phải điều kiện cần.
4. Kết luận: với bài toán lồi khả vi thỏa Slater, $`x^*`$ tối ưu ⇔ KKT ⇔ điểm yên ngựa ⇔ (cùng nghiệm đối ngẫu) đối ngẫu mạnh.

### Bẫy

- "Mọi ràng buộc đều affine nên Slater tinh luôn thỏa" chỉ đúng khi **không có** ràng buộc phi tuyến. Một ràng buộc phi tuyến cũng phải có điểm làm nó âm thực sự (câu 13.1).
- Ràng buộc đẳng thức phi tuyến ($`x_1^2 + x_2^2 = 1`$) làm bài toán **không lồi**, dù hàm mục tiêu lồi.
- KKT đủ cho bài lồi **không cần** Slater; Slater chỉ cần cho chiều "nghiệm ⇒ KKT" và cho đối ngẫu mạnh.
- Water-filling: phải kiểm tra $`L`$ nằm giữa $`\alpha_{(k)}`$ và $`\alpha_{(k+1)}`$ (sắp xếp tăng); quên kiểm tra sẽ ra $`x_i \lt 0`$.
- Khi viết $`q`$ cho hàm có miền (log), ghi rõ điều kiện $`\mu_i \gt 0`$; ngoài miền đó $`q = -\infty`$.
