# Tìm kiếm một chiều (câu 6, 7)

[← Tài liệu ôn thi](../README.md)

Số liệu trong các câu hỏi mẫu lấy từ output của [`check_unconstrained.m`](../check_unconstrained.m).

## Câu 6. Giải thuật chia đôi

### Lý thuyết cốt lõi

- **Hàm đơn mốt** trên $`[a, b]`$: có $`x^*`$ để $`f`$ giảm thực sự trên $`[a, x^*]`$ và tăng thực sự trên $`[x^*, b]`$. Hàm lồi (chặt) là đơn mốt; điều ngược lại không đúng.
- **Bổ đề loại bỏ:** với $`f`$ đơn mốt và $`x_1 \lt x_2`$: nếu $`f(x_1) \lt f(x_2)`$ thì $`x^* \in [a, x_2]`$; nếu $`f(x_1) \gt f(x_2)`$ thì $`x^* \in [x_1, b]`$; nếu bằng nhau thì $`x^* \in [x_1, x_2]`$.
- **Chia đôi dùng giá trị hàm (dichotomous):** $`c = (a+b)/2`$, $`x_{1,2} = c \mp \delta`$, với $`0 \lt \delta \lt \varepsilon/2`$. Độ dài khoảng

```math
L_{k+1} = \frac{L_k}{2} + \delta,\qquad L_k = \frac{L_0 - 2\delta}{2^k} + 2\delta,\qquad
k \ge \log_2 \frac{L_0 - 2\delta}{\varepsilon - 2\delta}\ \ \text{để } L_k \le \varepsilon .
```

  Mỗi bước tốn 2 lần tính $`f`$. Vì $`L_k \to 2\delta`$ chứ không về 0, nên bắt buộc $`\delta \lt \varepsilon/2`$.
- **Chia đôi dùng đạo hàm:** $`c = (a+b)/2`$; nếu $`f'(c) \gt 0`$ thì giữ $`[a, c]`$, nếu $`f'(c) \lt 0`$ thì giữ $`[c, b]`$, nếu $`f'(c) = 0`$ thì dừng. Ta có $`L_k = L_0/2^k`$, cần $`k \ge \log_2(L_0/\varepsilon)`$ bước, mỗi bước 1 lần tính $`f'`$. Đây chính là phương pháp chia đôi giải phương trình $`f'(x) = 0`$.
- **Sai số:** lấy $`\bar x`$ là trung điểm khoảng cuối thì $`|\bar x - x^*| \le L_k/2`$. Tốc độ hội tụ **tuyến tính** (hệ số 1/2).

### Thuật toán

1. Cho $`f`$, $`[a, b]`$, $`\varepsilon`$, $`\delta \lt \varepsilon/2`$.
2. Lặp khi $`b - a \gt \varepsilon`$: $`c = (a+b)/2`$, $`x_1 = c - \delta`$, $`x_2 = c + \delta`$. Nếu $`f(x_1) \lt f(x_2)`$ thì $`b \leftarrow x_2`$, ngược lại $`a \leftarrow x_1`$.
3. Trả về $`\bar x = (a + b)/2`$.

### Dạng câu hỏi có thể gặp

| Dạng | Cách làm nhanh | Matlab |
|---|---|---|
| Điền $`x^*`$ trên $`[a, b]`$ | Vẽ đồ thị; nếu $`f`$ đơn mốt thì đáp án là nghiệm của $`f'(x) = 0`$ | `fminbnd(f,a,b)`, `fzero(df,[a b])`, `bisection_min(f,a,b,struct('tol',1e-4))` |
| Số bước cần thiết | Công thức $`\log_2`$ ở trên, làm tròn **lên** | tính tay |
| Khoảng sau 1–2 bước | So sánh $`f(c \mp \delta)`$ hoặc xét dấu $`f'(c)`$ | `bisection_min(..., struct('display', true))` |
| So sánh với lát cắt vàng | Theo số lần tính $`f`$: chia đôi có hệ số $`(1/2)^{1/2} \approx 0.707`$ mỗi lần tính, lát cắt vàng $`0.618`$ | — |

### Câu hỏi mẫu

**6.1.** Dùng chia đôi (giá trị hàm) với $`\delta = 10^{-3}`$ trên $`[0, 2]`$. Số bước tối thiểu để độ dài khoảng $`\le 0.01`$ là:
A. 7 B. 8 C. 9 D. 10

<details><summary>Đáp án</summary>

**B.** $`k \ge \log_2 \frac{2 - 0.002}{0.01 - 0.002} = \log_2 249.75 = 7.96`$, nên $`k = 8`$.

</details>

**6.2.** Dùng chia đôi theo đạo hàm cho $`f(x) = x^2 - \sin x`$ trên $`[0, 1]`$. Khoảng chứa cực tiểu sau 2 bước là:
A. $`[0, 0.5]`$ B. $`[0.5, 0.75]`$ C. $`[0.25, 0.5]`$ D. $`[0.5, 1]`$

<details><summary>Đáp án</summary>

**C.** $`f'(x) = 2x - \cos x`$. Bước 1: $`f'(0.5) = 0.1224 \gt 0`$, giữ $`[0, 0.5]`$. Bước 2: $`f'(0.25) = -0.4689 \lt 0`$, giữ $`[0.25, 0.5]`$. (Nghiệm đúng $`x^* = 0.450184`$.)

</details>

**6.3.** $`f(x) = x^4 - 14x^3 + 60x^2 - 70x`$ có cực tiểu trên $`[0, 2]`$ (sai số $`10^{-5}`$) là $`x^* = \ ?`$ Ở bước 1 với $`\delta = 10^{-3}`$, giữ lại khoảng nào?

<details><summary>Đáp án</summary>

$`x^* \approx 0.78088`$ (`bisection_min`: 19 bước, 39 lần tính $`f`$; [câu 6, ví dụ 1](../../../solutions/q06-bisection-search/)). Bước 1: $`f(0.999) = -23.0120 \lt f(1.001) = -22.9880`$, nên giữ $`[0, 1.001]`$.

</details>

### Nếu ra tự luận

"Trình bày giải thuật chia đôi, chứng minh công thức độ dài khoảng, áp dụng cho ví dụ." Dàn ý:
1. Định nghĩa hàm đơn mốt và chứng minh bổ đề loại bỏ (nếu $`f(x_1) \lt f(x_2)`$ mà $`x^* \gt x_2`$ thì $`x_1 \lt x_2 \lt x^*`$ thuộc đoạn $`f`$ giảm, nên $`f(x_1) \gt f(x_2)`$, mâu thuẫn).
2. Các bước thuật toán (cả hai biến thể).
3. Định lý $`L_k = (L_0 - 2\delta)/2^k + 2\delta`$ chứng minh bằng quy nạp; suy ra số bước và điều kiện $`\delta \lt \varepsilon/2`$.
4. Ví dụ: lập bảng 3–4 bước ($`a_k, b_k, x_1, x_2, f(x_1), f(x_2)`$).

Báo cáo đầy đủ: [câu 6](../../../solutions/q06-bisection-search/).

### Bẫy

- Quên điều kiện $`\delta \lt \varepsilon/2`$, hoặc làm tròn số bước xuống.
- Áp dụng cho hàm không đơn mốt: thuật toán vẫn chạy nhưng có thể loại mất cực tiểu ([câu 6, ví dụ 3](../../../solutions/q06-bisection-search/)).
- Đếm nhầm "số bước" và "số lần tính $`f`$" (gấp đôi với biến thể dùng giá trị hàm).

---

## Câu 7. Giải thuật lát cắt vàng

### Lý thuyết cốt lõi

- Hai điểm thử chia đoạn theo tỉ lệ vàng:

```math
x_1 = b - \tau (b - a),\qquad x_2 = a + \tau (b - a),\qquad \tau = \frac{\sqrt5 - 1}{2} \approx 0.618,\qquad \tau^2 = 1 - \tau .
```

- **Tái dùng điểm:** nếu giữ $`[a, x_2]`$ thì $`x_1`$ cũ trở thành $`x_2`$ mới; nếu giữ $`[x_1, b]`$ thì $`x_2`$ cũ trở thành $`x_1`$ mới. Vì vậy mỗi bước chỉ tính $`f`$ tại **một** điểm mới (trừ bước đầu cần 2).
- $`L_k = \tau^k L_0`$, cần $`k \ge \dfrac{\ln(\varepsilon / L_0)}{\ln \tau}`$ bước. Tổng số lần tính $`f`$ là $`k + 1`$.
- **Vì sao chọn $`\tau`$:** yêu cầu (i) hai điểm đối xứng và (ii) tỉ lệ thu hẹp không đổi giữa các bước, đồng thời tái dùng được một điểm, cho ra $`\tau^2 + \tau - 1 = 0`$.
- So với Fibonacci: Fibonacci tối ưu khi biết trước số lần tính, còn lát cắt vàng là giới hạn của nó khi số bước lớn.
- Chỉ dùng giá trị $`f`$, nên độ chính xác bị chặn ở khoảng $`\sqrt{\epsilon_{\text{máy}}} \approx 10^{-8}`$ (tương đối): gần cực tiểu $`f`$ rất phẳng, hiệu $`f(x_1) - f(x_2)`$ chìm trong sai số làm tròn.

### Thuật toán

1. $`x_1 = b - \tau(b-a)`$, $`x_2 = a + \tau(b-a)`$; tính $`f_1, f_2`$.
2. Lặp khi $`b - a \gt \varepsilon`$:
   - Nếu $`f_1 \lt f_2`$: $`b \leftarrow x_2`$, $`x_2 \leftarrow x_1`$, $`f_2 \leftarrow f_1`$, $`x_1 = b - \tau(b - a)`$, tính $`f_1`$.
   - Ngược lại: $`a \leftarrow x_1`$, $`x_1 \leftarrow x_2`$, $`f_1 \leftarrow f_2`$, $`x_2 = a + \tau(b - a)`$, tính $`f_2`$.
3. Trả về $`(a + b)/2`$.

### Dạng câu hỏi có thể gặp

| Dạng | Cách làm nhanh | Matlab |
|---|---|---|
| Điền $`x^*`$ | Vẽ đồ thị, tìm điểm đổi dấu của $`f'`$, thu hẹp về đoạn đơn mốt | `golden(f,a,b,1e-4)` (toolkit), `golden_section`, `fminbnd` |
| Số bước | $`\lceil \ln(\varepsilon/L_0)/\ln 0.618 \rceil`$ | tính tay |
| Hai điểm thử đầu tiên / điểm mới ở bước 2 | Công thức $`x_{1,2}`$; điểm mới $`= b' - \tau(b' - a)`$ | `golden_section(..., struct('display', true))` |
| Cực tiểu địa phương "trong đoạn" khi $`f`$ không đơn mốt | Xét dấu $`f'(a)`$, $`f'(b)`$ và các điểm dừng (bẫy câu 4 Đề số 1) | `fzero(df, [..])` |

### Câu hỏi mẫu

**7.1.** Số bước lát cắt vàng tối thiểu để thu hẹp $`[0, 3]`$ về độ dài $`\le 10^{-3}`$ là:
A. 15 B. 16 C. 17 D. 18

<details><summary>Đáp án</summary>

**C.** $`k \ge \ln(10^{-3}/3)/\ln 0.618 = 16.64`$, nên $`k = 17`$.

</details>

**7.2.** Áp dụng lát cắt vàng cho $`f(x) = x^2 - \sin x`$ trên $`[0, 1]`$. Hai điểm thử đầu tiên, khoảng giữ lại sau bước 1, và điểm cần tính mới ở bước 2 là gì?

<details><summary>Đáp án</summary>

$`x_1 = 0.381966`$, $`f(x_1) = -0.226847`$; $`x_2 = 0.618034`$, $`f(x_2) = -0.197468`$. Vì $`f(x_1) \lt f(x_2)`$ nên giữ $`[0, 0.618034]`$. Điểm $`x_1`$ cũ thành $`x_2`$ mới; điểm mới là $`x_1 = 0.618034 - 0.618034 \cdot 0.618034 = 0.236068`$.

</details>

**7.3.** Dùng lát cắt vàng với sai số $`10^{-4}`$, cực tiểu của $`f(x) = e^x - 2x`$ trên $`[0, 2]`$ là $`x^* = \ ?`$

<details><summary>Đáp án</summary>

$`x^* \approx 0.6931`$. `golden_section` cho $`0.693145`$ sau 21 bước. Nghiệm đúng: $`f'(x) = e^x - 2 = 0`$ nên $`x^* = \ln 2 = 0.693147`$.

</details>

### Nếu ra tự luận

"Trình bày giải thuật lát cắt vàng, giải thích tỉ lệ $`\tau`$, so sánh với chia đôi, áp dụng cho ví dụ." Dàn ý:
1. Bổ đề loại bỏ (như câu 6).
2. Suy ra $`\tau`$ từ yêu cầu đối xứng và tái dùng điểm: $`\tau^2 = 1 - \tau`$.
3. Các bước thuật toán và $`L_k = \tau^k L_0`$.
4. So sánh theo số lần tính $`f`$: lát cắt vàng $`0.618`$ mỗi lần, chia đôi $`0.707`$ mỗi lần.
5. Ví dụ: bảng 3–4 bước ($`a_k, b_k, x_1, x_2, f_1, f_2`$), chỉ ra điểm được tái dùng.

Báo cáo đầy đủ: [câu 7](../../../solutions/q07-golden-section-search/). Bảng mẫu 23 bước của câu 4 Đề số 1: [lời giải](../../test1-paper1/).

### Bẫy

- Hàm không đơn mốt: câu 4 Đề số 1 ($`f`$ có cực đại tại 3.836 và cực tiểu ở biên $`x = 5`$). Lát cắt vàng ra đúng chỉ vì ngay bước đầu nó loại phần chứa điểm cực đại.
- Nhầm $`x_1 = a + \tau(b - a)`$ (đổi vai trò hai điểm thử).
- Quên rằng bước đầu tính $`f`$ hai lần.
