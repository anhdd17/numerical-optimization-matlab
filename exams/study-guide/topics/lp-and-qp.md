# Quy hoạch tuyến tính và quy hoạch toàn phương (câu 1, 2, 3, 18)

[← Tài liệu ôn thi](../README.md)

Mọi số trong các câu hỏi mẫu lấy từ output của [`check_lp_qp.m`](../check_lp_qp.m) (chạy `check_lp_qp` trong thư mục `exams/study-guide`).

**Quy ước dấu ước lượng.** File này dùng quy ước của giáo trình và của chương trình câu 1–3: $`\Delta_j = c_B^T z_j - c_j`$ với $`z_j = B^{-1}A_j`$, bài toán **min**. Với quy ước này, phương án là tối ưu khi mọi $`\Delta_j \le 0`$, và $`x_j`$ vào cơ sở khi $`\Delta_j \gt 0`$. Riêng hàm `basis_check` của `toolkit/` in ra $`c_j - y^T A_j`$, tức là **ngược dấu** (xem mục 3.2 của [README](../README.md)). Gặp bài **max** $`F`$ thì đổi thành bài $`\min (-F)`$.

---

## Câu 1. Thuật toán đơn hình

### Lý thuyết cốt lõi

- **Dạng chính tắc:** $`\min c^T x`$ với $`Ax = b`$, $`x \ge 0`$, trong đó $`A`$ cỡ $`m \times n`$, $`\mathrm{rank}\,A = m`$ và $`b \ge 0`$ (dòng nào có $`b_i \lt 0`$ thì nhân với $`-1`$). Ràng buộc $`\le`$ thì cộng biến bù, $`\ge`$ thì trừ biến bù. Biến tự do $`x_j`$ thì thay bằng $`x_j^+ - x_j^-`$.
- **Cơ sở** $`B`$ gồm $`m`$ cột độc lập tuyến tính. **Phương án cực biên** là $`x_B = B^{-1}b \ge 0`$, $`x_N = 0`$; về hình học, đó là một **đỉnh** của đa diện $`\{Ax = b,\ x \ge 0\}`$. Phương án cực biên **suy biến** khi có biến cơ sở bằng 0.
- **Định lý cơ bản:** nếu bài toán có phương án tối ưu thì có phương án cực biên tối ưu. Số đỉnh hữu hạn, không quá $`\binom{n}{m}`$.
- **Công thức then chốt:** khi tăng biến phi cơ sở $`x_s`$ từ 0 lên $`\theta`$, ta có $`x_B = B^{-1}b - \theta z_s`$ và $`f = f_0 - \theta \Delta_s`$. Từ đó suy ra ba định lý:

| Định lý | Điều kiện | Kết luận |
|---|---|---|
| Dấu hiệu tối ưu | mọi $`\Delta_j \le 0`$ | phương án hiện tại tối ưu |
| Dấu hiệu không bị chặn | có $`\Delta_s \gt 0`$ và mọi $`z_{is} \le 0`$ | $`f \to -\infty`$ theo hướng $`d`$ với $`d_s = 1`$, $`d_B = -z_s`$ |
| Cải tiến | có $`\Delta_s \gt 0`$ và có $`z_{is} \gt 0`$ | xoay thì $`f`$ không tăng (giảm thực sự nếu $`\theta \gt 0`$) |

- **Hữu hạn:** nếu không có suy biến thì $`\theta \gt 0`$, $`f`$ giảm thực sự ở mỗi bước, nên không lặp lại đỉnh nào. Khi có suy biến thì quy tắc Dantzig có thể **xoay vòng** (ví dụ của Beale); quy tắc **Bland** (chọn chỉ số nhỏ nhất cho cả biến vào và biến ra) bảo đảm thuật toán dừng.

### Thuật toán

1. Lập bảng $`T = B^{-1}[A \mid b]`$ từ cơ sở chấp nhận được ban đầu (thường là các biến bù).
2. Tính $`\Delta_j = c_B^T z_j - c_j`$. Nếu mọi $`\Delta_j \le 0`$ thì dừng: tối ưu, $`f^* = c_B^T x_B`$.
3. Chọn biến vào $`s`$: theo Dantzig là $`\Delta_s`$ lớn nhất, theo Bland là chỉ số nhỏ nhất có $`\Delta_s \gt 0`$.
4. Nếu mọi $`z_{is} \le 0`$ thì dừng: không bị chặn dưới.
5. Chọn biến ra theo tỉ số nhỏ nhất $`\theta = \min\{x_{B_i}/z_{is} : z_{is} \gt 0\}`$; phần tử trục là $`z_{rs}`$.
6. Xoay: chia dòng trục cho $`z_{rs}`$, rồi khử cột $`s`$ ở các dòng khác. Quay lại bước 2.

### Dạng câu hỏi có thể gặp

| Dạng câu | Cách làm nhanh | Lệnh Matlab |
|---|---|---|
| Điền $`x^*`$, $`f^*`$ "bằng thuật toán đơn hình" | Đáp án là nghiệm tối ưu, không phụ thuộc phương pháp. Giải bằng solver, chạy đơn hình chỉ để có bảng | `simplex(c, A, b)`, `solve_lp(c, [], [], A, b, zeros(n,1))` |
| Biến vào / biến ra ở bước đầu | Tính $`\Delta`$ ở bảng đầu, rồi tỉ số $`\theta`$ | `simplex(c, A, b, B, struct('display', true))` |
| $`B`$ là cơ sở / chấp nhận được / tối ưu? | Thang 4 mức: $`\det B`$, $`B^{-1}b \ge 0`$, dấu ước lượng | `basis_check(A, b, c, B, 'min')` (in ngược dấu) |
| Vector $`x`$ là phương án / cực biên / tối ưu? | $`Ax = b`$, $`x \ge 0`$; cột ứng với $`x_j \gt 0`$ độc lập tuyến tính; ước lượng | `point_check(A, b, c, x, 'min')` |
| Bài toán có bị chặn không? | Tìm cột có $`\Delta_s \gt 0`$ mà $`z_s \le 0`$ | `simplex` trả `exitflag = -3` cùng `out.direction` |
| Phương án tối ưu có suy biến không? | Đếm số biến cơ sở bằng 0 | `x(out.basis)` |

### Câu hỏi mẫu

**Mẫu 1.1 (điền).** Dùng thuật toán đơn hình giải $`\max F = 3x_1 + 5x_2`$ với $`x_1 \le 4`$, $`2x_2 \le 12`$, $`3x_1 + 2x_2 \le 18`$, $`x \ge 0`$. Cho biết biến vào, biến ra ở bước 1 (quy tắc Dantzig) và nghiệm tối ưu.

<details><summary>Đáp án</summary>

Thêm biến bù $`x_3, x_4, x_5`$ rồi giải $`\min -F`$, cơ sở ban đầu $`\{x_3, x_4, x_5\}`$.
- Bước 0: $`\Delta = (3, 5, 0, 0, 0)`$, nên $`x_2`$ vào. Tỉ số $`\theta = \min(12/2,\ 18/2) = 6`$, nên $`x_4`$ ra (phần tử trục 2).
- Bước 1: $`\Delta_1 = 3 \gt 0`$, nên $`x_1`$ vào. Tỉ số $`\theta = \min(4/1,\ 6/3) = 2`$, nên $`x_5`$ ra.
- Bước 2: $`\Delta = (0, 0, 0, -1.5, -1) \le 0`$, tối ưu.

**Đáp án:** $`x^* = (2, 6)`$ (đầy đủ $`(2, 6, 2, 0, 0)`$), $`F^* = 36`$, sau 2 bước. `glpk` cho cùng kết quả.

</details>

**Mẫu 1.2 (trắc nghiệm).** Bài toán $`\min f = -x_1 - 2x_2`$ với $`x_1 - x_2 + x_3 = 2`$, $`-x_1 + x_2 + x_4 = 3`$, $`x \ge 0`$:
A. Có phương án tối ưu $`x = (0, 3, 5, 0)`$. B. Không có phương án chấp nhận được. C. Hàm mục tiêu không bị chặn dưới. D. Có vô số phương án tối ưu.

<details><summary>Đáp án</summary>

Bước 0 (cơ sở $`\{x_3, x_4\}`$): $`x_2`$ vào ($`\Delta_2 = 2`$), $`x_4`$ ra ($`\theta = 3`$). Bước 1: phương án $`(0, 3, 5, 0)`$ với $`f = -6`$, nhưng $`\Delta_1 = 3 \gt 0`$ và cột $`z_1 = (0, -1)^T \le 0`$. Hướng $`d = (1, 1, 0, 0)`$ thỏa $`Ad = 0`$ và $`c^Td = -3 \lt 0`$, nên $`f \to -\infty`$. **Đáp án C.** Chương trình trả `exitflag = -3`; `glpk` không trả về nghiệm tối ưu.

</details>

**Mẫu 1.3 (trắc nghiệm).** Bài toán $`\min f = -3x_1 - 2x_2`$ với $`x_1 + x_2 + x_3 = 4`$, $`x_1 + 3x_2 + x_4 = 6`$, $`x_1 + x_5 = 3`$, $`x \ge 0`$ có phương án tối ưu $`x^* = (3, 1, 0, 0, 0)`$. Phương án này:
A. Không suy biến. B. Suy biến. C. Không phải phương án cực biên. D. Không chấp nhận được.

<details><summary>Đáp án</summary>

Cơ sở cuối là $`\{x_2, x_4, x_1\}`$ với $`x_B = (1, 0, 3)`$: biến cơ sở $`x_4 = 0`$. Chỉ có 2 thành phần dương, ít hơn $`m = 3`$. **Đáp án B** (suy biến); $`f^* = -11`$. Về hình học, ba đường biên cùng đi qua điểm $`(3, 1)`$.

</details>

### Nếu ra tự luận

*"Trình bày thuật toán đơn hình giải QHTT dạng chính tắc, lấy ví dụ minh họa."* Dàn ý:
1. Dạng chính tắc; định nghĩa cơ sở và phương án cực biên; định lý cơ bản.
2. Công thức $`f = f_0 - \theta\Delta_s`$, rồi ba định lý (tối ưu, không bị chặn, cải tiến), mỗi định lý chứng minh 1–2 dòng từ công thức này.
3. Các bước thuật toán, bảng đơn hình, quy tắc tỉ số; tính hữu hạn và quy tắc Bland.
4. Ví dụ: giải tay mẫu 1.1 (3 bảng) hoặc ví dụ 1 của báo cáo câu 1.

Báo cáo đầy đủ: [câu 1](../../../solutions/q01-simplex/).

### Bẫy

- Bài **max**: hoặc đổi thành $`\min(-F)`$, hoặc đảo tiêu chuẩn. Không được trộn hai cách.
- Quy ước $`\Delta_j = c_j - c_B^T z_j`$ (một số sách dùng) ngược dấu với file này: hãy đọc kỹ đề dùng quy ước nào.
- Tỉ số $`\theta`$ chỉ xét các $`z_{is} \gt 0`$. Khi hai tỉ số bằng nhau, phương án tiếp theo sẽ suy biến.
- "Không bị chặn" chỉ kết luận được khi **toàn bộ** cột $`z_s \le 0`$, không phải chỉ một phần tử.

---

## Câu 2. Đơn hình hai pha

### Lý thuyết cốt lõi

- Vấn đề: ràng buộc $`=`$ hoặc $`\ge`$ không cho sẵn cột đơn vị, nên chưa có cơ sở xuất phát và cũng chưa biết bài toán có phương án hay không.
- **Pha 1:** thêm biến giả $`x^a_i`$ cho các dòng thiếu cột đơn vị, rồi giải $`\min w = \sum_i x^a_i`$ với $`Ax + x^a = b`$. Bài toán phụ luôn có phương án cực biên xuất phát ($`x = 0`$, $`x^a = b`$) và $`w \ge 0`$, nên pha 1 luôn dừng ở một phương án tối ưu.
- **Định lý:** bài toán gốc có phương án $`\iff w^* = 0`$.
  *Chứng minh:* nếu $`\bar x`$ là phương án của bài toán gốc thì $`(\bar x, 0)`$ cho $`w = 0`$; ngược lại, $`w^* = 0`$ kéo theo $`x^a = 0`$ và $`Ax^* = b`$.
- **Chuyển pha** (khi $`w^* = 0`$ mà còn biến giả trong cơ sở, với giá trị 0): nếu dòng của nó có $`z_{rj} \ne 0`$ trên một cột gốc thì xoay để đưa biến giả ra (vì $`\theta = 0`$, phần tử trục có thể âm). Nếu mọi $`z_{rj} = 0`$ trên các cột gốc thì ràng buộc đó **thừa**, bỏ dòng này.
- **Pha 2:** chạy đơn hình với hàm mục tiêu gốc $`c`$, xuất phát từ cơ sở cuối của pha 1.

### Thuật toán

1. Đổi dấu các dòng có $`b_i \lt 0`$, rồi thêm biến giả cho các dòng chưa có cột đơn vị.
2. Pha 1: chạy đơn hình cho $`\min w`$.
3. Nếu $`w^* \gt 0`$ thì dừng: bài toán vô nghiệm.
4. Đưa biến giả ra khỏi cơ sở hoặc bỏ dòng thừa, rồi bỏ các cột biến giả.
5. Pha 2: chạy đơn hình cho $`\min c^Tx`$. Kết quả là phương án tối ưu, hoặc kết luận không bị chặn.

### Dạng câu hỏi có thể gặp

| Dạng câu | Cách làm nhanh | Lệnh Matlab |
|---|---|---|
| Điền $`x^*`$, $`f^*`$ | Như câu 1: đáp án là nghiệm tối ưu | `two_phase_simplex(c, A, b)`, `solve_lp` |
| Cần bao nhiêu biến giả? | Đếm số dòng không có sẵn cột đơn vị (ràng buộc $`\le`$ đã có biến bù làm cột đơn vị) | `out.artificial` |
| Kết luận sau pha 1 | $`w^* \gt 0`$ thì vô nghiệm; $`w^* = 0`$ thì sang pha 2 | `out.phase1.w`, `exitflag = -2` |
| Biến giả kẹt trong cơ sở với dòng gốc toàn 0 | Ràng buộc thừa, bỏ dòng | `out.removedRows` |

### Câu hỏi mẫu

**Mẫu 2.1 (điền).** Giải bằng phương pháp hai pha: $`\min f = 2x_1 + 3x_2`$ với $`x_1 + x_2 \ge 4`$, $`x_1 + 3x_2 \ge 6`$, $`x \ge 0`$. Cho biết số biến giả, giá trị $`w^*`$ và nghiệm.

<details><summary>Đáp án</summary>

Dạng chính tắc: $`x_1 + x_2 - x_3 = 4`$, $`x_1 + 3x_2 - x_4 = 6`$. Không có cột đơn vị nào, nên cần **2 biến giả** $`x_5, x_6`$.
- **Pha 1** ($`\min w = x_5 + x_6`$), ban đầu $`w = 10`$. Bước 0: $`x_2`$ vào ($`\Delta = 4`$), $`x_6`$ ra ($`\theta = 2`$). Bước 1: $`x_1`$ vào ($`\Delta = 2/3`$), $`x_5`$ ra ($`\theta = 3`$). Bước 2: $`w^* = 0`$, cơ sở $`\{x_1, x_2\}`$.
- **Pha 2:** bảng đầu đã có $`\Delta = (0, 0, -1.5, -0.5) \le 0`$, nên tối ưu ngay (0 bước).

**Đáp án:** $`x^* = (3, 1)`$, $`f^* = 9`$. `glpk` cho cùng kết quả.
Kiểm tra bằng đỉnh: các đỉnh $`(0, 4)`$ và $`(6, 0)`$ cho $`f = 12`$, đỉnh $`(3, 1)`$ cho $`f = 9`$.

</details>

**Mẫu 2.2 (trắc nghiệm).** Áp dụng pha 1 cho $`\min f = x_1 + x_2`$ với $`x_1 + x_2 \ge 5`$, $`x_1 + 2x_2 \le 4`$, $`x \ge 0`$ thì được $`w^* = 1`$. Kết luận:
A. Bài toán có phương án tối ưu với $`f^* = 1`$. B. Bài toán không bị chặn. C. Bài toán không có phương án chấp nhận được. D. Ràng buộc thứ hai là thừa.

<details><summary>Đáp án</summary>

Chỉ dòng 1 cần biến giả $`x_5`$ (dòng 2 có biến bù $`x_4`$). Pha 1: $`x_1`$ vào, $`x_4`$ ra, rồi dừng với $`x_5 = 1`$, $`w^* = 1 \gt 0`$. **Đáp án C.** Lý do trực tiếp: $`x_1 + x_2 \le x_1 + 2x_2 \le 4 \lt 5`$. Chương trình trả `exitflag = -2`; `glpk` cũng không có nghiệm.

</details>

### Nếu ra tự luận

*"Trình bày thuật toán đơn hình hai pha."* Dàn ý:
1. Vấn đề cơ sở xuất phát.
2. Bài toán phụ, kèm định lý $`w^* = 0`$ và chứng minh hai chiều.
3. Xử lý biến giả còn trong cơ sở và ràng buộc thừa.
4. Sơ đồ hai pha; so sánh với phương pháp M (không cần chọn $`M`$, ổn định số).
5. Ví dụ: mẫu 2.1, hoặc ví dụ 1 của báo cáo câu 2 ($`x^* = (0.4, 1.8, 0, 1, 0)`$, $`f^* = 3.4`$).

Báo cáo đầy đủ: [câu 2](../../../solutions/q02-two-phase-simplex/).

### Bẫy

- Ràng buộc $`\ge`$ có biến bù với hệ số $`-1`$, nên cột của nó **không** phải cột đơn vị: dòng đó vẫn cần biến giả.
- Pha 1 không bao giờ kết thúc ở trạng thái "không bị chặn". Nếu bạn thấy kết quả đó thì đã tính sai.
- $`w^* = 0`$ nhưng biến giả còn trong cơ sở **không** có nghĩa là vô nghiệm; đó là trường hợp suy biến hoặc có ràng buộc thừa.

---

## Câu 3. Phương pháp hàm phạt (bài toán M)

### Lý thuyết cốt lõi

- Bài toán M: $`\min c^Tx + M\sum_i x^a_i`$ với $`Ax + x^a = b`$, $`x, x^a \ge 0`$, trong đó $`M`$ đủ lớn. Chỉ cần **một** lần chạy đơn hình.
- **Định lý** (với $`M`$ đủ lớn):

| Kết quả của bài toán M | Kết luận cho bài toán gốc |
|---|---|
| Tối ưu với $`x^{a*} = 0`$ | $`x^*`$ tối ưu |
| Tối ưu với $`x^{a*} \ne 0`$ | Vô nghiệm |
| Không bị chặn, biến giả bằng 0 | Không bị chặn dưới |
| Không bị chặn, còn biến giả dương | Không có phương án tối ưu (vô nghiệm hoặc không bị chặn) |

- **Ước lượng chứa $`M`$:** $`\Delta_j = \alpha_j M + \beta_j`$, so sánh theo **thứ tự từ điển**: xét hệ số $`\alpha_j`$ của $`M`$ trước, nếu bằng nhau mới xét $`\beta_j`$. Bảng đơn hình có hai dòng ước lượng.
- **Ngưỡng của $`M`$:** $`M`$ phải lớn hơn $`\max_i \bar y_i`$, trong đó $`\bar y = B^{-T}c_B`$ là biến đối ngẫu tối ưu, lấy trên các dòng có biến giả (chứng minh qua đối ngẫu yếu, xem câu 5). Nếu $`M`$ nhỏ hơn ngưỡng này, chương trình có thể kết luận **sai** là bài toán vô nghiệm. Ngược lại, $`M`$ quá lớn (cỡ $`10^{17}`$) làm phần $`\beta_j`$ bị mất do sai số làm tròn.

| | Hai pha | Bài toán M |
|---|---|---|
| Số lần chạy đơn hình | 2 | 1 |
| Tham số | không cần | cần $`M`$ đủ lớn, hoặc viết $`\alpha M + \beta`$ hình thức |
| Sai số số học | tốt | kém nếu $`M`$ là số rất lớn |
| Số bước | bằng nhau trên ví dụ 1 ($`3 = 2 + 1`$) | |

### Thuật toán

1. Thêm biến giả cho các dòng thiếu cột đơn vị; hệ số của biến giả trong hàm mục tiêu là $`M`$.
2. Lập bảng với hai dòng ước lượng ($`\alpha`$ và $`\beta`$).
3. Chạy đơn hình, chọn biến vào theo thứ tự từ điển của $`(\alpha_j, \beta_j)`$. Khi mọi biến giả đã ra khỏi cơ sở thì chỉ còn dùng dòng $`\beta`$.
4. Đọc kết luận theo bảng định lý ở trên.

### Dạng câu hỏi có thể gặp

| Dạng câu | Cách làm nhanh | Lệnh Matlab |
|---|---|---|
| Điền $`x^*`$ | Đáp án là nghiệm tối ưu | `big_m_simplex(c, A, b)` |
| Biến vào khi $`\Delta`$ chứa $`M`$ | So sánh hệ số của $`M`$ trước | `big_m_simplex(..., struct('display', true))` |
| $`M`$ = số cụ thể, kết luận có đúng không | So với ngưỡng $`\max \bar y_i`$, hoặc so sánh $`F`$ của phương án có biến giả dương với $`f^*`$ | `big_m_simplex(c, A, b, struct('M', 1.4))` |
| Bài toán M tối ưu, biến giả dương | Bài toán gốc vô nghiệm | `exitflag = -2`, `out.xa` |

### Câu hỏi mẫu

**Mẫu 3.1 (trắc nghiệm).** Bảng đầu của bài toán M cho bài toán ở mẫu 2.1 có $`\Delta_1 = 2M - 2`$, $`\Delta_2 = 4M - 3`$, $`\Delta_3 = -M`$, $`\Delta_4 = -M`$. Theo quy tắc Dantzig, biến vào cơ sở là:
A. $`x_1`$. B. $`x_2`$. C. $`x_3`$. D. Không biến nào, vì đã tối ưu.

<details><summary>Đáp án</summary>

So sánh hệ số của $`M`$ trước: $`4 \gt 2 \gt -1`$, nên chọn $`x_2`$. **Đáp án B.** Chương trình in "Bien vao: x2 (Delta = 4M-3)". Tiếp theo $`x_6`$ ra ($`\theta = 2`$); bước 1 thì $`x_1`$ vào ($`\Delta_1 = \tfrac23 M - 1`$), $`x_5`$ ra. Sau 2 bước được $`x^* = (3, 1)`$, $`f^* = 9`$, đúng như phương pháp hai pha.

</details>

**Mẫu 3.2 (trắc nghiệm).** Giải bài toán ở mẫu 2.1 bằng bài toán M với giá trị $`M`$ cụ thể. Giá trị nào dưới đây cho kết luận **sai** (bài toán vô nghiệm)?
A. $`M = 1.4`$. B. $`M = 1.6`$. C. $`M = 2`$. D. $`M = 100`$.

<details><summary>Đáp án</summary>

Biến đối ngẫu tối ưu tại cơ sở $`\{x_1, x_2\}`$ là $`\bar y = B^{-T}c_B = (1.5, 0.5)`$, nên ngưỡng là $`M \gt 1.5`$.
Kiểm tra trực tiếp: điểm $`x = (0, 2, 0, 0)`$ với biến giả $`x_5 = 2`$ có $`F = 6 + 2M`$. Điểm này tốt hơn $`f^* = 9`$ khi $`M \lt 1.5`$.
Kết quả chạy: $`M = 1`$ và $`M = 1.4`$ cho `exitflag = -2`, $`x = (0, 2, 0, 0)`$, biến giả $`(2, 0)`$, tức kết luận sai. $`M = 1.6`$, $`2`$, $`100`$ cho đúng $`x^* = (3, 1, 0, 0)`$, $`f^* = 9`$. **Đáp án A.**

</details>

### Nếu ra tự luận

*"Trình bày phương pháp hàm phạt giải QHTT dạng chính tắc."* Dàn ý:
1. Ý tưởng và bài toán $`(P_M)`$.
2. Định lý quan hệ giữa $`(P_M)`$ và $`(P)`$ (ba ý, chứng minh ý 1–2 bằng so sánh giá trị mục tiêu, ý 3 bằng đối ngẫu yếu).
3. Ước lượng $`\alpha M + \beta`$ và thứ tự từ điển.
4. Các bước thuật toán; so sánh với hai pha.
5. Ví dụ: mẫu 2.1 giải bằng bài toán M (3 bảng), cộng ví dụ chọn $`M`$ (mẫu 3.2).

Báo cáo đầy đủ: [câu 3](../../../solutions/q03-big-m-method/).

### Bẫy

- Nhìn dấu $`\Delta_j`$ theo hệ số của $`M`$ trước: $`-2M + 1 \lt 0`$ dù phần không chứa $`M`$ dương.
- Kết luận "vô nghiệm" chỉ đúng khi bài toán M **tối ưu** mà biến giả dương. Nếu bài toán M không bị chặn thì phải xét thêm giá trị biến giả.
- Biến giả đã ra khỏi cơ sở thì không bao giờ cần đưa vào lại; có thể xóa cột của nó.

---

## Câu 18. Đơn hình Beale cho quy hoạch toàn phương

### Lý thuyết cốt lõi

- Bài toán: $`\min f = c^Tx + \tfrac12 x^T D x`$ với $`Ax = b`$, $`x \ge 0`$, trong đó $`D`$ đối xứng nửa xác định dương (nên $`f`$ lồi). Ràng buộc affine nên KKT là điều kiện cần và đủ:

```math
c + Dx^* - A^Ty - v = 0,\qquad Ax^* = b,\qquad x^* \ge 0,\qquad v \ge 0,\qquad v_j x^*_j = 0 .
```

- Khác QHTT: nghiệm **không nhất thiết là đỉnh**; nó có thể nằm giữa một cạnh, trên một mặt, hay ở trong miền.
- **Đạo hàm theo biến phi cơ sở** $`x_j`$ (hướng cạnh $`d_j`$ có $`(d_j)_j = 1`$, $`(d_j)_B = -B^{-1}A_j`$):

```math
\frac{\partial f}{\partial x_j} = d_j^T(c + D\bar x),\qquad \frac{\partial^2 f}{\partial x_j^2} = d_j^T D d_j \ge 0 .
```

Với $`D = 0`$ thì $`\partial f/\partial x_j = -\Delta_j`$, tức là trở về đơn hình.
- **Hai tỉ số** khi tăng biến vào $`x_k`$ (chọn $`x_k`$ có $`\partial f/\partial x_k \lt 0`$):
  - $`\theta_1 = \min\{\bar x_i / (-(d_k)_i) : (d_k)_i \lt 0\}`$: một biến cơ sở giảm về 0.
  - $`\theta_2 = -\partial_k f / \partial^2_{kk} f`$: cực tiểu của $`f`$ trên cạnh.

  Nếu $`\theta_1 \le \theta_2`$ thì **xoay** như đơn hình. Nếu $`\theta_2 \lt \theta_1`$ thì **tạo biến tự do** $`u = \tfrac12\, \partial f / \partial x_k`$ (không ràng buộc dấu), thêm dòng định nghĩa $`u`$ vào hệ, $`x_k`$ vào cơ sở trên dòng mới.
- **Tiêu chuẩn tối ưu:** $`\partial f/\partial x_j \ge 0`$ với mọi $`x_j`$ phi cơ sở và $`\partial f / \partial u_t = 0`$ với mọi biến tự do phi cơ sở. Khi đó $`v_j = \partial f/\partial x_j`$ chính là **nhân tử KKT** của $`x_j \ge 0`$.
  *Chứng minh:* bất đẳng thức gradient của hàm lồi cho $`f(x) \ge f(\bar x) + \sum_{x_j} \partial_j f\, x_j + \sum_{u_t} \partial_{u_t} f\, u_t \ge f(\bar x)`$.
- Biến tự do vào cơ sở thì bị loại cùng dòng định nghĩa của nó. Phương pháp là **hữu hạn** và cho nghiệm chính xác; giữ $`u = 0`$ buộc các hướng đi liên hợp đối với $`D`$.

### Thuật toán

1. Tìm cơ sở chấp nhận được ban đầu (pha 1 của câu 2, hoặc dùng biến bù).
2. Tính $`\partial f/\partial z_j`$ cho mọi biến phi cơ sở $`z_j`$ (cả $`x_j`$ lẫn $`u_t`$).
3. Dừng nếu thỏa tiêu chuẩn tối ưu.
4. Chọn biến vào có vi phạm lớn nhất (với $`u_t`$, chiều đi là $`-\mathrm{sign}(\partial f/\partial u_t)`$).
5. Tính $`\theta_1`$ và $`\theta_2`$. Nếu cả hai bằng $`\infty`$ thì $`f`$ không bị chặn dưới.
6. Nếu $`\theta_1 \le \theta_2`$ thì xoay; ngược lại tạo biến tự do mới. Loại các biến tự do đang ở trong cơ sở. Quay lại bước 2.

### Dạng câu hỏi có thể gặp

| Dạng câu | Cách làm nhanh | Lệnh Matlab |
|---|---|---|
| Điền $`x^*`$ "bằng đơn hình Beale" | Đáp án là nghiệm của QP lồi. Giải KKT tay (thường chỉ có 1 ràng buộc tích cực), kiểm tra bằng solver | `beale_qp(c, D, A, b)`, `solve_qp(H, f, A, b, Aeq, beq, lb)` |
| $`\theta_1`$, $`\theta_2`$ ở bước đầu, xoay hay tạo $`u`$ | Tính $`\partial f/\partial x_k`$ và $`d^TDd`$ tại $`x^0`$ | `beale_qp(..., struct('display', true, 'B0', B))` |
| Nhân tử của $`x_j \ge 0`$ tại nghiệm | $`\partial f/\partial x_j`$ theo biến phi cơ sở | `out.lambda` |
| Lý thuyết: khi nào tạo biến tự do / khi nào dừng | Bảng $`\theta_1`$, $`\theta_2`$ và tiêu chuẩn tối ưu | |

Lưu ý: `beale_qp` nhận $`f = c^Tx + \tfrac12 x^TDx`$, nên $`D`$ là **Hessian** (ví dụ $`x_1^2`$ cho $`D_{11} = 2`$).

### Câu hỏi mẫu

**Mẫu 18.1 (điền).** Dùng đơn hình Beale giải $`\min f = x_1^2 + x_2^2 - 4x_1 - 2x_2`$ với $`x_1 + x_2 + x_3 = 2`$, $`x \ge 0`$, xuất phát từ cơ sở $`\{x_3\}`$.

<details><summary>Đáp án</summary>

Ta có $`c = (-4, -2, 0)`$ và $`D = \mathrm{diag}(2, 2, 0)`$.
- **Bước 0:** $`\bar x = (0, 0, 2)`$, $`\partial f/\partial x_1 = -4`$, $`\partial f/\partial x_2 = -2`$, nên chọn $`x_1`$. Vì $`x_3 = 2 - x_1`$ nên $`\theta_1 = 2`$; $`f(\theta) = \theta^2 - 4\theta`$ nên $`\theta_2 = 2`$. Do $`\theta_1 \le \theta_2`$: **xoay**, $`x_1`$ vào, $`x_3`$ ra, được $`(2, 0, 0)`$ với $`f = -4`$.
- **Bước 1:** $`x_1 = 2 - x_2 - x_3`$. Hướng của $`x_2`$ là $`d = (-1, 1, 0)`$, cho $`\partial f/\partial x_2 = -1\cdot(2\cdot2 - 4) + (0 - 2) = -2`$ và $`d^TDd = 4`$. Suy ra $`\theta_2 = 0.5 \lt \theta_1 = 2`$: **tạo biến tự do** $`u_1`$, được $`(1.5, 0.5, 0)`$.
- **Bước 2:** tối ưu. Nhân tử $`\partial f/\partial x_3 = 1`$.

**Đáp án:** $`x^* = (1.5, 0.5)`$, $`f^* = -4.5`$, sau 2 bước.
Kiểm tra bằng KKT tay: $`2x_1 - 4 + \mu = 0`$, $`2x_2 - 2 + \mu = 0`$, $`x_1 + x_2 = 2`$ cho $`\mu = 1 \ge 0`$. Solver `qp` cho $`(1.5, 0.5)`$, $`f = -4.5`$.

</details>

**Mẫu 18.2 (trắc nghiệm).** Trong một bước của giải thuật Beale, tăng biến vào $`x_k`$ có $`\theta_1 = 1.25`$ (biến cơ sở về 0) và $`\theta_2 = 1.129`$ (cực tiểu của $`f`$ trên cạnh). Bước tiếp theo là:
A. Xoay như đơn hình, $`x_k`$ vào và biến cơ sở chặn đi ra. B. Tạo biến tự do $`u = \tfrac12\, \partial f / \partial x_k`$, $`x_k`$ vào cơ sở trên dòng mới. C. Dừng, phương án hiện tại tối ưu. D. Kết luận $`f`$ không bị chặn dưới.

<details><summary>Đáp án</summary>

$`\theta_2 \lt \theta_1`$: $`f`$ đạt cực tiểu trên cạnh trước khi chạm biên, nên dừng giữa cạnh và tạo biến tự do. **Đáp án B.**
Đây là bước 1 của ví dụ 1 trong câu 18: $`\min -4x_1 - 6x_2 + 2x_1^2 + 2x_2^2 - 2x_1x_2`$ với $`x_1 + x_2 + x_3 = 2`$, $`x_1 + 5x_2 + x_4 = 5`$. Bước 0: xoay ($`\theta_1 = 1 \lt \theta_2 = 1.5`$). Bước 1: tạo $`u_1`$. Nghiệm $`x^* = (35/31, 24/31, 3/31, 0)`$, $`f^* = -222/31 = -7.161290`$, sau 2 bước với 1 biến tự do.

</details>

### Nếu ra tự luận

*"Trình bày giải thuật đơn hình Beale giải quy hoạch toàn phương."* Dàn ý:
1. Bài toán QP lồi, điều kiện KKT; giải thích vì sao nghiệm không phải luôn là đỉnh.
2. Biểu diễn qua biến phi cơ sở, công thức đạo hàm riêng.
3. Bước Beale với $`\theta_1`$, $`\theta_2`$; biến tự do và quy tắc loại.
4. Tiêu chuẩn tối ưu, có chứng minh bằng bất đẳng thức gradient; nhân tử $`v_j = \partial f/\partial x_j`$.
5. Tính hữu hạn (hướng liên hợp).
6. Ví dụ: mẫu 18.1 (2 bước, gọn) hoặc ví dụ 1 của báo cáo câu 18.

Báo cáo đầy đủ: [câu 18](../../../solutions/q18-beale-quadratic-programming/).

### Bẫy

- $`D`$ là Hessian: số hạng $`2x_1^2`$ cho $`D_{11} = 4`$, còn $`-2x_1x_2`$ cho $`D_{12} = D_{21} = -2`$.
- $`\theta_1 = \theta_2`$ thì **xoay** (điều kiện là $`\theta_1 \le \theta_2`$), như bước 0 của mẫu 18.1.
- Với biến tự do, tiêu chuẩn tối ưu là đạo hàm **bằng 0**, không phải $`\ge 0`$.
- Nếu đề cho QP dạng bất đẳng thức thì phải thêm biến bù trước khi lập bảng. Khi đối chiếu bằng `solve_qp`, dùng dạng bất đẳng thức với ít biến: ở câu 18, `qp` đã trả kết quả sai với Hessian suy biến $`4 \times 4`$ ở dạng chính tắc.
