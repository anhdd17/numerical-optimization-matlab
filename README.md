# numerical-optimization-matlab

Bài tập lớn môn **Tối ưu hóa** (Thạc sĩ AI Engineering, học kỳ 1): trình bày lý thuyết và cài đặt bằng MATLAB/GNU Octave các thuật toán tối ưu.

Đề bài đầy đủ: [`problems/`](problems/README.md).

## Cấu trúc

```
.
├── problems/            # Đề bài (README.md) và ảnh chụp đề gốc
├── tools/               # Script kiểm tra công thức trên GitHub
└── solutions/
    └── qNN-<topic>/     # Mỗi câu một thư mục
        ├── README.md    # Báo cáo: lý thuyết, giải thuật, ví dụ
        ├── *.m          # Mã nguồn MATLAB/Octave
        └── figures/     # Hình minh họa
```

## Tiến độ

| Câu | Chủ đề | Thư mục | Trạng thái |
|---|---|---|---|
| 1 | Thuật toán đơn hình | [`q01-simplex`](solutions/q01-simplex/) | ✅ |
| 2 | Đơn hình hai pha | [`q02-two-phase-simplex`](solutions/q02-two-phase-simplex/) | ✅ |
| 3 | Phương pháp hàm phạt (bài toán M) | [`q03-big-m-method`](solutions/q03-big-m-method/) | ✅ |
| 4 | Điều kiện chính quy, KKT | [`q04-kkt-conditions`](solutions/q04-kkt-conditions/) | ✅ |
| 5 | Đối ngẫu Lagrange | [`q05-lagrange-duality`](solutions/q05-lagrange-duality/) | ✅ |
| 6 | Chia đôi | [`q06-bisection-search`](solutions/q06-bisection-search/) | ✅ |
| 7 | Lát cắt vàng | [`q07-golden-section-search`](solutions/q07-golden-section-search/) | ✅ |
| 8 | Gauss–Seidel | [`q08-gauss-seidel`](solutions/q08-gauss-seidel/) | ✅ |
| 9 | Tìm kiếm trực tiếp (Hooke–Jeeves) | [`q09-hooke-jeeves`](solutions/q09-hooke-jeeves/) | ✅ |
| 10 | Gradient và Gradient quán tính | [`q10-gradient-momentum`](solutions/q10-gradient-momentum/) | ✅ |
| 11 | Newton | [`q11-newton-method`](solutions/q11-newton-method/) | ✅ |
| 12 | Tựa Newton | [`q12-quasi-newton`](solutions/q12-quasi-newton/) | ✅ |
| 13 | Quy hoạch lồi: chính quy, KKT | [`q13-convex-programming`](solutions/q13-convex-programming/) | ✅ |
| 14 | Frank–Wolfe | [`q14-frank-wolfe`](solutions/q14-frank-wolfe/) | ✅ |
| 15 | Rosen, Murtagh–Sargent | [`q15-gradient-projection`](solutions/q15-gradient-projection/) | ✅ |
| 16 | Newton cho quy hoạch lồi ràng buộc tuyến tính | [`q16-newton-equality-constrained`](solutions/q16-newton-equality-constrained/) | ✅ |
| 17 | Hàm chắn Newton | — | ⬜ |
| 18 | Đơn hình Beale cho quy hoạch toàn phương | — | ⬜ |

## Chạy chương trình

Yêu cầu: MATLAB (R2016b trở lên) hoặc GNU Octave.

```matlab
>> cd solutions/q01-simplex
>> run_examples
```

## Kiểm tra công thức trên GitHub

Báo cáo dùng cú pháp công thức an toàn với GitHub (` ```math ` cho công thức khối, $`...`$ cho công thức inline). Sau khi viết hoặc sửa báo cáo:

```bash
python3 tools/github_math.py fix   solutions/*/README.md   # chuyển cú pháp
python3 tools/github_math.py check solutions/*/README.md   # render thử qua GitHub API
```

## Quy ước đặt tên

- Thư mục câu: `qNN-<topic>` (số có 2 chữ số, chủ đề tiếng Anh, kebab-case), ví dụ `q02-two-phase-simplex`.
- File MATLAB: `snake_case.m` (tên hàm trùng tên file); script chạy ví dụ luôn là `run_examples.m`.
- Hình: `figures/<example>_<content>.png`.
