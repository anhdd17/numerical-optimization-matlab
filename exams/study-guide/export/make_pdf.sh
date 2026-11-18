#!/usr/bin/env bash
# Xuat tai lieu on thi thanh 1 file HTML tu chua (cong thuc MathML, khong can mang)
# va 1 file PDF de doc tren may thi. Can: pandoc >= 3, Google Chrome.
#   bash export/make_pdf.sh        (chay tu exams/study-guide)
set -euo pipefail
cd "$(dirname "$0")/.."
CHROME=${CHROME:-"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"}
FILES=(README.md topics/lp-and-qp.md topics/kkt-duality-convex.md topics/line-search.md \
       topics/unconstrained.md topics/constrained-methods.md practice.md)
pandoc "${FILES[@]}" -f gfm -t html5 -s --mathml --file-scope \
  --lua-filter export/offline.lua --css export/print.css --embed-resources \
  --toc --toc-depth=2 --metadata title="Tài liệu ôn thi Tối ưu hóa" --metadata lang=vi \
  -o study-guide.html
"$CHROME" --headless=new --disable-gpu --no-pdf-header-footer \
  --print-to-pdf="$PWD/study-guide.pdf" "file://$PWD/study-guide.html" 2>/dev/null
ls -l study-guide.html study-guide.pdf
