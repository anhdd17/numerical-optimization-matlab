-- Bo loc pandoc cho ban PDF offline:
--  * so do Mermaid (khong render duoc offline) -> danh sach cac buoc
--  * <details> -> <details open> de dap an an hien ra khi in
local flow = [[
1. **Đọc đề**: không ràng buộc / một biến / QHTT / có ràng buộc?
2. **Bài toán lồi (chặt)?**
   - Có → tìm nghiệm nhanh: `-H\q`, KKT, `fminunc`, `linprog`/`glpk`, `fmincon`/`sqp`.
   - Không / một biến → vẽ đồ thị, xét dấu f', thử nhiều điểm xuất phát.
3. **Chạy phương pháp đề yêu cầu** (`toolkit/` hoặc `solutions/`), ghi bảng lặp.
4. **Khớp đến 4 chữ số?**
   - Có → ghi đáp án + tiêu chuẩn dừng.
   - Không → kiểm tra lại gradient, hệ số, dấu ràng buộc, quay lại bước 1.
]]

function CodeBlock(el)
  if el.classes:includes('mermaid') then
    return pandoc.read(flow, 'gfm').blocks
  end
end

local function open_details(el)
  if el.format == 'html' and el.text:match('^%s*<details') and not el.text:match('<details[^>]*open') then
    el.text = el.text:gsub('<details', '<details open', 1)
    return el
  end
end
RawBlock = open_details
RawInline = open_details

-- Link tuong doi viet tu topics/*.md lech mot cap khi gop vao study-guide.html
-- (dat o study-guide/): neu dich khong ton tai thi bo bot mot "../".
-- Link toi README.md cua tai lieu -> neo noi bo.
local function exists(p) return os.rename(p, p) ~= nil end
function Link(el)
  local t = el.target
  if t:match('^%a+:') or t:match('^#') then return nil end
  local path, frag = t:match('^([^#]*)(#?.*)$')
  if path ~= '' and not exists(path) and path:match('^%.%./') and exists(path:sub(4)) then
    path = path:sub(4)
  end
  if path == 'README.md' then
    el.target = (frag ~= '' and ('#readmemd__' .. frag:sub(2))) or '#TOC'
  else
    el.target = path .. frag
  end
  return el
end

-- Link toi ca mot file (#topics__lp-and-qpmd) -> tieu de cap 1 cua file do.
function Pandoc(doc)
  local first = {}
  doc.blocks:walk({ Header = function(h)
    local pre = h.identifier:match('^(.*md)__')
    if h.level == 1 and pre and not first[pre] then first[pre] = h.identifier end
  end })
  return doc:walk({ Link = function(l)
    local k = l.target:match('^#(.*md)$')
    if k and first[k] then l.target = '#' .. first[k]; return l end
  end })
end
