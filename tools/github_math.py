#!/usr/bin/env python3
"""Make LaTeX math in Markdown files render correctly on GitHub, and verify it.

Usage:
    python3 tools/github_math.py fix   solutions/qNN-*/README.md
    python3 tools/github_math.py check solutions/qNN-*/README.md

fix   - rewrites the files in place:
        * $$ ... $$ blocks (on their own lines)  -> ```math fenced blocks
        * every inline $...$                     -> $`...`$
        * '<' and '>' inside math                -> \\lt and \\gt
check - renders each file through GitHub's Markdown API and verifies that every
        formula in the source comes back as a math element, unchanged.
        Exits with status 1 if any formula is lost or altered.

Why: GitHub runs its Markdown parser before MathJax. In plain $...$ / $$...$$ it
eats backslash escapes (\\; \\, \\{ \\} \\\\), pairs '*' (e.g. x^*) into emphasis
across formulas, and double-escapes '<' / '>'. ```math and $`...`$ are protected.
Disallowed macros (\\operatorname, \\DeclareMathOperator, \\newcommand, \\def) are
reported by `check` as well; use \\mathrm{...} instead.
"""
import html
import json
import re
import sys
import urllib.request

INLINE = re.compile(r'(?<![\$`\\])\$(?![\$`])([^$\n]+?)(?<![`\\])\$(?!\$)')
BANNED = re.compile(r'\\(operatorname|DeclareMathOperator|newcommand|renewcommand|def)\b')


def lt_gt(tex):
    tex = re.sub(r'\s*>\s*', r' \\gt ', tex)
    return re.sub(r'\s*<\s*', r' \\lt ', tex)


def fix(path):
    lines, out = open(path, encoding='utf-8').read().split('\n'), []
    fence = in_dollar_block = False
    n_block = n_inline = 0
    for line in lines:
        s = line.strip()
        if s.startswith('```') and not in_dollar_block:
            fence = not fence
            out.append(line)
            continue
        if fence:
            out.append(line)
            continue
        if s == '$$':
            out.append('```' if in_dollar_block else '```math')
            n_block += not in_dollar_block
            in_dollar_block = not in_dollar_block
            continue
        if in_dollar_block:
            out.append(line)
            continue
        # leave inline code spans and already-converted $`...`$ untouched
        parts = re.split(r'(\$`.*?`\$|`[^`]*`)', line)
        for j in range(0, len(parts), 2):
            parts[j], k = INLINE.subn(lambda m: '$`' + m.group(1) + '`$', parts[j])
            n_inline += k
        out.append(''.join(parts))
    text = '\n'.join(out)
    text = re.sub(r'(```math\n)(.*?)(```)', lambda m: m.group(1) + lt_gt(m.group(2)) + m.group(3),
                  text, flags=re.S)
    text = re.sub(r'\$`(.*?)`\$', lambda m: '$`' + lt_gt(m.group(1)) + '`$', text)
    open(path, 'w', encoding='utf-8').write(text)
    print(f'{path}: {n_block} display blocks, {n_inline} inline formulas converted')


def norm(tex):
    tex = re.sub(r'^\$\$|\$\$$|^\$`?|`?\$$', '', tex.strip()).strip()
    return tex.replace('\\\\\\\n', '\\\\\n')  # GitHub adds a '\' before newlines that follow '\\'


def check(path):
    src = open(path, encoding='utf-8').read()
    req = urllib.request.Request(
        'https://api.github.com/markdown',
        data=json.dumps({'text': src, 'mode': 'gfm'}).encode(),
        headers={'Accept': 'application/vnd.github+json', 'Content-Type': 'application/json'})
    rendered = urllib.request.urlopen(req).read().decode()
    got = {norm(html.unescape(m)) for m in
           re.findall(r'<math-renderer[^>]*>(.*?)</math-renderer>', rendered, re.S)}
    blocks = [b.strip() for b in re.findall(r'```math\n(.*?)```', src, re.S)]
    prose = re.sub(r'```.*?```', '', src, flags=re.S)
    inline = [norm(p) for p in re.findall(r'\$`.*?`\$', prose)]
    bare = INLINE.findall(prose) + re.findall(r'^\s*\$\$\s*$', prose, re.M)
    missing = [p for p in blocks + inline if p not in got]
    banned = sorted(set(BANNED.findall(src)))
    ok = not missing and not bare and not banned
    print(f'{path}: {len(blocks)} blocks + {len(inline)} inline; '
          f'not rendered: {len(missing)}; bare $...$/$$: {len(bare)}; banned macros: {banned or 0}'
          f' -> {"OK" if ok else "FAIL"}')
    for p in missing[:10]:
        print('    not rendered:', p[:100].replace('\n', ' '))
    return ok


if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('fix', 'check'):
        sys.exit(__doc__)
    if sys.argv[1] == 'fix':
        for f in sys.argv[2:]:
            fix(f)
    else:
        results = [check(f) for f in sys.argv[2:]]
        sys.exit(0 if all(results) else 1)
