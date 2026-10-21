#!/usr/bin/env python3
"""R 脚本括号/引号配平自检。

本机没有 R，无法做真正的 parse。这个脚本做的是"能做的那一半"：
逐字符扫描，正确处理字符串与 # 注释，检查 () [] {} 是否配平、
字符串是否闭合。能抓出绝大多数手误（漏括号、多括号、引号不闭合）。

用法：
    python check_syntax.py [目录]
"""

import os
import sys

PAIRS = {")": "(", "]": "[", "}": "{"}
OPENERS = "([{"
CLOSERS = ")]}"


def check_file(path):
    with open(path, encoding="utf-8") as fh:
        src = fh.read()

    stack = []          # [(char, line)]
    in_str = None       # None / '"' / "'"
    i = 0
    line = 1
    n = len(src)
    errors = []

    while i < n:
        ch = src[i]

        if ch == "\n":
            line += 1
            i += 1
            continue

        if in_str is not None:
            if ch == "\\":
                i += 2
                continue
            if ch == in_str:
                in_str = None
            i += 1
            continue

        if ch in ('"', "'"):
            in_str = ch
            i += 1
            continue

        if ch == "#":
            while i < n and src[i] != "\n":
                i += 1
            continue

        if ch in OPENERS:
            stack.append((ch, line))
        elif ch in CLOSERS:
            if not stack:
                errors.append(f"line {line}: 多余的 '{ch}'")
            elif stack[-1][0] != PAIRS[ch]:
                errors.append(
                    f"line {line}: '{ch}' 与 line {stack[-1][1]} 的 '{stack[-1][0]}' 不匹配"
                )
                stack.pop()
            else:
                stack.pop()
        i += 1

    if in_str is not None:
        errors.append(f"字符串未闭合（以 {in_str} 开始）")
    for ch, ln in stack:
        errors.append(f"line {ln}: '{ch}' 未闭合")

    return errors


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    files = sorted(
        os.path.join(root, f) for f in os.listdir(root) if f.endswith(".R")
    )
    if not files:
        print(f"未在 {root} 找到 .R 文件")
        return 1

    bad = 0
    for path in files:
        errs = check_file(path)
        name = os.path.basename(path)
        if errs:
            bad += 1
            print(f"[FAIL] {name}")
            for e in errs[:6]:
                print(f"       {e}")
        else:
            n_lines = sum(1 for _ in open(path, encoding="utf-8"))
            print(f"[ OK ] {name}  ({n_lines} 行)")

    print("-" * 46)
    print(f"共 {len(files)} 个文件，{bad} 个有问题" if bad else f"共 {len(files)} 个文件，配平检查全部通过")
    print("注意：配平通过 ≠ R 能运行。first run 请用 synthetic 模式。")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
