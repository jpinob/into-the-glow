"""
Build the page: insert data/gfp.json into src/template.html and write index.html.

Run from the repository root:  python pipeline/build.py
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
template = (ROOT / "src" / "template.html").read_text(encoding="utf-8")
data = (ROOT / "data" / "gfp.json").read_text(encoding="utf-8")
assert "__DATA__" in template, "placeholder __DATA__ not found in src/template.html"
(ROOT / "index.html").write_text(template.replace("__DATA__", data), encoding="utf-8")
print("index.html written")
