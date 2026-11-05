"""
Contrast check for the signature footer ("made with love and Ai by Pin0").

Opens the page in light and dark colour schemes and measures each coloured part
of the footer against the footer's real background, using the WCAG formula:
4.5:1 for the text and the letters A, i, 0; 3:1 for the heart.
It also checks that screen readers get "love" instead of the heart.

Setup:  pip install playwright
        playwright install chromium
Run:    python tests/firma_contrast.py
Exit code 1 if any check fails.
"""
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
URL = (ROOT / "index.html").as_uri()
# Each part of the footer and the minimum contrast it needs
PARTS = {
    ".firma": 4.5,
    ".firma-a": 4.5,
    ".firma-i": 4.5,
    ".firma-cero": 4.5,
    ".firma-corazon": 3.0,
}


def parse_rgb(css: str) -> tuple[float, float, float]:
    """Turn 'rgb(r, g, b)' or 'rgba(r, g, b, a)' into three 0-255 floats."""
    inner = css[css.index("(") + 1 : css.index(")")]
    r, g, b = (float(v) for v in inner.replace("/", ",").split(",")[:3])
    return r, g, b


def luminance(rgb: tuple[float, float, float]) -> float:
    def channel(v: float) -> float:
        v /= 255
        return v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4

    r, g, b = (channel(v) for v in rgb)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    la, lb = sorted((luminance(a), luminance(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


def run() -> int:
    failures = 0
    with sync_playwright() as pw:
        browser = pw.chromium.launch()
        for scheme in ("light", "dark"):
            page = browser.new_page(color_scheme=scheme)
            page.goto(URL)
            bg_css = page.eval_on_selector(".firma", "el => getComputedStyle(el).backgroundColor")
            # Reason: a transparent background would make the measurement meaningless
            if "rgba" in bg_css and bg_css.replace(" ", "").endswith(",0)"):
                print(f"FAIL {scheme}: footer background is transparent")
                failures += 1
                continue
            bg = parse_rgb(bg_css)
            for selector, minimum in PARTS.items():
                ratio = contrast(parse_rgb(page.eval_on_selector(selector, "el => getComputedStyle(el).color")), bg)
                ok = ratio >= minimum
                failures += not ok
                print(f"{'PASS' if ok else 'FAIL'} {scheme}: {selector} {ratio:.2f}:1 (needs {minimum}:1)")
            # Reason: the visually hidden "love" adds layout whitespace; screen readers collapse it
            name = " ".join(page.get_by_role("contentinfo").inner_text().split())
            ok = name == "made with love and Ai by Pin0"
            failures += not ok
            print(f"{'PASS' if ok else 'FAIL'} {scheme}: footer reads '{name}'")
            page.close()
        browser.close()
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(run())
