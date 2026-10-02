"""
Visual check: open the page in a headless browser and save a screenshot of
every stop, on a desktop-sized and a phone-sized screen.

Setup:  pip install playwright
        playwright install chromium
Run:    python tests/screenshots.py                     (local index.html)
        python tests/screenshots.py https://jpinob.github.io/into-the-glow/
Output: tests/output/desktop_stop1.png ... tests/output/phone_stop6.png
Exit code 1 if the page throws a JavaScript error.
"""
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
URL = sys.argv[1] if len(sys.argv) > 1 else (ROOT / "index.html").as_uri()
OUT = ROOT / "tests" / "output"
VIEWPORTS = {
    "desktop": {"width": 1280, "height": 800},
    "phone": {"width": 390, "height": 844},
}
STOPS = 6
# Software WebGL so the check also works on machines without a usable GPU
LAUNCH_ARGS = ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"]


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch(args=LAUNCH_ARGS)
        for name, viewport in VIEWPORTS.items():
            page = browser.new_page(viewport=viewport, has_touch=(name == "phone"))
            page.on("pageerror", lambda e, n=name: errors.append(f"{n}: {e}"))
            page.goto(URL)
            page.wait_for_timeout(4000)
            for stop in range(1, STOPS + 1):
                if stop > 1:
                    page.click("#next")
                    page.wait_for_timeout(4500)
                if stop == 5:
                    page.click("#shine")
                    page.wait_for_timeout(3000)
                page.screenshot(path=str(OUT / f"{name}_stop{stop}.png"))
            page.close()
        browser.close()
    print(f"Screenshots saved in {OUT.relative_to(ROOT)}")
    if errors:
        print("Page errors:", *errors, sep="\n  ")
        sys.exit(1)


if __name__ == "__main__":
    run()
