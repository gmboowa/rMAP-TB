"""Figure S1: the three report states the test-suite asserts.

Each panel is a screenshot of the phylogeny section of a report built by
gen_states.py from synthetic inputs, cropped from the section header to the
last element that carries the state's evidence. The sticky navigation bar is
hidden so that it does not overlay the section heading.

Run:  python3 fig_s1.py <out_dir>
"""
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from playwright.sync_api import sync_playwright

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "out")
SCALE = 2
WIDTH = 1180

# file, last child selector to include, panel letter, caption
PANELS = [
    ("state_A_rendering_failure.html", ".tree-failure", "A",
     "Rendering failed: the error, the renderer status keys and the log are "
     "shown in place of an empty image."),
    ("state_B_annotations_incomplete.html", ".tree-panel", "B",
     "Rendered, annotations incomplete: the topology is valid but no tip "
     "matched the metadata, and the report says which section to use instead."),
    ("state_C_legacy_placeholder.html", ".tree-failure", "C",
     "Legacy 68-byte placeholder image: recognised as unusable and not "
     "embedded as if it were a tree."),
]

HIDE = """
  document.querySelectorAll('nav, .sec-nav, .topnav, .rm-nav, [class*="sticky"]')
    .forEach(e => { e.style.position = 'static'; e.style.display = 'none'; });
  document.querySelectorAll('*').forEach(e => {
    const p = getComputedStyle(e).position;
    if (p === 'sticky' || p === 'fixed') e.style.display = 'none';
  });
"""

shots = []
with sync_playwright() as pw:
    br = pw.chromium.launch()
    pg = br.new_page(viewport={"width": WIDTH, "height": 1400},
                     device_scale_factor=SCALE)
    for fname, last_sel, tag, _cap in PANELS:
        pg.goto((OUT / fname).resolve().as_uri())
        pg.wait_for_timeout(900)
        pg.evaluate(HIDE)
        pg.wait_for_timeout(250)
        box = pg.evaluate("""(sel) => {
            const sec = document.querySelector('#sec-tree');
            const head = sec.querySelector('.sec-head');
            const last = sec.querySelector(sel);
            const a = head.getBoundingClientRect();
            const b = last.getBoundingClientRect();
            const y = a.top + window.scrollY;
            return {x: Math.floor(a.left), y: Math.floor(y),
                    width: Math.ceil(a.width),
                    height: Math.ceil(b.bottom + window.scrollY - y)};
        }""", last_sel)
        path = OUT / ("s1_%s.png" % tag)
        pg.screenshot(path=str(path), clip=box, full_page=True)
        im = Image.open(path).convert("RGB")
        shots.append((tag, im))
        print("%s  %-42s %d x %d px" % (tag, fname, im.width, im.height))
    br.close()

# ------------------------------------------------------------- compose
PAD, GAP, LABEL_W = 26, 34, 70
try:
    f_lab = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 46)
    f_cap = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 25)
except OSError:
    f_lab = f_cap = ImageFont.load_default()

inner = max(im.width for _t, im in shots)
W = LABEL_W + inner + 2 * PAD


def wrap(draw, text, font, maxw):
    words, lines, cur = text.split(), [], ""
    for w in words:
        t = (cur + " " + w).strip()
        if draw.textlength(t, font=font) <= maxw:
            cur = t
        else:
            lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    return lines


probe = ImageDraw.Draw(Image.new("RGB", (10, 10)))
caps = [wrap(probe, c, f_cap, inner) for _f, _s, _t, c in PANELS]
cap_h = [len(l) * 32 + 12 for l in caps]

H = PAD + sum(im.height + cap_h[i] + GAP for i, (_t, im) in enumerate(shots)) - GAP + PAD
canvas = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(canvas)

y = PAD
for i, (tag, im) in enumerate(shots):
    d.text((PAD, y - 4), tag, font=f_lab, fill="#16202b")
    canvas.paste(im, (LABEL_W + PAD, y))
    yy = y + im.height + 10
    for line in caps[i]:
        d.text((LABEL_W + PAD, yy), line, font=f_cap, fill="#39424e")
        yy += 32
    y += im.height + cap_h[i] + GAP

png = OUT / "Figure_S1_report_states.png"
canvas.save(png, dpi=(600, 600))
canvas.save(OUT / "Figure_S1_report_states.pdf", "PDF", resolution=600.0)
print("\nFigure S1: %d x %d px (%.1f x %.1f mm at 600 dpi)"
      % (canvas.width, canvas.height,
         canvas.width / 600 * 25.4, canvas.height / 600 * 25.4))
