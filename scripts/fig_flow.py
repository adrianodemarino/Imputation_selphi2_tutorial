"""The tutorial at a glance (img/flow_{dark,light}.svg): the five steps and the files each one makes.
A highlight walks through the steps in order. SMIL animation, which GitHub plays in a README.
usage: fig_flow.py OUT_DIR"""
import sys

W, H = 1000, 250
T = 10.0
FONT = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif"
MONO = "ui-monospace, SFMono-Regular, Menlo, Consolas, monospace"
THEMES = {
    "dark": dict(bg="#0d1117", box="#161b22", line="#30363d", fg="#e6edf3", sub="#8b949e",
                 steps=["#8b949e", "#8b949e", "#58a6ff", "#d29922", "#3fb950"]),
    "light": dict(bg="#ffffff", box="#f6f8fa", line="#d0d7de", fg="#1f2328", sub="#59636e",
                  steps=["#59636e", "#59636e", "#0969da", "#9a6700", "#1a7f37"]),
}
STEPS = [
    ("1", "Download", ["1000 Genomes panel", "genetic map"]),
    ("2", "Prepare", ["panel  3,102 people", "array  100 people", "truth  30x"]),
    ("3", "Selphi 2", [".srp2 panel", "imputed VCF", "accuracy"]),
    ("4", "Beagle 5.5", ["optional:", "same job,", "to compare"]),
    ("5", "1x reads", ["optional:", "impute from", "a BAM/CRAM"]),
]
BW, BH, GAP, X0, Y0 = 168, 150, 30, 20, 56


def show_glow(k):
    """Highlight step k during its slot, every loop."""
    a, b = 0.4 + 1.8 * k, 0.4 + 1.8 * (k + 1)
    ka, kb = a / T, b / T
    return (f'<animate attributeName="opacity" values="0;0;1;1;0;0" '
            f'keyTimes="0;{ka:.4f};{ka + 0.03:.4f};{kb - 0.02:.4f};{kb + 0.01:.4f};1" dur="{T}s" repeatCount="indefinite"/>')


def flow(name):
    t = THEMES[name]
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" '
         f'aria-label="The tutorial in five steps: download, prepare, Selphi 2, Beagle (optional), 1x reads (optional)">',
         f'<rect width="{W}" height="{H}" rx="12" fill="{t["bg"]}"/>',
         f'<text x="{X0}" y="34" font-family="{FONT}" font-size="15" font-weight="700" fill="{t["fg"]}">'
         f'The tutorial in five steps</text>']
    for k, (num, title, files) in enumerate(STEPS):
        x = X0 + k * (BW + GAP)
        col = t["steps"][k]
        dash = ' stroke-dasharray="5 4"' if k >= 3 else ""
        o.append(f'<rect x="{x}" y="{Y0}" width="{BW}" height="{BH}" rx="10" fill="{t["box"]}" '
                 f'stroke="{t["line"]}" stroke-width="1.5"{dash}/>')
        o.append(f'<rect x="{x}" y="{Y0}" width="{BW}" height="{BH}" rx="10" fill="none" stroke="{col}" '
                 f'stroke-width="2.5" opacity="0">{show_glow(k)}</rect>')
        o.append(f'<circle cx="{x + 24}" cy="{Y0 + 28}" r="13" fill="{col}"/>')
        o.append(f'<text x="{x + 24}" y="{Y0 + 33}" text-anchor="middle" font-family="{FONT}" font-size="13" '
                 f'font-weight="700" fill="{t["bg"]}">{num}</text>')
        o.append(f'<text x="{x + 46}" y="{Y0 + 33}" font-family="{FONT}" font-size="15" font-weight="700" '
                 f'fill="{t["fg"]}">{title}</text>')
        for i, f in enumerate(files):
            o.append(f'<text x="{x + 16}" y="{Y0 + 70 + 22 * i}" font-family="{MONO}" font-size="12" '
                     f'fill="{t["sub"]}">{f}</text>')
        if k < len(STEPS) - 1:
            ax = x + BW + 4
            o.append(f'<path d="M{ax} {Y0 + BH / 2} h{GAP - 10} m-6 -5 l6 5 l-6 5" fill="none" '
                     f'stroke="{t["line"]}" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>')
    o.append(f'<text x="{W - 20}" y="{H - 16}" text-anchor="end" font-family="{FONT}" font-size="12" '
             f'fill="{t["sub"]}">dashed = optional</text>')
    o.append("</svg>")
    return "\n".join(o) + "\n"


if __name__ == "__main__":
    for name in THEMES:
        with open(f"{sys.argv[1]}/flow_{name}.svg", "w") as f:
            f.write(flow(name))
