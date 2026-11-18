"""Result charts (img/r2_*, img/speed_*, img/lowcov_*; dark and light) from results/summary.json.
Lines draw and bars grow once when the page opens (SMIL, fill=freeze).
usage: fig_results.py results/summary.json OUT_DIR"""
import json
import sys

FONT = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif"
BINS = ["0.5-1%", "1-2%", "2-5%", "5-10%", "10-20%", "20-50%"]
THEMES = {
    "light": dict(bg="#ffffff", fg="#1f2328", sub="#59636e", grid="#e6eaef", axis="#afb8c1",
                  selphi="#2a78d6", beagle="#eb6834", cov=["#a9c6ee", "#5f9be0", "#1f5fae"]),
    "dark": dict(bg="#0d1117", fg="#e6edf3", sub="#8b949e", grid="#21262d", axis="#484f58",
                 selphi="#3987e5", beagle="#d95926", cov=["#1f4f8f", "#3987e5", "#9cc3f5"]),
}


def grow(attr, frm, to, delay=0.0, dur=1.1):
    """Animate attr from `frm` to `to` once on load. The element's own attribute holds the final value, so a
    renderer without SMIL still shows the finished chart; the animation holds `frm` during `delay`."""
    total = delay + dur
    if delay > 0:
        vals, kt, spl = f"{frm};{frm};{to}", f"0;{delay / total:.4f};1", "0 0 1 1;.2 .7 .2 1"
    else:
        vals, kt, spl = f"{frm};{to}", "0;1", ".2 .7 .2 1"
    return (f'<animate attributeName="{attr}" values="{vals}" keyTimes="{kt}" calcMode="spline" '
            f'keySplines="{spl}" dur="{total:.2f}s"/>')


def frame(t, w, h, title, subtitle):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="{FONT}">',
            f'<rect width="{w}" height="{h}" rx="12" fill="{t["bg"]}"/>',
            f'<text x="24" y="32" font-size="16" font-weight="700" fill="{t["fg"]}">{title}</text>',
            f'<text x="24" y="52" font-size="12" fill="{t["sub"]}">{subtitle}</text>']


def line_chart(t, series, title, subtitle, w=920, h=380):
    """series: [(label, colour, [r2 per bin or None], dashed)]"""
    x0, x1, y0, y1 = 70, w - 150, h - 56, 80          # plot box (y0 bottom)
    o = frame(t, w, h, title, subtitle)
    keep = [i for i in range(len(BINS)) if any(v[i] is not None for _, _, v, _ in series)]
    bins = [BINS[i] for i in keep]
    series = [(a, b, [v[i] for i in keep], d) for a, b, v, d in series]
    xs = [x0 + (x1 - x0) * i / (len(bins) - 1) for i in range(len(bins))]
    allv = [v for _, _, vals, _ in series for v in vals if v is not None]
    lo = max(0.0, int(min(allv) * 10) / 10 - 0.1)
    hi = min(1.0, int(max(allv) * 10 + 1) / 10 + 0.05)
    ys = lambda v: y0 - (y0 - y1) * (v - lo) / (hi - lo)
    ticks = [round(lo + 0.1 * i, 1) for i in range(int(round((hi - lo) / 0.1)) + 1) if lo + 0.1 * i <= hi + 1e-9]
    for v in ticks:
        o.append(f'<line x1="{x0}" x2="{x1}" y1="{ys(v):.1f}" y2="{ys(v):.1f}" stroke="{t["grid"]}"/>')
        o.append(f'<text x="{x0 - 10}" y="{ys(v) + 4:.1f}" text-anchor="end" font-size="11" fill="{t["sub"]}">{v:g}</text>')
    for x, b in zip(xs, bins):
        o.append(f'<text x="{x:.1f}" y="{y0 + 20}" text-anchor="middle" font-size="11" fill="{t["sub"]}">{b}</text>')
    o.append(f'<text x="{(x0 + x1) / 2:.0f}" y="{h - 12}" text-anchor="middle" font-size="12" fill="{t["sub"]}">'
             f'minor allele frequency in the test people</text>')
    o.append(f'<text x="18" y="{(y0 + y1) / 2:.0f}" font-size="12" fill="{t["sub"]}" '
             f'transform="rotate(-90 18 {(y0 + y1) / 2:.0f})" text-anchor="middle">r² with the 30x truth</text>')
    ends = []
    for k, (label, col, vals, dashed) in enumerate(series):
        pts = [(x, ys(v)) for x, v in zip(xs, vals) if v is not None]
        if not pts:
            continue
        d = "M" + " L".join(f"{x:.1f} {y:.1f}" for x, y in pts)
        length = sum(((pts[i + 1][0] - pts[i][0]) ** 2 + (pts[i + 1][1] - pts[i][1]) ** 2) ** 0.5
                     for i in range(len(pts) - 1)) + 2
        dash = ' stroke-dasharray="6 5"' if dashed else ""
        delay = 0.15 * k
        if dashed:
            o.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="2"{dash}>'
                     f'{grow("opacity", 0, 1, delay, 0.8)}</path>')
        else:
            o.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="2.5" stroke-linejoin="round" '
                     f'stroke-dasharray="{length:.0f}" stroke-dashoffset="0">'
                     f'{grow("stroke-dashoffset", f"{length:.0f}", 0, delay, 1.3)}</path>')
        for i, (x, y) in enumerate(pts):
            o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" fill="{col}" stroke="{t["bg"]}" stroke-width="2" '
                     f'>{grow("opacity", 0, 1, delay + 0.2 * i, 0.3)}</circle>')
        ends.append([pts[-1][1], pts[-1][0], label, vals[-1]])
    ends.sort()
    for i in range(1, len(ends)):                       # keep end labels at least 20 px apart
        ends[i][0] = max(ends[i][0], ends[i - 1][0] + 20)
    for ly, lx, label, v in ends:
        o.append(f'<text x="{lx + 14:.1f}" y="{ly + 4:.1f}" font-size="13" font-weight="600" fill="{t["fg"]}">'
                 f'{label} <tspan font-weight="400" fill="{t["sub"]}">{v:.3f}</tspan></text>')
    # legend (identity never by colour alone: the end labels name each line too)
    lx = x0
    for label, col, _, dashed in series:
        dash = ' stroke-dasharray="6 5"' if dashed else ""
        o.append(f'<line x1="{lx}" x2="{lx + 22}" y1="66" y2="66" stroke="{col}" stroke-width="3"{dash}/>')
        o.append(f'<text x="{lx + 28}" y="70" font-size="12" fill="{t["fg"]}">{label}</text>')
        lx += 40 + 7 * len(label)
    o.append("</svg>")
    return "\n".join(o) + "\n"


def bars(t, title, subtitle, rows, unit, w=450, h=210):
    """rows: [(label, colour, value, text)] horizontal bars, one measure per chart."""
    o = frame(t, w, h, title, subtitle)
    vmax = max(v for _, _, v, _ in rows)
    x0, wmax = 100, w - 100 - 110
    for i, (label, col, v, txt) in enumerate(rows):
        y = 78 + 54 * i
        bw = max(wmax * v / vmax, 3)
        o.append(f'<text x="{x0 - 12}" y="{y + 22}" text-anchor="end" font-size="13" font-weight="600" '
                 f'fill="{t["fg"]}">{label}</text>')
        o.append(f'<rect x="{x0}" y="{y}" width="{bw:.1f}" height="34" rx="4" fill="{col}">'
                 f'{grow("width", 0, f"{bw:.1f}", 0.1 + 0.15 * i)}</rect>')
        o.append(f'<text x="{x0 + bw + 10:.1f}" y="{y + 22}" font-size="14" font-weight="700" fill="{t["fg"]}" '
                 f'>{txt}{grow("opacity", 0, 1, 0.9 + 0.15 * i, 0.3)}</text>')
    o.append("</svg>")
    return "\n".join(o) + "\n"


def main():
    s = json.load(open(sys.argv[1]))
    out = sys.argv[2]
    for name, t in THEMES.items():
        def save(fn, text):
            with open(f"{out}/{fn}_{name}.svg", "w") as f:
                f.write(text)
        series = [("Selphi 2", t["selphi"], [s["selphi2"]["bins"].get(b) for b in BINS], False)]
        if "beagle" in s:
            series.append(("Beagle 5.5", t["beagle"], [s["beagle"]["bins"].get(b) for b in BINS], False))
        save("r2", line_chart(t, series, "Imputation accuracy by allele frequency",
                              f'{s["n_test"]} people imputed from {s["n_array"]:,} chip sites on chr22 · '
                              f'imputed sites only · panel of {s["n_panel"]:,} people'))
        if "beagle" in s:
            sp = [("Selphi 2", t["selphi"], s["selphi2"]["sec"], f'{s["selphi2"]["sec"]:.1f} s'),
                  ("Beagle 5.5", t["beagle"], s["beagle"]["sec"], f'{s["beagle"]["sec"]:.1f} s')]
            save("time", bars(t, "Phase + impute: wall time", f'{s["n_test"]} people · chr22 · {s["threads"]} threads', sp, "s"))
            mem = [("Selphi 2", t["selphi"], s["selphi2"]["gb"], f'{s["selphi2"]["gb"]:.1f} GB'),
                   ("Beagle 5.5", t["beagle"], s["beagle"]["gb"], f'{s["beagle"]["gb"]:.1f} GB')]
            save("memory", bars(t, "Phase + impute: peak memory", "resident set size, whole process", mem, "GB"))
        if "lowcov" in s:
            covs = sorted(s["lowcov"], key=float)
            ser = [(f"{c}x reads", t["cov"][i % 3], [s["lowcov"][c]["bins"].get(b) for b in BINS], False)
                   for i, c in enumerate(covs)]
            ser.append(("chip", t["sub"], [s["chip_lowcov_people"]["bins"].get(b) for b in BINS], True))
            save("lowcov", line_chart(t, ser, "Imputing from low-coverage reads",
                                      f'{s["n_lowcov"]} people · chr22 · reads thinned from 30x · '
                                      f'dashed: the same people from the chip'))


if __name__ == "__main__":
    main()
