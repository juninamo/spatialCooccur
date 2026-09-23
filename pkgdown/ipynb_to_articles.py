#!/usr/bin/env python3
"""Convert the executed Jupyter tutorials in vignettes/ into pkgdown articles.

The notebooks stay the single source of truth. Each vignettes/<name>.ipynb
becomes vignettes/articles/<name>.Rmd, a *static* article: code is shown in
plain ```r blocks (knitr does not run it) and the stored outputs and figures
are embedded. This keeps heavy tutorials (e.g. the 10x Xenium one, which
needs large external data) on the website and in its search index without
re-running them during the site build.

Usage (from the package root):
    python3 pkgdown/ipynb_to_articles.py

Re-run after re-executing a notebook, then pkgdown::build_site().
"""
import base64
import json
import pathlib
import re
import shutil

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "vignettes"
OUT = SRC / "articles"
REPO = "https://github.com/juninamo/spatialCooccur/blob/master/vignettes"

# stderr lines that are package-loading noise rather than tutorial content
NOISE = [
    r"^Warning messages?:?$",
    r"replacing previous import",
    r"^\d+: replacing",
    r"^Loading spatialCooccur",
    r"^To cite this package",
    r"^\s+Inamo J, et al\.",
    r"^\s*(juvenile idiopathic|arthritis|immune-stromal|doi:)",
    r"^Developed by:",
    r"^Attaching SeuratObject",
    r"^Loading required package",
    r"^The legacy packages maptools",
    r"^Attaching package",
    r"^The following objects? (is|are) masked",
    r"^\s{4}\S",
    r"^\s*$",
    r"was built under R version",
    r"^SCP python environment",
    r"^If you have already created an SCP",
    r"^\u201cNot validating \w+ objects\u201d$",
    r"^Found more than one class",
    r"^Also defined by",
    r"^\*+$",
    r"^\|$",
    r"^0%\s+10\s+20",
    r"^\[-+\|",
]
NOISE_RE = [re.compile(p) for p in NOISE]
PROGRESS_RE = re.compile(r"^\s*\|[=\s]*\|\s*\d+%$")


def text(x):
    return "".join(x) if isinstance(x, list) else x


def fence(body, info=""):
    body = body.rstrip("\n")
    ticks = "````" if "```" in body else "```"
    return f"{ticks}{info}\n{body}\n{ticks}\n"


ANSI_RE = re.compile(r"\x1b\[[0-9;]*m|\[[0-9;]*m(?=\S)")


def strip_ansi(s):
    return ANSI_RE.sub("", s)


def clean_stderr(s):
    lines = strip_ansi(s).splitlines()
    # drop package start-up banners framed by "=====" rulers
    out, in_banner = [], False
    for l in lines:
        if re.match(r"^=+$", l.strip()):
            in_banner = not in_banner
            continue
        if not in_banner:
            out.append(l)
    keep = []
    for l in out:
        if any(r.search(l) for r in NOISE_RE):
            continue
        if keep and keep[-1] == l:   # collapse repeated warnings
            continue
        keep.append(l)
    return "\n".join(keep)


def guard_inline_r(s, where):
    # knitr would evaluate `r ...` inline code even outside chunks.
    return re.sub(r"`r ", "`&#8203;r ", s)


def convert(nb_path):
    name = nb_path.stem
    nb = json.loads(nb_path.read_text())
    fig_dir = OUT / "figures" / name
    if fig_dir.exists():
        shutil.rmtree(fig_dir)
    fig_dir.mkdir(parents=True)

    title, parts, n_fig = None, [], 0
    for cell in nb["cells"]:
        src = text(cell["source"]).strip("\n")
        if cell["cell_type"] == "markdown":
            if title is None:
                m = re.match(r"^#\s+(.+?)\s*$", src.splitlines()[0]) if src else None
                if m:
                    title = m.group(1)
                    src = "\n".join(src.splitlines()[1:]).strip("\n")
            if src:
                parts.append(guard_inline_r(src, name) + "\n")
            continue
        if cell["cell_type"] != "code" or not src.strip():
            continue
        parts.append(fence(src, "r"))
        for out in cell.get("outputs", []):
            kind = out["output_type"]
            if kind == "stream":
                body = text(out["text"])
                if out.get("name") == "stderr":
                    body = clean_stderr(body)
                else:
                    body = "\n".join(l for l in strip_ansi(body).splitlines()
                                      if not PROGRESS_RE.match(l))
                if body.strip():
                    parts.append(fence(guard_inline_r(body, name), "{.output}"))
            elif kind in ("display_data", "execute_result"):
                data = out.get("data", {})
                if "image/png" in data:
                    n_fig += 1
                    fn = f"fig-{n_fig:02d}.png"
                    (fig_dir / fn).write_bytes(base64.b64decode(text(data["image/png"])))
                    parts.append(f'![](figures/{name}/{fn}){{width="100%"}}\n')
                elif "text/html" in data:
                    parts.append('<div class="nb-output">\n' + text(data["text/html"]).strip() + "\n</div>\n")
                elif "text/plain" in data:
                    parts.append(fence(guard_inline_r(text(data["text/plain"]), name), "{.output}"))
            elif kind == "error":
                parts.append(fence("\n".join(out.get("traceback", [])), "{.output}"))

    # pkgdown uses the YAML title as the page <h1> and indexes <h2>/<h3>
    # sections for search, so demote headings if the body still has an <h1>.
    md_idx = [i for i, p in enumerate(parts) if not p.startswith(("```", "![", "<div"))]
    if any(re.search(r"^# ", parts[i], re.M) for i in md_idx):
        for i in md_idx:
            parts[i] = re.sub(r"^(#{1,5}) ", r"#\1 ", parts[i], flags=re.M)

    title = (title or name).replace('"', '\\"')
    header = (
        "---\n"
        f'title: "{title}"\n'
        "---\n\n"
        '<div class="nb-source">\n'
        f"This article is a rendered copy of the Jupyter notebook "
        f"[`vignettes/{nb_path.name}`]({REPO}/{nb_path.name}); "
        "download it to run the code yourself.\n"
        "</div>\n\n"
    )
    (OUT / f"{name}.Rmd").write_text(header + "\n".join(parts))
    print(f"{nb_path.name} -> articles/{name}.Rmd ({n_fig} figures)")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for nb in sorted(SRC.glob("*.ipynb")):
        convert(nb)
