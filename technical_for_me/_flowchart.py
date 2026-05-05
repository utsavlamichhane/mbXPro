"""
Flowchart helpers used by the technical_for_me docx builders.

Designed for clear black-on-white diagrams that print well on A4 paper.
Each shape kind matches a standard flowchart symbol.
"""

from __future__ import annotations

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MplPath
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


# ── Color palette (match the public docs) ─────────────────────────────────────
NAVY     = "#2C3E50"
BLUE     = "#1F6FEB"
GREEN    = "#1F8A4F"
ORANGE   = "#D35400"
RED      = "#C0392B"
GREY     = "#5A6470"
LIGHTBLUE = "#E7F1FF"
LIGHTGREY = "#F2F4F7"
LIGHTGREEN = "#E5F4EC"
LIGHTRED   = "#FBE6E3"
LIGHTORANGE = "#FCEEDC"


def _box(ax, x, y, w, h, label, kind="proc", color=None, fontsize=9):
    """Draw a single labeled flowchart box and return its center coords."""
    if kind == "proc":
        face = LIGHTBLUE if color is None else color
        edge = BLUE
        rounding = "round,pad=0.05,rounding_size=0.04"
    elif kind == "io":
        face = LIGHTGREY if color is None else color
        edge = GREY
        rounding = "round,pad=0.05,rounding_size=0.10"
    elif kind == "decision":
        face = LIGHTORANGE if color is None else color
        edge = ORANGE
        rounding = "round,pad=0.05,rounding_size=0.20"
    elif kind == "term":
        face = LIGHTGREEN if color is None else color
        edge = GREEN
        rounding = "round,pad=0.05,rounding_size=0.50"
    elif kind == "fail":
        face = LIGHTRED if color is None else color
        edge = RED
        rounding = "round,pad=0.05,rounding_size=0.50"
    else:
        face = LIGHTBLUE
        edge = NAVY
        rounding = "round,pad=0.05,rounding_size=0.04"

    patch = FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                           boxstyle=rounding,
                           linewidth=1.2, edgecolor=edge, facecolor=face)
    ax.add_patch(patch)
    ax.text(x, y, label, ha="center", va="center",
            fontsize=fontsize, color=NAVY,
            fontfamily="DejaVu Sans", wrap=True)
    return (x, y)


def _arrow(ax, p_from, p_to, label=None, dx=0.0, dy=0.0):
    """Draw an arrow between two box centers."""
    arrow = FancyArrowPatch(p_from, p_to,
                            arrowstyle="-|>",
                            mutation_scale=14,
                            linewidth=1.0,
                            color=NAVY,
                            shrinkA=14, shrinkB=14)
    ax.add_patch(arrow)
    if label:
        midx = (p_from[0] + p_to[0]) / 2 + dx
        midy = (p_from[1] + p_to[1]) / 2 + dy
        ax.text(midx, midy, label, fontsize=8, color=GREY,
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.9))


def render_flowchart(nodes, edges, out_path,
                     title=None, figsize=(8.0, 9.0)):
    """
    Render a flowchart and save to PNG.

    nodes : list of dicts with keys
              id, x, y, w, h, label, kind  ('proc'|'io'|'decision'|'term'|'fail')
              optional: fontsize, color
    edges : list of dicts with keys
              from, to,  optional: label, dx, dy
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=160)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.set_aspect("equal")
    ax.axis("off")

    if title:
        ax.text(5.0, 11.6, title, fontsize=12, fontweight="bold",
                ha="center", va="top", color=NAVY)

    centers = {}
    for n in nodes:
        c = _box(ax, n["x"], n["y"], n["w"], n["h"], n["label"],
                 kind=n.get("kind", "proc"),
                 color=n.get("color"),
                 fontsize=n.get("fontsize", 9))
        centers[n["id"]] = c

    for e in edges:
        _arrow(ax, centers[e["from"]], centers[e["to"]],
               label=e.get("label"),
               dx=e.get("dx", 0.0),
               dy=e.get("dy", 0.0))

    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    return out_path


def render_legend(out_path):
    """Generate the symbol-legend image once, embedded by every doc."""
    fig, ax = plt.subplots(figsize=(7.5, 1.6), dpi=160)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 2)
    ax.axis("off")

    items = [
        (1.0, 1.0, 1.4, 0.6, "Process",       "proc"),
        (3.0, 1.0, 1.4, 0.6, "Input/output",  "io"),
        (5.0, 1.0, 1.4, 0.6, "Decision",      "decision"),
        (7.0, 1.0, 1.4, 0.6, "Start/end",     "term"),
        (9.0, 1.0, 1.4, 0.6, "Failure",       "fail"),
    ]
    for x, y, w, h, label, kind in items:
        _box(ax, x, y, w, h, label, kind=kind, fontsize=9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path
