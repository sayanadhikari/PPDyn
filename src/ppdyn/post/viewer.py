"""Interactive, self-contained HTML viewer (plotly) with a time slider."""
from __future__ import annotations

import numpy as np

from .analysis import Run
from .plots import Units, color_values


def viewer_html(run: Run, path, color="speed", max_points=20_000, max_frames=60, si=False, objects=True):
    try:
        import plotly.graph_objects as go
    except ImportError as exc:  # pragma: no cover
        raise ImportError("the interactive viewer needs plotly: pip install 'ppdyn[viz]'") from exc
    u = Units(run, si)
    frames_k = np.unique(np.linspace(0, len(run) - 1, min(len(run), max_frames)).astype(int))
    rng = np.random.default_rng(0)
    idx = np.arange(run.n) if run.n <= max_points else np.sort(rng.choice(run.n, max_points, replace=False))
    allc = [color_values(run, k, color)[0][idx] for k in (frames_k[0], frames_k[-1])]
    cmin, cmax = float(min(c.min() for c in allc)), float(max(c.max() for c in allc))
    clab = color_values(run, 0, color)[1]
    size = 2.5 if len(idx) < 5000 else 1.5

    def trace(k):
        p = run.position(k)[idx] * u.x
        keep = run.state(k)[idx] != 2
        c = color_values(run, k, color)[0][idx]
        return go.Scatter3d(x=p[keep, 0], y=p[keep, 1], z=p[keep, 2], mode="markers",
                            marker=dict(size=size, color=c[keep], colorscale="Viridis", cmin=cmin, cmax=cmax,
                                        colorbar=dict(title=clab), opacity=0.9),
                            hoverinfo="skip", name="particles")

    lo, hi = run.lo * u.x, run.hi * u.x
    frames = [go.Frame(data=[trace(k)], name=str(i)) for i, k in enumerate(frames_k)]
    fig = go.Figure(data=[trace(frames_k[0])], frames=frames)
    steps = [dict(method="animate", label=f"{run.times[k] * u.t:.3g}",
                  args=[[str(i)], dict(mode="immediate", frame=dict(duration=0, redraw=True), transition=dict(duration=0))])
             for i, k in enumerate(frames_k)]
    span = np.maximum(hi - lo, 1e-12)
    fig.update_layout(
        template="plotly_white", height=760, margin=dict(l=0, r=0, t=40, b=0),
        title=f"PPDyn: {run.n} particles ({len(idx)} shown)",
        scene=dict(xaxis=dict(range=[lo[0], hi[0]], title=f"x [{u.xl}]"),
                   yaxis=dict(range=[lo[1], hi[1]], title=f"y [{u.xl}]"),
                   zaxis=dict(range=[lo[2], hi[2]], title=f"z [{u.xl}]"),
                   aspectmode="manual", aspectratio=dict(x=1, y=span[1] / span[0], z=span[2] / span[0])),
        updatemenus=[dict(type="buttons", x=0.02, y=0.02, xanchor="left", yanchor="bottom", direction="left",
                          buttons=[dict(label="Play", method="animate",
                                        args=[None, dict(frame=dict(duration=80, redraw=True), fromcurrent=True)]),
                                   dict(label="Pause", method="animate",
                                        args=[[None], dict(mode="immediate", frame=dict(duration=0, redraw=False))])])],
        sliders=[dict(active=0, x=0.15, len=0.83, y=0.02, currentvalue=dict(prefix=f"t [{u.tl}] = "), steps=steps)],
    )
    fig.write_html(path, include_plotlyjs="cdn", auto_play=False)
    return path
