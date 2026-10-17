# Sphinx configuration for PPDyn (Markdown sources via MyST)
import os
import sys

sys.path.insert(0, os.path.abspath("../../src"))

project = "PPDyn"
author = "Sayan Adhikari"
copyright = "2021-2026, Sayan Adhikari and PPDyn contributors"
try:
    from ppdyn import __version__ as release
except Exception:  # docs build without the compiled deps
    release = "2.0.0"
version = release

extensions = ["myst_parser"]
myst_enable_extensions = ["dollarmath", "colon_fence"]
source_suffix = {".md": "markdown"}
master_doc = "index"
html_theme = "furo"
html_static_path = ["_static"]
html_logo = "_static/ppdyn_logo.png"
html_title = f"PPDyn {release}"
