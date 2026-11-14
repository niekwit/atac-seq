"""Sphinx configuration of the atac-seq workflow documentation"""

project = "atac-seq"
author = "Niek Wit"
copyright = "2026, Niek Wit"

extensions = [
    "myst_parser",
    "sphinx_copybutton",
]

# Markdown (MyST) pages
source_suffix = {".md": "markdown"}
myst_enable_extensions = ["colon_fence", "deflist"]
# Anchors for headings, so pages can link to sections (page.md#section)
myst_heading_anchors = 3

# Scripts and configs of the TE alignment comparison are not pages
exclude_patterns = ["_build", "te_alignment_comparison/config"]

html_theme = "furo"
html_title = "atac-seq"
html_theme_options = {
    "source_repository": "https://github.com/niekwit/atac-seq",
    "source_branch": "main",
    "source_directory": "docs/",
}

# Copy only the commands, not the prompt or output, from code blocks
copybutton_prompt_text = r"\$ "
copybutton_prompt_is_regexp = True
