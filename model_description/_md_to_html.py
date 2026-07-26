"""Convert MODEL_DESCRIPTION.md to MODEL_DESCRIPTION.html with styling."""
import markdown
import re

with open("MODEL_DESCRIPTION.md", encoding="utf-8") as f:
    md_text = f.read()

# Convert markdown to HTML body
body = markdown.markdown(
    md_text,
    extensions=["tables", "fenced_code", "toc"],
    extension_configs={"toc": {"permalink": False}},
)

# Wrap KaTeX-style $...$ and $$...$$ for MathJax rendering
# $$...$$ blocks → \[...\]
body = re.sub(r'\$\$(.+?)\$\$', r'\\[\1\\]', body, flags=re.DOTALL)
# inline $...$ → \(...\)  (but not $$)
body = re.sub(r'(?<!\$)\$(?!\$)(.+?)(?<!\$)\$(?!\$)', r'\\(\1\\)', body)

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>AMR Simulation — Model Description</title>
<style>
body {{
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  line-height: 1.6;
  max-width: 1000px;
  margin: 0 auto;
  padding: 20px;
  color: #333;
}}
h1 {{ border-bottom: 2px solid #2c3e50; padding-bottom: 10px; }}
h2 {{ border-bottom: 1px solid #ddd; padding-bottom: 6px; margin-top: 40px; }}
h3 {{ margin-top: 30px; }}
table {{
  border-collapse: collapse;
  margin-bottom: 20px;
  width: 100%;
  display: block;
  overflow-x: auto;
}}
th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: left; }}
th {{ background-color: #f2f2f2; font-weight: 600; }}
tr:nth-child(even) {{ background-color: #fafafa; }}
code {{ background-color: #f4f4f4; padding: 2px 6px; border-radius: 4px; font-size: 0.9em; }}
pre code {{ display: block; padding: 12px; overflow-x: auto; background-color: #f8f8f8; border: 1px solid #e0e0e0; border-radius: 6px; }}
blockquote {{ border-left: 4px solid #2c3e50; margin: 1em 0; padding: 0.5em 1em; background-color: #f9f9f9; }}
strong {{ color: #2c3e50; }}
hr {{ border: none; border-top: 1px solid #ddd; margin: 30px 0; }}
a {{ color: #2980b9; text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
</style>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
</head>
<body>
{body}
</body>
</html>
"""

with open("MODEL_DESCRIPTION.html", "w", encoding="utf-8") as f:
    f.write(html)

print("MODEL_DESCRIPTION.html generated successfully.")
