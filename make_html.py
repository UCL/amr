import markdown
import codecs

with open('MODEL_DESCRIPTION.md', 'r', encoding='utf-8') as f:
    text = f.read()

html = markdown.markdown(text, extensions=['tables', 'fenced_code', 'toc', 'md_in_html'])

mathjax_script = '''
<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
'''

with codecs.open('MODEL_DESCRIPTION.html', 'w', encoding='utf-8') as f:
    f.write('<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="UTF-8">\n<style>body { font-family: sans-serif; line-height: 1.6; max-width: 1000px; margin: 0 auto; padding: 20px; } table { border-collapse: collapse; margin-bottom: 20px; } th, td { border: 1px solid #ddd; padding: 8px; text-align: left; } th { background-color: #f2f2f2; } code { background-color: #f4f4f4; padding: 2px 4px; border-radius: 4px; } pre code { display: block; padding: 10px; overflow-x: auto; }</style>\n')
    f.write(mathjax_script)
    f.write('</head>\n<body>\n')
    f.write(html)
    f.write('\n</body>\n</html>')
print('HTML Regenerated.')
