import sys
import os
import re
from pathlib import Path

def convert_txt_to_html(input_file: str):
    input_path = Path(input_file)
    if not input_path.exists():
        print(f"Error: {input_file} does not exist.")
        sys.exit(1)
        
    # Create a directory to hold the split table files
    output_dir = input_path.parent / f"{input_path.stem}_tables"
    output_dir.mkdir(exist_ok=True, parents=True)
        
    with open(input_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    def get_html_head(title):
        return [
            "<!DOCTYPE html>",
            "<html lang='en'>",
            "<head>",
            "    <meta charset='UTF-8'>",
            f"    <title>{title}</title>",
            "    <style>",
            "        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; max-width: 1400px; margin: 0 auto; padding: 20px; color: #333; }",
            "        h1 { border-bottom: 2px solid #ddd; padding-bottom: 10px; color: #2c3e50; }",
            "        h2 { margin-top: 30px; color: #2c3e50; border-bottom: 1px solid #eee; padding-bottom: 5px; }",
            "        p { margin-bottom: 10px; }",
            "        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; font-size: 14px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }",
            "        th, td { border: 1px solid #ddd; padding: 8px 12px; text-align: right; }",
            "        th:first-child, td:first-child { text-align: left; }",
            "        th { background-color: #f8f9fa; font-weight: 600; text-align: right!important; }",
            "        th:first-child { text-align: left!important; }",
            "        tr:nth-child(even) td { background-color: #fcfcfc; }",
            "        tr:hover td { background-color: #f1f7fd; }",
            "        .summary-block { background: #f8f9fa; padding: 15px; border-left: 4px solid #3498db; margin-bottom: 20px; border-radius: 0 4px 4px 0; }",
            "        pre { background: #eee; padding: 10px; border-radius: 5px; overflow-x: auto; font-family: Consolas, monospace; }",
            "        .note { font-style: italic; color: #666; font-size: 0.9em; margin-bottom: 10px; }",
            "    </style>",
            "</head>",
            "<body>",
            f"    <h1>{title}</h1>"
        ]

    i = 0
    in_table = False
    table_count = 0
    
    current_content = []
    headers_matches = []
    header_spans = []

    def save_table_file(content_arr, t_num):
        # We find the table title (H2 or H1) if possible to make the filename clear
        title_tag = None
        for block in reversed(content_arr):
            if '<h2>' in block:
                title_tag = block.replace('<h2>', '').replace('</h2>', '').strip()
                break
        
        safe_title = re.sub(r'[^a-zA-Z0-9_\-]', '_', title_tag) if title_tag else ""
        filename = output_dir / f"table_{t_num:02d}_{safe_title[:30] if safe_title else 'data'}.html"
        
        full_html = get_html_head(f"{input_path.stem} | Table {t_num} | {title_tag if title_tag else ''}")
        full_html.extend(content_arr)
        full_html.extend(["</body>", "</html>"])
        with open(filename, 'w', encoding='utf-8') as fout:
            fout.write('\n'.join(full_html))
        print(f"Saved: {filename}")

    while i < len(lines):
        line = lines[i].rstrip('\n')
        
        # Skip completely empty lines if we aren't in a table
        if not line.strip():
            if in_table:
                current_content.append("    </tbody>\n</table>")
                in_table = False
                table_count += 1
                save_table_file(current_content, table_count)
                current_content = []  # Reset for the next table
            i += 1
            continue

        # Look for section headers (lines without double spaces that precede tables)
        if not in_table and '  ' not in line and not line.startswith(' ') and len(line) > 3 and not line.startswith('-') and ' ' in line and ':' not in line[:15] and not line.startswith('Note:'):
            peek = i + 1
            while peek < len(lines) and not lines[peek].strip():
                peek += 1
                
            if peek < len(lines):
                next_line = lines[peek]
                if '  ' in next_line or '\t' in next_line:
                    current_content.append(f"    <h2>{line.strip()}</h2>")
                    i += 1
                    continue

        # Check if line looks like it belongs to a table (multiple spaces between items)
        looks_like_table_row = '  ' in line
        if not in_table and looks_like_table_row and not line.strip().startswith('-'):
            in_table = True
            
            headers_matches = list(re.finditer(r'\S+(?: \S+)*', line))
            headers = [m.group() for m in headers_matches]
            header_spans = [(m.start(), m.end()) for m in headers_matches]
               
            current_content.append("<table>\n    <thead>\n        <tr>")
            for h in headers:
                current_content.append(f"            <th>{h}</th>")
            current_content.append("        </tr>\n    </thead>\n    <tbody>")
            i += 1
            continue
            
        elif in_table and looks_like_table_row:
            row_items = [v for v in re.split(r'\s{2,}', line.strip()) if v]
            
            if len(row_items) != len(headers_matches) and len(headers_matches) > 1:
                splits = []
                for j in range(len(headers_matches)-1):
                    c1 = (header_spans[j][0] + header_spans[j][1]) // 2
                    c2 = (header_spans[j+1][0] + header_spans[j+1][1]) // 2
                    chunk = line[c1:c2]
                    spaces = list(re.finditer(r' +', chunk))
                    if not spaces:
                        splits.append(c1 + len(chunk)//2)
                    else:
                        best = max(spaces, key=lambda m: m.end()-m.start())
                        splits.append(c1 + best.start() + (best.end()-best.start())//2)
                
                row_items = [line[:splits[0]].strip()]
                for j in range(len(splits)-1):
                    row_items.append(line[splits[j]:splits[j+1]].strip())
                row_items.append(line[splits[-1]:].strip())
            
            if any(row_items):
                current_content.append("        <tr>")
                for cell in row_items:
                    current_content.append(f"            <td>{cell}</td>")
                current_content.append("        </tr>")
            else:
                current_content.append("    </tbody>\n</table>")
                in_table = False
                table_count += 1
                save_table_file(current_content, table_count)
                current_content = []
                
            i += 1
            continue
        
        elif in_table and not looks_like_table_row:
            current_content.append("    </tbody>\n</table>")
            in_table = False
            table_count += 1
            save_table_file(current_content, table_count)
            current_content = []

        # Regular text handling outside tables
        if not in_table:
            if line.startswith('Note:'):
                current_content.append(f"    <p class='note'>{line}</p>")
            elif line.startswith('- Mean'):
                 # Wrap 'Summary' bullet points
                 if len(current_content) > 0 and (current_content[-1].startswith('    <p>') or current_content[-1] == "    <div class='summary-block'>"):
                     pass
                 else:
                     current_content.append("    <div class='summary-block'>")
                 current_content.append(f"        <div>{line}</div>")
                 
                 # If next line is not a bullet, close the block
                 if i + 1 < len(lines) and not lines[i+1].startswith('-'):
                     current_content.append("    </div>")
            elif ':' in line:
                parts = line.split(':', 1)
                current_content.append(f"    <p><strong>{parts[0]}:</strong>{parts[1]}</p>")
            else:
                if line.strip():
                    current_content.append(f"    <p>{line}</p>")

        i += 1

    if in_table:
        current_content.append("    </tbody>\n</table>")
        table_count += 1
        save_table_file(current_content, table_count)

    print(f"\nCompleted! Generated {table_count} table files in folder: {output_dir.absolute()}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python txt_to_html_summary.py <path_to_txt_file>")
        sys.exit(1)
        
    in_target = sys.argv[1]
    convert_txt_to_html(in_target)