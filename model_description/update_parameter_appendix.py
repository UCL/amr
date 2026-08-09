from __future__ import annotations

import subprocess
from pathlib import Path


DOC_DIR = Path(__file__).resolve().parent
REPO_ROOT = DOC_DIR.parent
MODEL_DESCRIPTION_PATH = DOC_DIR / "MODEL_DESCRIPTION.md"
APPENDIX_B_HEADING = "## Appendix B — Parameter Reference"
APPENDIX_C_HEADING = "## Appendix C — Output Specification"


def generate_appendix_markdown(newline: str) -> str:
    result = subprocess.run(
        ["cargo", "run", "--quiet", "--bin", "dump_parameter_appendix"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return (result.stdout.strip() + "\n\n").replace("\n", newline)


def replace_appendix_block(document_text: str, appendix_markdown: str) -> str:
    appendix_b_start = document_text.index(APPENDIX_B_HEADING)
    appendix_c_start = document_text.index(APPENDIX_C_HEADING)
    return (
        document_text[:appendix_b_start]
        + appendix_markdown
        + document_text[appendix_c_start:]
    )


def main() -> None:
    existing_text = MODEL_DESCRIPTION_PATH.read_text(encoding="utf-8")
    newline = "\r\n" if "\r\n" in existing_text else "\n"
    appendix_markdown = generate_appendix_markdown(newline)
    updated_text = replace_appendix_block(existing_text, appendix_markdown)
    with MODEL_DESCRIPTION_PATH.open("w", encoding="utf-8", newline="") as file:
        file.write(updated_text)


if __name__ == "__main__":
    main()
