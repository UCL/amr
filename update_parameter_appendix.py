from __future__ import annotations

import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
MODEL_DESCRIPTION_PATH = REPO_ROOT / "MODEL_DESCRIPTION.md"
APPENDIX_B_HEADING = "## Appendix B — Parameter Reference"
APPENDIX_C_HEADING = "## Appendix C — Output Specification"


def generate_appendix_markdown() -> str:
    result = subprocess.run(
        ["cargo", "run", "--quiet", "--bin", "dump_parameter_appendix"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() + "\n\n"


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
    appendix_markdown = generate_appendix_markdown()
    updated_text = replace_appendix_block(existing_text, appendix_markdown)
    MODEL_DESCRIPTION_PATH.write_text(updated_text, encoding="utf-8")


if __name__ == "__main__":
    main()