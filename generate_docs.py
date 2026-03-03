import pdoc
from pathlib import Path
import os

if __name__ == "__main__":

    curr_path = os.path.dirname(os.path.abspath(__file__))

    pdoc.render.configure(
        show_source=False,
        mermaid=True,
        math=True,
        search=True,
        edit_url_map={"IMP_Toolbox": "https://github.com/isblab/IMP_Toolbox/tree/main/IMP_Toolbox/"}
    )
    doc = pdoc.pdoc(
        *[
            os.path.join(curr_path, "IMP_Toolbox"),
        ],
        output_directory=Path(os.path.join(curr_path, "docs")),
    )

    # out = pdoc.render.html_module(
    #     module=doc,
    #     all_modules={
    #         IMP_Toolbox: doc.modules[IMP_Toolbox],
    #     }
    # )

    # with open("pdoc_output.html", "w") as f:
    #     f.write(out)
