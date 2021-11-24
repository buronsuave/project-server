import uuid
import pdflatex as ptex
import base64
import os
from datetime import date

def text_to_pdf(latex):
    # Create .tex file
    id = str(uuid.uuid1())
    filepath = f'compiler/files/{id}.tex'
    today = date.today()
    dateLatex = today.strftime("%B %d, %Y")

    with open(filepath, "w") as file_tex:
        file_tex.write("\\documentclass{article}\n")
        file_tex.write("\\usepackage{amsmath}\n")
        file_tex.write("\\title{ODE Solution}\n")
        file_tex.write("\\author{Solved by PWA for ODE}\n")
        file_tex.write("\\date{"+dateLatex+"}\n")

        file_tex.write("\\begin{document}\n")
        file_tex.write("\\maketitle\n")
        file_tex.write("\\noindent\n")
        print(latex)
        file_tex.write(latex + "\n")
        file_tex.write("\\end{document}\n")
    
    # Create .pdf file
    os.system(f"pdflatex -output-directory=compiler/files {filepath}")

    # Read .pdf file
    encoded_string = None
    with open(f'compiler/files/{id}.pdf', 'rb') as file_pdf:
        encoded_string = str(base64.b64encode(file_pdf.read()))

    # Remove .tex, .pdf, .aux, .log files
    os.remove(filepath)
    os.remove(f'compiler/files/{id}.pdf')
    os.remove(f'compiler/files/{id}.aux')
    os.remove(f'compiler/files/{id}.log')
    encoded_string = encoded_string[2:]
    encoded_string = encoded_string[:-1]

    return encoded_string