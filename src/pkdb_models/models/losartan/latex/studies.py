"""Module for creating latex table from study table."""

from pathlib import Path
import pandas as pd
from sbmlutils.console import console


def create_latex_table(df: pd.DataFrame) -> str:
    """Transform DataFrame into LaTeX table."""

    # Preamble:
    latex_header = r"""
\begin{table}[H]
\centering
\tabcolsep=3.5pt\relax
\scriptsize
\begin{threeparttable} 

\caption{\textbf{Summary of studies for modeling.} 
Overview of study identifiers, PK-DB IDs, administred substance and administration route, dosing regimens, doses [mg], 
and subject characteristics, including health status, renal functional impairment (\emph{RFI}), 
hepatic functional impairment (\emph{HFI}), and the studied genotypes (\emph{CYP2C9}, \emph{ABCB1}).}
\label{table:curated_data_overview}

\begin{tabularx}{\textwidth}{
  p{2.4cm}  % Study
  p{1.7cm}  % PK-DB ID
  p{1.5cm}  % Substance
  p{1.0cm}  % Route
  p{1.0cm}  % Dosing
  p{1.0cm}  % Dose [mg]
  c{1.1cm}  % Healthy
  c{0.8cm}  % HRI
  c{0.8cm}  % RFI
  c{1.1cm}  % CYP2C9
  c{1.1cm}  % ABCB1
}
\toprule
""".strip('\n')

    # Rename columns
    df = df.rename(columns={
        "study": "Study",
        "pkdb": "PK-DB",
        "substance": "Substance",
        "route": "Route",
        "dosing": "Dosing",
        "dose": "Dose [mg]",
        "healthy": "Healthy",
        "renal impairment": "RFI",
        "hepatic impairment": "HFI",
    })

    # Convert True/False to checkmarks
    for col in ["Healthy", "RFI", "HFI", "CYP2C9", "ABCB1"]:
        df[col] = df[col].apply(lambda x: r"\checkmark" if str(x).strip().upper() == "TRUE" else "")

    console.print(df)

    # Column headers row
    column_names = df.columns.tolist()
    column_headers = " & ".join([f"\\textbf{{{col}}}" for col in column_names]) + r" \\"

    # Genotypes columns "*" -> "\textasteriskcentered"
    # if "Genotype" in df.columns:
    #     df["Genotype"] = df["Genotype"].apply(
    #         lambda x: x.replace("*", r"\textasteriskcentered") if isinstance(x, str) else x
    #     )

    # LaTeX body
    latex_body = f"{column_headers}\n\\midrule\n"

    for _, row in df.iterrows():
        values = list(row.astype(str).values)
        values[0] = f"{values[0]} \\cite{{{values[0]}}}"

        # PK-DB ID to clickable link
        values[1] = f"\\href{{https://identifiers.org/pkdb:{values[1]}}}{{{values[1]}}}"

        # Converting leftover "TRUE" or "-" strings
        row_str = " & ".join(v.replace("TRUE", r"\checkmark").replace("-", "")
                             for v in values) + r" \\"
        latex_body += row_str + "\n"

    # End of the tabular portion
    latex_footer = r"""
\bottomrule
\end{tabularx}

\end{threeparttable} 
\end{table}
""".strip('\n')

    # Combine all parts
    full_latex = latex_header + "\n" + latex_body + latex_footer
    return full_latex


if __name__ == "__main__":
    # Input and output file paths
    tsv_path = Path(__file__).parent / 'losartan_studies.tsv'
    latex_path = Path(__file__).parent / 'losartan_studies.tex'

    # Load TSV file
    df = pd.read_csv(
        tsv_path, sep="\t", skiprows=1, skipfooter=48,
        usecols=[
            "study",
            "pkdb",
            "substance",
            "route",
            "dosing",
            "dose",
            "healthy",
            "renal impairment",
            "hepatic impairment",
            "CYP2C9",
            "ABCB1",
        ]
    )[[
        "study",
        "pkdb",
        "substance",
        "route",
        "dosing",
        "dose",
        "healthy",
        "renal impairment",
        "hepatic impairment",
        "CYP2C9",
        "ABCB1",
    ]]

    # Replace NaN values with an empty string
    df = df.fillna("")

    # Create LaTeX string
    latex_str = create_latex_table(df)

    # Save LaTeX table to file
    with open(latex_path, 'w') as f_tex:
        f_tex.write(latex_str)
