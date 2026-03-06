#!/usr/bin/env python3
"""
Submit a job to https://comet.lih.lu/index.php and capture the results
in the same Python script.

Supports:
- paste FASTA into the form
- upload a FASTA file

Tested as a generic form-submitter by discovering the form structure at runtime,
so it does not hardcode fragile field names.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup

import pandas as pd


BASE_URL = "https://comet.lih.lu/index.php?cat=hiv1"  # switch cat=hiv2 or cat=hcv if needed


def pick_form(soup: BeautifulSoup, mode: str):
    """
    Heuristically pick the correct form:
    - paste mode: form containing a <textarea>
    - upload mode: form containing <input type=file>
    """
    forms = soup.find_all("form")
    if not forms:
        raise RuntimeError("No forms found on page.")

    if mode == "paste":
        for form in forms:
            if form.find("textarea") is not None:
                return form
    elif mode == "upload":
        for form in forms:
            if form.find("input", {"type": "file"}) is not None:
                return form

    raise RuntimeError(f"Could not find a suitable form for mode={mode!r}.")


def build_payload(form, fasta_text: str | None):
    """
    Collect hidden/default inputs and populate the textarea/checkbox fields.
    """
    data: dict[str, str] = {}
    files = None

    # Copy hidden inputs and preserve default values where present
    for inp in form.find_all("input"):
        name = inp.get("name")
        if not name:
            continue

        typ = (inp.get("type") or "text").lower()
        value = inp.get("value", "")

        if typ in {"hidden", "submit"}:
            data[name] = value

    # Find and tick required confirmation checkbox(es)
    for cb in form.find_all("input", {"type": "checkbox"}):
        name = cb.get("name")
        if name:
            # Most PHP forms accept the checkbox value if checked; default to "on" if absent
            data[name] = cb.get("value", "on")

    # Fill textarea for paste mode
    if fasta_text is not None:
        textarea = form.find("textarea")
        if textarea is None:
            raise RuntimeError("Paste form has no textarea.")
        textarea_name = textarea.get("name")
        if not textarea_name:
            raise RuntimeError("Textarea has no name attribute.")
        data[textarea_name] = fasta_text

    return data, files


def submit_paste(session: requests.Session, page_url: str, fasta_text: str, timeout: int = 120):
    r = session.get(page_url, timeout=timeout)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "html.parser")
    form = pick_form(soup, "paste")

    method = (form.get("method") or "post").lower()
    action = urljoin(page_url, form.get("action") or page_url)
    data, _ = build_payload(form, fasta_text)

    if method == "get":
        resp = session.get(action, params=data, timeout=timeout)
    else:
        resp = session.post(action, data=data, timeout=timeout)

    resp.raise_for_status()
    return resp


def submit_upload(session: requests.Session, page_url: str, fasta_path: Path, timeout: int = 120):
    r = session.get(page_url, timeout=timeout)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "html.parser")
    form = pick_form(soup, "upload")

    method = (form.get("method") or "post").lower()
    action = urljoin(page_url, form.get("action") or page_url)

    data, _ = build_payload(form, fasta_text=None)

    file_input = form.find("input", {"type": "file"})
    if file_input is None or not file_input.get("name"):
        raise RuntimeError("Upload form has no usable file input.")

    file_field = file_input["name"]

    with fasta_path.open("rb") as fh:
        files = {
            file_field: (fasta_path.name, fh, "application/octet-stream"),
        }
        if method == "get":
            raise RuntimeError("Upload form unexpectedly uses GET.")
        resp = session.post(action, data=data, files=files, timeout=timeout)

    resp.raise_for_status()
    return resp


def extract_results(html: str) -> str:
    """
    Best-effort extraction:
    - keep tables
    - keep <pre>
    - otherwise return visible text
    """
    soup = BeautifulSoup(html, "html.parser")

    chunks = []

    # Tables often hold COMET output
    for table in soup.find_all("table"):
        rows = []
        for tr in table.find_all("tr"):
            cells = [td.get_text(" ", strip=True) for td in tr.find_all(["th", "td"])]
            if cells:
                rows.append("\t".join(cells))
        if rows:
            chunks.append("\n".join(rows))

    # Preformatted output
    for pre in soup.find_all("pre"):
        txt = pre.get_text("\n", strip=True)
        if txt:
            chunks.append(txt)

    if chunks:
        return "\n\n".join(chunks)

    # Fallback: visible body text
    body = soup.get_text("\n", strip=True)
    return body


def format_hivtyper(text: str, infilename: str):
    name1 = infilename.rsplit("/")[-1] # gives a file name.fasta
    name2 = name1.split("_")[1] # gives a middle part after splitting by "_"
    name3 = name1.rsplit(".")[-2] # gives a file name (cuts .fasta)

    # skip first lines containing the version change info
    # write the rest to a .csv file (it is separated by tab, but we will read it with pandas and write it again with comma separation)
    with open("comet_" + name3 + ".csv", "w") as f:
        dataline = False
        for line in text.splitlines():
            if line == "name	virus	subtype	support (1)":
                dataline = True
            if dataline:
                f.write(line + "\n")

    # Read .csv (it is separated by tab)
    df = pd.read_csv("comet_" + name3 + ".csv", sep="\t")

    # Rename some columns (as done for stanford df)
    df.rename(columns = {"name":"SequenceName", "subtype": "Comet_" + name2 + "_Subtype"}, inplace = True)

    # Add to the "Comment" column bootstrap support info
    df["Comet_" + name2 + "_Comment"] = df["support (1)"].astype(str)

    # Delete undesired columns
    df.drop(columns=["virus", "support (1)"], axis = 1,  inplace = True)

    # Replace some patterns so they look Stanford-like (add CRF)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"^(\w{2})_(\D)(\d?)(\w)(\w{0,1}?)(\d?)$", r"CRF\1_\2\4\5", regex=True)

    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"^(\d{2,3})_(\w{2,4})$", r"CRF\1_\2", regex=True)

    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"^(\w{1,2})(\s\(check for\s)(\d{2,4}_\w{2,4}\))$", r"\1\2CRF\3", regex=True)

    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"^(\d{2,3}_\w{2,4})(\s\(check for\s)(\d{2,4}_\w{2,4}\))$", r"CRF\1\2CRF\3", regex=True)

    # Replace according to LANL
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF03_AB", r"CRF03_A6B", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF12_BF", r"CRF12_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF16_AD", r"CRF16_A2D", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF29_BF$", r"CRF29_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF35_AD", r"CRF35_A1D", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF38_BF$", r"CRF38_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF39_BF$", r"CRF39_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF40_BF$", r"CRF40_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF44_BF$", r"CRF44_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF47_BF$", r"CRF47_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF50_AD", r"CRF50_A1D", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF70_BF$", r"CRF70_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF71_BF$", r"CRF71_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF72_BF$", r"CRF72_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF89_BF$", r"CRF89_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF90_BF$", r"CRF90_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF99_BF$", r"CRF99_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF122_BF$", r"CRF122_BF1", regex=True)
    df["Comet_" + name2 + "_Subtype"] = df["Comet_" + name2 + "_Subtype"].replace(r"CRF141_BF$", r"CRF141_BF1", regex=True)

    # Replace "unassigned_" group with "Unassigned"
    df.loc[df["Comet_" + name2 + "_Subtype"].str.contains("unassigned"), "Comet_" + name2 + "_Subtype"] = "Unassigned"

    # Replace "nan" with "0"
    df["Comet_" + name2 + "_Comment"] = df["Comet_" + name2 + "_Comment"].replace("nan", "0")

    # Sort df by SequenceName
    df = df.sort_values(by=["SequenceName"])

    # Prepare a clean .csv file
    df.to_csv("comet_" + name3 + ".csv", sep=",", index=False, encoding="utf-8")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default=BASE_URL, help="COMET URL")
    ap.add_argument("--fasta-file", type=Path, required=True, help="Upload this FASTA file")
    ap.add_argument("--save-html", type=Path, default=Path("comet_result.html"))
    ap.add_argument("--save-text", type=Path, default=Path("comet_result.txt"))
    args = ap.parse_args()

    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) Python requests",
        }
    )

    try:
        resp = submit_upload(session, args.url, args.fasta_file)

        args.save_html.write_text(resp.text, encoding="utf-8")
        parsed = extract_results(resp.text)
        args.save_text.write_text(parsed, encoding="utf-8")
        format_hivtyper(parsed, args.fasta_file.name)

        print(f"Saved raw HTML to: {args.save_html}")
        print(f"Saved parsed text to: {args.save_text}")
        print("\n=== Parsed result preview ===\n")
        print(parsed[:4000])

    except requests.HTTPError as e:
        print(f"HTTP error: {e}", file=sys.stderr)
        if e.response is not None:
            print(e.response.text[:2000], file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()