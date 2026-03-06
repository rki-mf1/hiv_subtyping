#!/usr/bin/env python3
from __future__ import annotations

import argparse
import io
import re
import sys
import time
from pathlib import Path
from urllib.parse import urljoin, urlparse, parse_qs

import pandas as pd
import requests
from bs4 import BeautifulSoup


DEFAULT_COMET_URL = "https://comet.lih.lu/index.php?cat=hiv1"


class CometError(RuntimeError):
    pass


def find_upload_form(html: str, base_url: str) -> tuple[str, str, dict[str, str], str]:
    """
    Discover the upload form dynamically.

    Returns:
        method, action_url, form_data, file_field_name
    """
    soup = BeautifulSoup(html, "html.parser")
    forms = soup.find_all("form")
    if not forms:
        raise CometError("No forms found on COMET page.")

    for form in forms:
        file_input = form.find("input", {"type": "file"})
        if not file_input or not file_input.get("name"):
            continue

        method = (form.get("method") or "post").lower()
        action_url = urljoin(base_url, form.get("action") or base_url)

        form_data: dict[str, str] = {}

        for inp in form.find_all("input"):
            name = inp.get("name")
            if not name:
                continue

            input_type = (inp.get("type") or "text").lower()
            value = inp.get("value", "")

            if input_type == "hidden":
                form_data[name] = value
            elif input_type == "checkbox":
                # Tick all checkboxes in the upload form, including the required
                # non-commercial-use confirmation.
                form_data[name] = value or "on"
            elif input_type == "submit" and value:
                # Keep one submit value if present.
                form_data.setdefault(name, value)

        return method, action_url, form_data, file_input["name"]

    raise CometError("Could not find an upload form with a file input.")


def extract_job_id(text: str, response_url: str | None = None) -> str:
    """
    Extract job id from HTML or URL.
    """
    candidates = []

    if response_url:
        parsed = urlparse(response_url)
        qs = parse_qs(parsed.query)
        if "job" in qs and qs["job"]:
            candidates.append(qs["job"][0])

    patterns = [
        r"[?&]job=([A-Za-z0-9_-]+)",
        r'csv\.php\?job=([A-Za-z0-9_-]+)',
        r'"job"\s*:\s*"([A-Za-z0-9_-]+)"',
        r"'job'\s*:\s*'([A-Za-z0-9_-]+)'",
    ]
    for pat in patterns:
        m = re.search(pat, text)
        if m:
            candidates.append(m.group(1))

    for job_id in candidates:
        if job_id:
            return job_id

    raise CometError("Could not extract COMET job ID from submission response.")


def looks_like_csv(text: str) -> bool:
    """
    Heuristic: ready CSV should contain at least one non-empty line with commas/semicolons/tabs.
    """
    stripped = text.strip()
    if not stripped:
        return False

    lines = [line for line in stripped.splitlines() if line.strip()]
    if not lines:
        return False

    head = "\n".join(lines[:5]).lower()
    if "error" in head and "job" in head:
        return False
    if "not found" in head or "invalid" in head:
        return False
    if "," in lines[0] or ";" in lines[0] or "\t" in lines[0]:
        return True

    # Sometimes CSV-like output may only show delimiter after the first line
    return any(("," in line or ";" in line or "\t" in line) for line in lines[:5])


def fetch_csv_when_ready(
    session: requests.Session,
    csv_url: str,
    timeout_seconds: int = 300,
    poll_interval_seconds: float = 3.0,
) -> str:
    """
    Poll until CSV becomes available or timeout is reached.
    """
    deadline = time.time() + timeout_seconds
    last_status = None
    last_body = ""

    while time.time() < deadline:
        try:
            resp = session.get(csv_url, timeout=60)
            last_status = resp.status_code
            last_body = resp.text

            if resp.ok and looks_like_csv(resp.text):
                return resp.text

        except requests.RequestException as e:
            last_body = f"{type(e).__name__}: {e}"

        time.sleep(poll_interval_seconds)

    raise CometError(
        f"Timed out waiting for CSV.\n"
        f"URL: {csv_url}\n"
        f"Last HTTP status: {last_status}\n"
        f"Last response excerpt: {last_body[:500]!r}"
    )


def read_csv_to_dataframe(csv_text: str) -> pd.DataFrame:
    """
    Load CSV robustly, trying common delimiters.
    """
    for sep in [",", ";", "\t"]:
        try:
            df = pd.read_csv(io.StringIO(csv_text), sep=sep)
            if df.shape[1] >= 2:
                return df
        except Exception:
            pass

    # Final fallback: let pandas sniff with python engine
    try:
        return pd.read_csv(io.StringIO(csv_text), sep=None, engine="python")
    except Exception as e:
        raise CometError(f"Downloaded content does not parse as CSV: {e}") from e


def submit_one_multifasta(
    fasta_path: Path,
    comet_url: str = DEFAULT_COMET_URL,
    wait_timeout: int = 300,
    poll_interval: float = 3.0,
) -> tuple[str, Path, pd.DataFrame]:
    """
    Submit exactly one multi-FASTA file, wait for CSV, save it, and return DataFrame.
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"Input file not found: {fasta_path}")
    if not fasta_path.is_file():
        raise CometError(f"Input path is not a file: {fasta_path}")

    out_csv = fasta_path.with_suffix(fasta_path.suffix + ".comet.csv")

    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) Python requests",
        }
    )

    landing = session.get(comet_url, timeout=60)
    landing.raise_for_status()

    method, action_url, form_data, file_field_name = find_upload_form(
        landing.text, comet_url
    )

    if method != "post":
        raise CometError(f"Unexpected upload form method: {method!r}")

    with fasta_path.open("rb") as fh:
        files = {
            file_field_name: (
                fasta_path.name,
                fh,
                "application/octet-stream",
            )
        }
        submit_resp = session.post(action_url, data=form_data, files=files, timeout=120)

    submit_resp.raise_for_status()

    job_id = extract_job_id(submit_resp.text, submit_resp.url)
    csv_url = f"https://comet.lih.lu/csv.php?job={job_id}"

    csv_text = fetch_csv_when_ready(
        session=session,
        csv_url=csv_url,
        timeout_seconds=wait_timeout,
        poll_interval_seconds=poll_interval,
    )

    out_csv.write_text(csv_text, encoding="utf-8")
    df = read_csv_to_dataframe(csv_text)   

    return job_id, out_csv, df


def format_hivtyper(df: pd.DataFrame, infilename: str):
    name1 = infilename.rsplit("/")[-1] # gives a file name.fasta
    name2 = name1.split("_")[1] # gives a middle part after splitting by "_"
    name3 = name1.rsplit(".")[-2] # gives a file name (cuts .fasta)

    # Rename some columns (as done for stanford df)
    df.rename(columns = {"name":"SequenceName", "subtype": "Comet_" + name2 + "_Subtype"}, inplace = True)

    # Add to the "Comment" column bootstrap support info
    df["Comet_" + name2 + "_Comment"] = df["bootstrap support"].astype(str)

    # Delete undesired columns
    df.drop(columns=["virus", "bootstrap support"], axis = 1,  inplace = True)

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

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Submit exactly one multi-FASTA file to COMET, wait for CSV, save it, and create a pandas DataFrame."
    )
    parser.add_argument(
        "multifasta",
        type=Path,
        help="Path to one multi-FASTA file (.fa/.fasta/.gz also works if COMET accepts it).",
    )
    parser.add_argument(
        "--comet-url",
        default=DEFAULT_COMET_URL,
        help="COMET endpoint, e.g. https://comet.lih.lu/index.php?cat=hiv1",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=300,
        help="Maximum seconds to wait for the CSV.",
    )
    parser.add_argument(
        "--poll-interval",
        type=float,
        default=5.0,
        help="Polling interval in seconds.",
    )

    args = parser.parse_args()

    try:
        job_id, csv_path, df = submit_one_multifasta(
            fasta_path=args.multifasta,
            comet_url=args.comet_url,
            wait_timeout=args.timeout,
            poll_interval=args.poll_interval,
        )
        format_hivtyper(df, args.multifasta.name)

        print(f"COMET job_id: {job_id}")
        print(f"CSV saved to: {csv_path}")
        print(f"DataFrame shape: {df.shape}")
        print(df.head().to_string(index=False))
        return 0

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())