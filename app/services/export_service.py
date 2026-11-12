"""
Export service for alignment results.
Provides CSV and JSON formatted file downloads and streaming responses using Pandas.
"""

import json
from typing import Any, List, Union
import pandas as pd
from app.models import Alignment


def alignments_to_dataframe(alignments: List[Any]) -> pd.DataFrame:
    """Convert a list of Alignment ORM objects, AlignmentHits, or dictionaries into a Pandas DataFrame."""
    rows = []
    for item in alignments:
        if isinstance(item, dict):
            row = item.copy()
            if isinstance(row.get("mismatch_positions"), list):
                row["mismatch_positions"] = json.dumps(row["mismatch_positions"])
        elif hasattr(item, "to_dict"):
            row = item.to_dict()
            row["id"] = getattr(item, "id", "")
            row["job_id"] = getattr(item, "job_id", "")
            if isinstance(row.get("mismatch_positions"), list):
                row["mismatch_positions"] = json.dumps(row["mismatch_positions"])
        else:
            row = {
                "id": str(getattr(item, "id", "")),
                "job_id": str(getattr(item, "job_id", "")),
                "chromosome": getattr(item, "chromosome", ""),
                "start_position": getattr(item, "start_position", 0),
                "end_position": getattr(item, "end_position", 0),
                "strand": getattr(item, "strand", "+"),
                "matched_sequence": getattr(item, "matched_sequence", ""),
                "pam_sequence": getattr(item, "pam_sequence", ""),
                "mismatch_count": getattr(item, "mismatch_count", 0),
                "mismatch_positions": str(getattr(item, "mismatch_positions", "[]")),
                "seed_mismatch_count": getattr(item, "seed_mismatch_count", 0),
                "score": getattr(item, "score", 0.0),
                "risk_category": getattr(item, "risk_category", ""),
            }
        rows.append(row)

    if not rows:
        return pd.DataFrame(
            columns=[
                "id", "job_id", "chromosome", "start_position", "end_position",
                "strand", "matched_sequence", "pam_sequence", "mismatch_count",
                "mismatch_positions", "seed_mismatch_count", "score", "risk_category"
            ]
        )

    df = pd.DataFrame(rows)
    return df


def export_to_csv(alignments: List[Union[Alignment, dict]]) -> str:
    """Generate CSV string from alignment records."""
    df = alignments_to_dataframe(alignments)
    return df.to_csv(index=False)


def export_to_json(alignments: List[Union[Alignment, dict]]) -> str:
    """Generate pretty-printed JSON string from alignment records."""
    df = alignments_to_dataframe(alignments)
    # Parse back mismatch_positions as JSON lists in the json export
    records = df.to_dict(orient="records")
    for r in records:
        if isinstance(r.get("mismatch_positions"), str):
            try:
                r["mismatch_positions"] = json.loads(r["mismatch_positions"])
            except Exception:
                pass
    return json.dumps(records, indent=2)
