"""
Unit and API integration tests for CSV and JSON result export.
"""

import json
from app.services.export_service import export_to_csv, export_to_json
from app.core.aligner_cpu import AlignmentHit


def test_export_to_csv():
    hits = [
        AlignmentHit(
            chromosome="chr1",
            start_position=100,
            end_position=123,
            strand="+",
            matched_sequence="GAGTCCGAGCAGAAGAAGAA",
            pam_sequence="TGG",
            mismatch_count=0,
            mismatch_positions=[],
            seed_mismatch_count=0,
            score=100.0,
            risk_category="Exact match",
        ),
        AlignmentHit(
            chromosome="chr2",
            start_position=500,
            end_position=523,
            strand="-",
            matched_sequence="GCGTCCGAGCAGAAGAAGAA",
            pam_sequence="CGG",
            mismatch_count=1,
            mismatch_positions=[2],
            seed_mismatch_count=0,
            score=90.0,
            risk_category="High-risk off-target",
        ),
    ]

    csv_output = export_to_csv(hits)
    assert "chromosome,start_position,end_position" in csv_output
    assert "chr1,100,123,+,GAGTCCGAGCAGAAGAAGAA,TGG,0" in csv_output
    assert "chr2,500,523,-,GCGTCCGAGCAGAAGAAGAA,CGG,1" in csv_output


def test_export_to_json():
    hits = [
        AlignmentHit(
            chromosome="chr1",
            start_position=100,
            end_position=123,
            strand="+",
            matched_sequence="GAGTCCGAGCAGAAGAAGAA",
            pam_sequence="TGG",
            mismatch_count=0,
            mismatch_positions=[],
            seed_mismatch_count=0,
            score=100.0,
            risk_category="Exact match",
        )
    ]

    json_output = export_to_json(hits)
    parsed = json.loads(json_output)
    assert isinstance(parsed, list)
    assert len(parsed) == 1
    assert parsed[0]["chromosome"] == "chr1"
    assert parsed[0]["score"] == 100.0
    assert parsed[0]["risk_category"] == "Exact match"

