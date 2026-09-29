def test_oracle_layout(oracle_dir):
    assert (oracle_dir / "PROVENANCE.sha256").is_file()
    assert any((oracle_dir / "runs" / "agg").iterdir())
    assert (oracle_dir / "archive_outputs" / "crystal").is_dir()
