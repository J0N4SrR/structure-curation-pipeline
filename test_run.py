from curation.app import run_backend_pipeline, WIZARD_STAGES
import pandas as pd
raw_bytes = b"CCO\nC(C)C"
report = run_backend_pipeline(raw_bytes, "test.smi", max_mw=1000.0, max_ha=100, deduplicate=True)
print(report.approved)
