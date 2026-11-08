"""Run with Pixi after setting SEA_API_KEY and SEA_DATASET_ID.

SEA_URL optionally selects a server (default: http://127.0.0.1:8000).
SEA_QUERY optionally replaces the example SMILES.
"""
import json
import os

from sea import SEAClient


with SEAClient() as client:
    result = client.predict(
        os.getenv("SEA_QUERY", "COc1ccc2[nH]cc(c2c1)C(CC)CN"),
        dataset_id=os.environ["SEA_DATASET_ID"], top_k=10,
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))
