"""Describe the frozen holdout by anchor family, without selecting parameters."""
import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from sea.parameter_validation import paired_bootstrap, summarize


def summarize_families(study, annotations):
    groups = defaultdict(list)
    for row in study['test_query_results']:
        info = annotations.get(row['anchor_target_id'], {})
        family = (info.get('families') or ['Unclassified'])[0]
        groups[family].append(row)
    result = {}
    for name, rows in sorted(groups.items()):
        baseline = np.array([row['baseline'] for row in rows])
        candidate = np.array([row['candidate'] for row in rows])
        # Tiny strata do not support informative percentile-bootstrap intervals.
        comparison = paired_bootstrap(candidate, baseline) if len(rows) >= 10 else None
        result[name] = {
            'queries': len(rows), 'baseline': summarize(baseline), 'candidate': summarize(candidate),
            'mrr_gain_95': [comparison['lower_95'][0], comparison['upper_95'][0]] if comparison else None,
            'interval_scope': 'Descriptive pointwise paired bootstrap; no multiplicity correction.' if comparison
                              else 'Insufficient cohort (fewer than 10 queries); interval omitted.',
        }
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--annotations', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    study = json.loads(args.study.read_text())
    annotation_bytes = args.annotations.read_bytes()
    annotation = json.loads(annotation_bytes)
    if annotation['dataset_id'] != study['dataset_id']:
        raise ValueError('Annotations and holdout must refer to the same source cohort.')
    report = {
        'source_dataset_id': study['dataset_id'], 'budget': study['null_trials'],
        'annotations_sha256': hashlib.sha256(annotation_bytes).hexdigest(),
        'grouping': 'Anchor target primary ChEMBL family; descriptive subgroup analysis, not used for selection or adoption. Not a family-specific validation claim.',
        'families': summarize_families(study, annotation['targets']),
    }
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
