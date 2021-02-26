from collections import defaultdict
from pathlib import Path

from ruamel.yaml import YAML
from sklearn.metrics import auc, precision_recall_curve
import numpy as np


def load_benchmarks_yml():
    with open(Path(__file__).parent.parent / "benchmarks.yml") as benchmarks_yml:
        return YAML(typ="safe").load(benchmarks_yml.read())


def adjusted_aupr(y_true, y_score, force_monotonic=False):
    # precision_recall_curve wants nonzero probas_pred, so adjust any zeroes to epsilon
    adjusted_y_score = [(1e-100 if score == 0 else score) for score in y_score]
    # Adapted from https://github.com/yesimon/metax_bakeoff_2019/blob/master/plotting/Metagenomics%20Bench.ipynb
    original_precision, recall, thresholds = precision_recall_curve(y_true, adjusted_y_score)

    precision = original_precision
    if force_monotonic:
        # adjusts precision per each recall valu on the curve
        # it guarantees that the curve is monotonic decreasing
        precision_max_per_recall = defaultdict(float)
        for r, p in zip(recall, original_precision):
            precision_max_per_recall[r] = max(precision_max_per_recall[r], p)
        adjusted_precision = []
        for r, p in zip(recall, original_precision):
            adjusted_precision.append(precision_max_per_recall[r])
        precision = adjusted_precision

    # force start at zero
    if thresholds[0] == 0:
        precision[0] = 0
        recall[0] = recall[1]
        recall = np.insert(recall, 0, 1)
        precision = np.insert(precision, 0, 0)

    aupr = auc(recall, precision)
    # figure which point on the P/R curve maximizes F1
    argmax_f1 = np.argmax((2 * p * r / (p + r) for (p, r) in zip(precision, recall)))
    return {
        k: v
        for k, v in {
            "aupr": aupr,
            "original_precision": original_precision if force_monotonic else None,
            "recall": recall,
            "precision": precision,
            "max_f1_recall": recall[argmax_f1],
            "max_f1_precision": precision[argmax_f1],
            "thresholds": thresholds,
        }.items()
        if v is not None
    }
