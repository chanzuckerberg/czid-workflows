from sklearn.metrics import auc, precision_recall_curve
from collections import defaultdict
import numpy as np


def truth_aupr(classified_taxa, truth_taxa):
    """
    Compute AUPR (and other metrics) of taxon classifications vs. truth data
    """
    missed_taxa = [tax_id for tax_id in truth_taxa if tax_id not in classified_taxa]
    correctness_labels = [
        1 if tax_id in truth_taxa else 0 for tax_id in classified_taxa
    ]
    correctness_labels += [1] * len(missed_taxa)
    # Using raw abundances as proxies for confidence score per Ye2009 methodology
    # https://www.cell.com/cell/fulltext/S0092-8674(19)30775-5#fig2
    confidence_scores = list(classified_taxa.values())
    confidence_scores += [0] * len(missed_taxa)
    return adjusted_aupr(correctness_labels, confidence_scores, force_monotonic=False)


def truth_l2_norm(classified_taxa, truth_taxa):
    """
    Measure accuracy of taxa relative abundance (L2 norm of vector difference from truth vector)
    """
    truth_sum = sum(truth_taxa.values())
    relative_abundances_diff = [
        truth_taxa[taxon] / truth_sum - classified_taxa.get(taxon, 1e-100)
        for taxon in truth_taxa
    ]
    return np.linalg.norm(relative_abundances_diff, ord=2)


def adjusted_aupr(y_true, y_score, force_monotonic=False):
    # precision_recall_curve wants nonzero probas_pred, so adjust any zeroes to epsilon
    adjusted_y_score = [(1e-100 if score == 0 else score) for score in y_score]
    # Adapted from https://github.com/yesimon/metax_bakeoff_2019/blob/master/plotting/Metagenomics%20Bench.ipynb
    original_precision, recall, thresholds = precision_recall_curve(
        y_true, adjusted_y_score
    )

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
            "recall": list(recall),
            "precision": list(precision),
            "max_f1_recall": recall[argmax_f1],
            "max_f1_precision": precision[argmax_f1],
            "thresholds": list(thresholds),
        }.items()
        if v is not None
    }
