import numpy as np

from morpheus.src.eval.retrieval_metrics import paired_retrieval_metrics


def test_perfect_retrieval_metrics():
    x = np.eye(4, dtype=np.float32)
    metrics = paired_retrieval_metrics(x, x, (1, 2))
    assert metrics["recall_at_1"] == 1.0
    assert metrics["recall_at_2"] == 1.0
    assert metrics["median_rank"] == 1.0


def test_retrieval_rejects_unpaired_lengths():
    x = np.eye(2, dtype=np.float32)
    y = np.eye(3, dtype=np.float32)
    try:
        paired_retrieval_metrics(x, y)
    except ValueError as exc:
        assert "equal row counts" in str(exc)
    else:
        raise AssertionError("Expected ValueError")
