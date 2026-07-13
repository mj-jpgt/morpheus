import numpy as np

from morpheus.src.eval.eval_molecular_prompting import soft_knn_predict


def test_soft_knn_returns_weighted_targets():
    query = np.array([[1.0, 0.0]], dtype=np.float32)
    reference = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    targets = np.array([[2.0, 0.0], [0.0, 2.0]], dtype=np.float32)
    pred = soft_knn_predict(query, reference, targets, k=1, tau=0.1)
    assert pred.shape == (1, 2)
    assert np.allclose(pred[0], [2.0, 0.0])
