import numpy as np
import pytest

from morpheus.src.encoders.bulkformer_encoder import write_embedding_store


def test_write_embedding_store_requires_512_dimensions(tmp_path):
    with pytest.raises(ValueError) as exc:
        write_embedding_store(["TCGA-AA-0001"], np.zeros((1, 128), dtype=np.float32), tmp_path / "x.parquet")
    assert "512D" in str(exc.value)
