import numpy as np
import pytest
import json
from hashlib import sha256
from types import SimpleNamespace

from morpheus.v2.pbs import LegibilityOperator, ReferenceDictionary
from morpheus.v2.build_pbs_targets import fit_development_expression_transform
from morpheus.v2.runner import _attach_programme_matrix, attach_external_programme_targets


def test_reference_dictionary_requires_exact_gene_identity_for_image_free_codes():
    rng = np.random.default_rng(0)
    dictionary = ReferenceDictionary.fit(rng.normal(size=(20, 8)), [f"G{i}" for i in range(8)], [f"A{i}" for i in range(20)], n_components=4)
    codes = dictionary.encode_expression(rng.normal(size=(6, 8)), [f"G{i}" for i in range(8)])
    assert codes.shape == (6, 4)
    with pytest.raises(ValueError, match="exactly match"):
        dictionary.encode_expression(rng.normal(size=(6, 8)), [f"G{i}" for i in reversed(range(8))])


def test_legibility_operator_is_cross_cancer_fit_and_nonnegative():
    rng = np.random.default_rng(1)
    wsi = rng.normal(size=(60, 6)); codes = wsi[:, :3] @ rng.normal(size=(3, 4)) + .1 * rng.normal(size=(60, 4))
    operator = LegibilityOperator.fit(wsi, codes, np.repeat(["A", "B", "C"], 20))
    assert operator.weights.shape == (4,)
    assert np.all((0 <= operator.weights) & (operator.weights <= 1))


def test_pbs_bulk_expression_transform_is_fit_on_development_rows_only():
    expression = np.asarray([[1., 3.], [3., 7.], [100., -50.]], dtype=np.float32)
    transformed, mean, scale = fit_development_expression_transform(expression, np.asarray([True, True, False]))
    assert np.allclose(mean, [2., 5.])
    assert np.allclose(scale, [1., 2.])
    assert np.allclose(transformed[:2], [[-1., -1.], [1., 1.]])
    # Altering held-out data cannot change the fitted transform.
    altered = expression.copy(); altered[2] = [-999., 999.]
    _, repeat_mean, repeat_scale = fit_development_expression_transform(altered, np.asarray([True, True, False]))
    assert np.allclose(mean, repeat_mean) and np.allclose(scale, repeat_scale)


def test_pbs_target_artifact_requires_matching_schema_and_provenance(tmp_path):
    patient_ids = np.asarray(["TCGA-AA-0001", "TCGA-BB-0002"])
    split = np.asarray(["train", "train"])
    scores = np.arange(128, dtype=np.float32).reshape(2, 64) + 1.0
    names = np.asarray([f"PBS_{i:03d}" for i in range(64)])
    genes, singular = np.asarray(["A", "B", "C"]), np.arange(1, 65, dtype=np.float32)
    gene_basis = np.arange(192, dtype=np.float32).reshape(3, 64) + 1.0
    atom_ids = np.asarray(["atom_a", "atom_b"])
    atom_coordinates = np.arange(128, dtype=np.float32).reshape(2, 64) + 1.0
    string_digest = lambda value: sha256("\n".join(map(str, value)).encode()).hexdigest()
    array_digest = lambda value: sha256(np.ascontiguousarray(value).view(np.uint8)).hexdigest()
    manifest = {"target_kind": "external_perturbation_dictionary_coordinates", "patient_id_digest": string_digest(patient_ids),
                "split_digest": string_digest(split), "overlap_gene_digest": string_digest(genes),
                "scores_sha256": array_digest(scores), "singular_values_sha256": array_digest(singular),
                "gene_basis_sha256": array_digest(gene_basis), "atom_coordinates_sha256": array_digest(atom_coordinates),
                "atom_id_digest": string_digest(atom_ids)}
    manifest["fit_patient_id_digest"] = string_digest(patient_ids)
    path = tmp_path / "pbs.npz"
    np.savez_compressed(path, patient_ids=patient_ids, split=split, scores=scores, target_names=names,
                        target_groups=np.asarray(["PBS"] * 64), genes=genes, singular_values=singular,
                        gene_basis=gene_basis, atom_coordinates=atom_coordinates, atom_ids=atom_ids,
                        manifest_json=np.asarray(json.dumps(manifest)))
    data = SimpleNamespace(patient_ids=patient_ids.tolist(), split=split, cancers=np.asarray(["A", "A"]))
    result = attach_external_programme_targets(data, str(path), np.asarray([True, True]))
    assert result["target_dimension"] == 64 and data._v2_programmes.shape == (2, 64)
    manifest["target_kind"] = "hallmark"
    np.savez_compressed(path, patient_ids=patient_ids, split=split, scores=scores, target_names=names,
                        target_groups=np.asarray(["PBS"] * 64), genes=genes, singular_values=singular,
                        gene_basis=gene_basis, atom_coordinates=atom_coordinates, atom_ids=atom_ids,
                        manifest_json=np.asarray(json.dumps(manifest)))
    with pytest.raises(ValueError, match="not a PBS"):
        attach_external_programme_targets(data, str(path), np.asarray([True, True]))


def test_programme_targets_reject_development_constant_axis_even_if_test_varies():
    data = SimpleNamespace(patient_ids=["A", "B", "C"], split=np.asarray(["train", "train", "test"]),
                           cancers=np.asarray(["A", "A", "B"]))
    values = np.asarray([[1., 1.], [1., 2.], [99., 3.]], dtype=np.float32)
    with pytest.raises(ValueError, match="development-constant"):
        _attach_programme_matrix(data, values, np.asarray([True, True, True]), np.asarray([True, True, False]))
