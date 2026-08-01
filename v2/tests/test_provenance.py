from __future__ import annotations

from morpheus.v2.provenance import source_manifest, source_tree_digest


def test_source_manifest_is_stable_for_unchanged_tree():
    assert source_tree_digest() == source_tree_digest()
    manifest = source_manifest(configuration={"seed": 42})
    assert manifest["package"] == "morpheus.v2"
    assert len(manifest["source_tree_sha256"]) == 64
    assert len(manifest["configuration_sha256"]) == 64
