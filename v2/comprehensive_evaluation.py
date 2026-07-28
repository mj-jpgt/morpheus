"""Coverage-gated, shared downstream evaluation for frozen MORPHEUS artifacts.

Runs only tasks supported by aligned TCGA WSI/RNA/genome/clinical records and
writes explicit unavailable rows for tasks requiring absent external assays.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from .evaluation import evaluate_few_shot, evaluate_direct_zero_shot
from .tasks import (binary_metrics, cancer_controlled_discordance, harrell_cindex,
                    mean_column_correlation, retrieval_task, train_ridge_predict)
from .contracts import require_state, validate_artifact


def _state(z, view="patient", order=None):
    # A task must consume the state trained for that task.  In particular, do
    # not substitute a retrieval state for biology or a legacy random fusion
    # head for a patient representation.
    choices = {"wsi": ("wsi_biology",), "rna": ("rna_biology",),
               "patient": ("full_patient",), "wsi_identity": ("wsi_identity",),
               "rna_identity": ("rna_identity",), "semantic": ("semantic_wsi",)}
    for key in choices[view]:
        value = require_state(z, key)
        if value is not None: return value if order is None else value[order]
    return None


def _safe_ridge(x, y, tr, te):
    if tr.sum() < 3 or te.sum() < 2: return None
    return train_ridge_predict(x, y, tr, te)


def _cox_risk(x, time, event, tr, te):
    if tr.sum() < 20 or te.sum() < 10 or event[tr].sum() < 3: return None
    try:
        from sksurv.linear_model import CoxPHSurvivalAnalysis
        scaler=StandardScaler().fit(x[tr]); xt=scaler.transform(x); n=min(64, xt.shape[1], tr.sum()-1)
        pca=PCA(n_components=n, random_state=42).fit(xt[tr]); xp=pca.transform(xt)
        target=np.asarray([(bool(e),float(t)) for e,t in zip(event[tr],time[tr])],dtype=[("event","?"),("time","<f8")])
        return CoxPHSurvivalAnalysis(alpha=1.0).fit(xp[tr],target).predict(xp[te])
    except Exception: return None


def _rows_for(path, root, text_labels, text_vectors, hallmark, master, snv, cnv, canonical_ids):
    z=np.load(path,allow_pickle=False); source_ids=z["patient_ids"].astype(str)
    order=np.asarray([dict(zip(source_ids, range(len(source_ids))))[pid] for pid in canonical_ids], dtype=np.int64)
    ids=source_ids[order]; cancer=z["cancers"].astype(str)[order]; split=z["split"].astype(str)[order]
    tr,te=split=="train",split=="test"; method=Path(path).stem; rows=[]
    wsi=_state(z,"wsi",order); rna=_state(z,"rna",order); patient=_state(z,"patient",order); semantic=_state(z,"semantic",order)
    wsi_identity, rna_identity = _state(z,"wsi_identity",order), _state(z,"rna_identity",order)
    # A baseline may have no biology head.  Its explicitly trained identity
    # state is still a legitimate frozen feature baseline for biological probes.
    wsi_task, rna_task = (wsi if wsi is not None else wsi_identity), (rna if rna is not None else rna_identity)
    if wsi_identity is not None and rna_identity is not None:
        r=retrieval_task(wsi_identity[te],rna_identity[te],cancer[te])
        rows += [{"method":method,"task":"retrieval","metric":k,"value":float(v)} for k,v in r.global_metrics.items()]
        rows.append({"method":method,"task":"retrieval","metric":"within_cancer_r10","value":r.within_cancer_r10})
    h=hallmark.reindex(ids); cols=h.columns.tolist(); y=h.to_numpy(float); valid=np.isfinite(y).all(1)
    for view,x in (("wsi",wsi_task),("rna",rna_task),("patient",patient)):
        if x is None: continue
        pred=_safe_ridge(x,y,tr&valid,te&valid)
        if pred is not None:
            rows.append({"method":method,"task":"molecular_prompting","view":view,"metric":"hallmark_mean_pearson","value":mean_column_correlation(y[te&valid],pred)})
    seen,unseen,allc=np.unique(cancer[tr]),np.unique(cancer[te]),np.unique(cancer)
    lookup={c:i for i,c in enumerate(text_labels)} if text_labels is not None else {}
    if semantic is not None and text_vectors is not None and set(allc).issubset(lookup):
        for name,candidates in (("unseen_heldout",unseen),("all_cancers",allc)):
            ordered=np.asarray(sorted(candidates)); m=evaluate_direct_zero_shot(semantic,cancer,split,ordered,text_vectors[[lookup[c] for c in ordered]])
            rows += [{"method":method,"task":"zero_shot_"+name,"view":"semantic_wsi","metric":k,"value":float(v),"protocol":"direct_frozen_plip"} for k,v in m.items()]
    else:
        rows.append({"method":method,"task":"zero_shot","metric":"status","value":float("nan"),"note":"requires_declared_semantic_wsi_and_plip_prompt_cache"})
    for view,x in (("wsi_identity",wsi_identity),("rna_identity",rna_identity),("wsi_biology",wsi),("rna_biology",rna),("patient",patient)):
        if x is None: continue
        for episode in evaluate_few_shot(x,cancer,split,unseen,episodes=100):
            rows += [{"method":method,"task":"few_shot_unseen","view":view,"k":episode["k"],"metric":k,"value":float(v)} for k,v in episode.items() if k.startswith("accuracy_")]
    immune=[c for c in cols if any(t in c for t in ("INTERFERON","ALLOGRAFT","INFLAMMATORY","COMPLEMENT"))]
    iv=np.isfinite(h[immune].to_numpy(float)).all(1); it=h[immune].mean(1).to_numpy(float)
    if wsi_task is not None and rna_task is not None and immune and (tr&iv).sum()>=3 and (te&iv).sum()>=2:
        wp=_safe_ridge(wsi_task,it[:,None],tr&iv,te&iv); rp=_safe_ridge(rna_task,it[:,None],tr&iv,te&iv)
        rows.append({"method":method,"task":"immune_programme","metric":"wsi_pearson","value":mean_column_correlation(it[te&iv,None],wp)})
        d=cancer_controlled_discordance(wp[:,0],rp[:,0],cancer[te&iv])
        rows += [{"method":method,"task":"immune_discordance","metric":k,"value":float(v)} for k,v in d.items() if np.isscalar(v)]
    else: rows.append({"method":method,"task":"immune_programme","metric":"coverage","value":float("nan"),"note":"insufficient_hallmark_coverage"})
    surv=master.reindex(ids); time=surv.survival_time.to_numpy(float); event=surv.survival_event.to_numpy(float); sv=np.isfinite(time)&np.isfinite(event)&(time>0)
    risk=None if patient is None else _cox_risk(patient,time,event.astype(bool),tr&sv,te&sv)
    if risk is not None: rows.append({"method":method,"task":"survival","metric":"harrell_cindex","value":harrell_cindex(time[te&sv],event[te&sv].astype(bool),risk)})
    for name,table in (("snv_burden",snv),("cnv_burden",cnv)):
        values=table.reindex(ids).fillna(0).to_numpy(float); burden=np.abs(values).sum(1); threshold=np.median(burden[tr]); label=(burden>=threshold).astype(int)
        # This is intentionally a morphology-side task: using a full state
        # could include the very genomic modality from which the label was made.
        if wsi_task is not None and len(np.unique(label[tr]))==2 and len(np.unique(label[te]))==2:
            scale=StandardScaler().fit(wsi_task[tr]); clf=LogisticRegression(max_iter=2000,class_weight="balanced").fit(scale.transform(wsi_task[tr]),label[tr]); prob=clf.predict_proba(scale.transform(wsi_task[te]))[:,1]
            rows += [{"method":method,"task":"molecular_phenotype_"+name,"metric":k,"value":v} for k,v in binary_metrics(label[te],prob).items()]
    rows.append({"method":method,"task":"missing_modality_next_test","metric":"status","value":float("nan"),
                 "note":"unavailable: no trained uncertainty head or predeclared uncertainty protocol"})
    return rows


def _shared_cohort(paths: list[str]) -> np.ndarray:
    """Fail closed unless every candidate artifact names exactly the same cohort/split."""
    reference: dict[str, tuple[str, str]] | None = None
    for path in paths:
        validate_artifact(path)
        with np.load(path, allow_pickle=False) as raw:
            ids, cancers, split = raw["patient_ids"].astype(str), raw["cancers"].astype(str), raw["split"].astype(str)
        if len(set(ids)) != len(ids):
            raise ValueError(f"duplicate patient IDs in {path}")
        current = dict(zip(ids, zip(cancers, split)))
        if reference is None:
            reference = current
        elif current != reference:
            raise ValueError(f"{path} does not have the exact same patient/cancer/split cohort as the first artifact")
    if not reference:
        raise ValueError("no artifacts supplied")
    return np.asarray(sorted(reference), dtype=str)


def main():
    p=argparse.ArgumentParser(); p.add_argument("--root",required=True); p.add_argument("--inputs",nargs="+",required=True); p.add_argument("--output",required=True); p.add_argument("--plip-prompts",default=""); a=p.parse_args(); root=Path(a.root); out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
    if a.plip_prompts:
        t=np.load(a.plip_prompts,allow_pickle=False)
        if str(t["model_id"].item()) != "vinid/plip" or not bool(t["direct_image_text_space"].item()):
            raise ValueError("--plip-prompts must be a direct vinid/plip prompt cache")
        text_labels=t["cancers"].astype(str); text_vectors=t["embeddings"].astype(np.float32)
    else:
        text_labels=text_vectors=None
    hallmark=pd.read_parquet(root/"data/processed/genesets/msigdb_hallmark_scores.parquet").set_index("patient_id"); master=pd.read_parquet(root/"morpheus/data/processed/master_patient_table.parquet").set_index("patient_id")
    snv=pd.read_parquet(root/"meta-intersurv/data/omics/processed/snv_features.parquet").rename(columns={"Patient ID":"patient_id"}).set_index("patient_id"); cnv=pd.read_parquet(root/"meta-intersurv/data/omics/processed/cnv_features.parquet").rename(columns={"Patient ID":"patient_id"}).set_index("patient_id")
    canonical_ids=_shared_cohort(a.inputs)
    protocol={"canonical_patient_count":int(len(canonical_ids)),"artifact_paths":[str(Path(x)) for x in a.inputs],
              "few_shot":{"k":[1,2,5,10,20],"episodes":100,"seed":42,"shared_canonical_order":True},
              "bootstrap":"not yet reported; no method-specific resampling"}
    (out/"evaluation_protocol.json").write_text(json.dumps(protocol,indent=2),encoding="utf-8")
    rows=[r for x in a.inputs for r in _rows_for(x,root,text_labels,text_vectors,hallmark,master,snv,cnv,canonical_ids)]
    for task in ("protein_phosphoprotein","spatial","perturbation_drug","organoid"):
        for path in a.inputs:
            rows.append({"method":Path(path).stem,"task":task,"metric":"status","value":float("nan"),"note":"external_paired_data_not_yet_coverage_gated"})
    pd.DataFrame(rows).to_csv(out/"comprehensive_task_suite.csv",index=False); (out/"comprehensive_task_suite.json").write_text(json.dumps(rows,indent=2,default=str)); print(out/"comprehensive_task_suite.csv")
if __name__=="__main__": main()
