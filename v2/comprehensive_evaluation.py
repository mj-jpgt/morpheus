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
from .descriptions import load_tcga_descriptions
from .evaluation import evaluate_few_shot, evaluate_zero_shot
from .tasks import (binary_metrics, cancer_controlled_discordance, harrell_cindex,
                    mean_column_correlation, retrieval_task, train_ridge_predict)
from .text_prototypes import embed_descriptions


def _state(z, view="patient"):
    choices = {"wsi": ("wsi_biology", "wsi_identity"), "rna": ("rna_identity",),
               "patient": ("full_patient", "wsi_biology", "wsi_identity")}
    for key in choices[view]:
        if key in z.files: return z[key].astype(np.float32)
    raise ValueError(f"artifact missing {view} state")


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


def _rows_for(path, root, text_labels, text_vectors, hallmark, master, snv, cnv):
    z=np.load(path,allow_pickle=False); ids=z["patient_ids"].astype(str); cancer=z["cancers"].astype(str); split=z["split"].astype(str)
    tr,te=split=="train",split=="test"; method=Path(path).stem; rows=[]
    wsi=_state(z,"wsi"); rna=_state(z,"rna"); patient=_state(z,"patient")
    if "wsi_identity" in z.files and "rna_identity" in z.files:
        r=retrieval_task(z["wsi_identity"][te],z["rna_identity"][te],cancer[te])
        rows += [{"method":method,"task":"retrieval","metric":k,"value":float(v)} for k,v in r.global_metrics.items()]
        rows.append({"method":method,"task":"retrieval","metric":"within_cancer_r10","value":r.within_cancer_r10})
    h=hallmark.reindex(ids); cols=h.columns.tolist(); y=h.to_numpy(float); valid=np.isfinite(y).all(1)
    for view,x in (("wsi",wsi),("rna",rna),("patient",patient)):
        pred=_safe_ridge(x,y,tr&valid,te&valid)
        if pred is not None:
            rows.append({"method":method,"task":"molecular_prompting","view":view,"metric":"hallmark_mean_pearson","value":mean_column_correlation(y[te&valid],pred)})
    seen,unseen,allc=np.unique(cancer[tr]),np.unique(cancer[te]),np.unique(cancer)
    lookup={c:i for i,c in enumerate(text_labels)}
    if set(seen).issubset(lookup) and set(allc).issubset(lookup):
        for view,x in (("wsi",wsi),("rna",rna),("patient",patient)):
            for name,candidates in (("unseen_heldout",unseen),("all_cancers",allc)):
                m=evaluate_zero_shot(x,cancer,split,seen,text_vectors[[lookup[c] for c in seen]],candidates,text_vectors[[lookup[c] for c in candidates]])
                rows += [{"method":method,"task":"zero_shot_"+name,"view":view,"metric":k,"value":float(v)} for k,v in m.items()]
            for episode in evaluate_few_shot(x,cancer,split,unseen,episodes=100):
                rows += [{"method":method,"task":"few_shot_unseen","view":view,"k":episode["k"],"metric":k,"value":float(v)} for k,v in episode.items() if k.startswith("accuracy_")]
    immune=[c for c in cols if any(t in c for t in ("INTERFERON","ALLOGRAFT","INFLAMMATORY","COMPLEMENT"))]
    iv=np.isfinite(h[immune].to_numpy(float)).all(1); it=h[immune].mean(1).to_numpy(float)
    if immune and (tr&iv).sum()>=3 and (te&iv).sum()>=2:
        wp=_safe_ridge(wsi,it[:,None],tr&iv,te&iv); rp=_safe_ridge(rna,it[:,None],tr&iv,te&iv)
        rows.append({"method":method,"task":"immune_programme","metric":"wsi_pearson","value":mean_column_correlation(it[te&iv,None],wp)})
        d=cancer_controlled_discordance(wp[:,0],rp[:,0],cancer[te&iv])
        rows += [{"method":method,"task":"immune_discordance","metric":k,"value":float(v)} for k,v in d.items() if np.isscalar(v)]
    else: rows.append({"method":method,"task":"immune_programme","metric":"coverage","value":float("nan"),"note":"insufficient_hallmark_coverage"})
    surv=master.reindex(ids); time=surv.survival_time.to_numpy(float); event=surv.survival_event.to_numpy(float); sv=np.isfinite(time)&np.isfinite(event)&(time>0)
    risk=_cox_risk(patient,time,event.astype(bool),tr&sv,te&sv)
    if risk is not None: rows.append({"method":method,"task":"survival","metric":"harrell_cindex","value":harrell_cindex(time[te&sv],event[te&sv].astype(bool),risk)})
    for name,table in (("snv_burden",snv),("cnv_burden",cnv)):
        values=table.reindex(ids).fillna(0).to_numpy(float); burden=np.abs(values).sum(1); threshold=np.median(burden[tr]); label=(burden>=threshold).astype(int)
        if len(np.unique(label[tr]))==2 and len(np.unique(label[te]))==2:
            scale=StandardScaler().fit(patient[tr]); clf=LogisticRegression(max_iter=2000,class_weight="balanced").fit(scale.transform(patient[tr]),label[tr]); prob=clf.predict_proba(scale.transform(patient[te]))[:,1]
            rows += [{"method":method,"task":"molecular_phenotype_"+name,"metric":k,"value":v} for k,v in binary_metrics(label[te],prob).items()]
    if valid[te].sum()>=2:
        pw=_safe_ridge(wsi,y,tr&valid,te&valid); pf=_safe_ridge(np.c_[wsi,rna],y,tr&valid,te&valid)
        if pw is not None and pf is not None:
            ew=np.abs(y[te&valid]-pw).mean(1); ef=np.abs(y[te&valid]-pf).mean(1); train_norm=np.linalg.norm(wsi[tr&valid],axis=1); q=np.linalg.norm(wsi[te&valid],axis=1)
            high=q>=np.quantile(q,.8); rows.append({"method":method,"task":"missing_modality_next_test","metric":"rna_value_over_wsi_high_uncertainty","value":float((ew[high]-ef[high]).mean())})
    return rows


def main():
    p=argparse.ArgumentParser(); p.add_argument("--root",required=True); p.add_argument("--inputs",nargs="+",required=True); p.add_argument("--output",required=True); a=p.parse_args(); root=Path(a.root); out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
    cache=out/"biomedbert_prompts.npz"; desc=load_tcga_descriptions(root/"code/morpheus/v2/tcga_cancer_descriptions.yaml")
    if not cache.exists(): embed_descriptions(desc,cache)
    t=np.load(cache); hallmark=pd.read_parquet(root/"data/processed/genesets/msigdb_hallmark_scores.parquet").set_index("patient_id"); master=pd.read_parquet(root/"morpheus/data/processed/master_patient_table.parquet").set_index("patient_id")
    snv=pd.read_parquet(root/"meta-intersurv/data/omics/processed/snv_features.parquet").rename(columns={"Patient ID":"patient_id"}).set_index("patient_id"); cnv=pd.read_parquet(root/"meta-intersurv/data/omics/processed/cnv_features.parquet").rename(columns={"Patient ID":"patient_id"}).set_index("patient_id")
    rows=[r for x in a.inputs for r in _rows_for(x,root,t["cancers"].astype(str),t["embeddings"].astype(np.float32),hallmark,master,snv,cnv)]
    for task in ("protein_phosphoprotein","spatial","perturbation_drug","organoid"):
        rows.append({"method":"all","task":task,"metric":"status","value":float("nan"),"note":"external_paired_data_not_yet_coverage_gated"})
    pd.DataFrame(rows).to_csv(out/"comprehensive_task_suite.csv",index=False); (out/"comprehensive_task_suite.json").write_text(json.dumps(rows,indent=2,default=str)); print(out/"comprehensive_task_suite.csv")
if __name__=="__main__": main()
