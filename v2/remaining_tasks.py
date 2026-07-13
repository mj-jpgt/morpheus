"""Frozen-representation zero/few-shot, survival, and discordance evaluation."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
from .descriptions import load_tcga_descriptions
from .text_prototypes import embed_descriptions
from .evaluation import evaluate_few_shot, evaluate_zero_shot
from .tasks import cancer_controlled_discordance, harrell_cindex, mean_column_correlation, train_ridge_predict

def state(x):
    for key in ("full_patient","wsi_biology","wsi_identity"):
        if key in x.files: return x[key].astype(np.float32)
    raise ValueError("artifact has no usable state")

def main():
    p=argparse.ArgumentParser(); p.add_argument("--root",required=True); p.add_argument("--inputs",nargs="+",required=True); p.add_argument("--output",required=True); a=p.parse_args()
    root,out=Path(a.root),Path(a.output); out.mkdir(parents=True,exist_ok=True)
    desc=load_tcga_descriptions(root/"code/morpheus/v2/tcga_cancer_descriptions.yaml")
    cache=out/"biomedbert_tcga_prompts.npz"
    if not cache.exists(): embed_descriptions(desc,cache)
    text=np.load(cache); text_labels=text["cancers"].astype(str); text_vectors=text["embeddings"].astype(np.float32)
    hallmark=pd.read_parquet(root/"data/processed/genesets/msigdb_hallmark_scores.parquet").set_index("patient_id")
    master=pd.read_parquet(root/"morpheus/data/processed/master_patient_table.parquet").set_index("patient_id")
    rows=[]
    for path in map(Path,a.inputs):
        z=np.load(path,allow_pickle=False); ids=z["patient_ids"].astype(str); cancer=z["cancers"].astype(str); split=z["split"].astype(str); tr=split=="train"; te=split=="test"; f=state(z); method=path.stem
        seen=np.unique(cancer[tr]); unseen=np.unique(cancer[te]); allc=np.unique(cancer)
        for name,candidates in (("unseen_22",unseen),("all_33",allc)):
            si=np.asarray([np.where(text_labels==c)[0][0] for c in seen]); ci=np.asarray([np.where(text_labels==c)[0][0] for c in candidates])
            metrics=evaluate_zero_shot(f,cancer,split,seen,text_vectors[si],candidates,text_vectors[ci])
            rows += [{"method":method,"task":"zero_shot_"+name,"metric":k,"value":float(v)} for k,v in metrics.items()]
        for row in evaluate_few_shot(f,cancer,split,unseen,episodes=100):
            rows += [{"method":method,"task":"few_shot_unseen","metric":k,"value":float(v)} for k,v in row.items() if k not in {"k","episodes"}]
        y=hallmark.reindex(ids); cols=[c for c in hallmark if c.startswith("HALLMARK_")]; immune=[c for c in cols if any(q in c for q in ("INTERFERON","ALLOGRAFT","INFLAMMATORY","COMPLEMENT"))]
        valid=np.isfinite(y[immune].to_numpy(float)).all(1); target=y[immune].mean(1).to_numpy(float); pred=train_ridge_predict(f,target[:,None],tr&valid,te&valid)[:,0]
        rows.append({"method":method,"task":"immune_programme","metric":"pearson","value":mean_column_correlation(target[te&valid,None],pred[:,None])})
        if "wsi_identity" in z.files and "rna_identity" in z.files:
            wp=train_ridge_predict(z["wsi_identity"],target[:,None],tr&valid,te&valid)[:,0]; rp=train_ridge_predict(z["rna_identity"],target[:,None],tr&valid,te&valid)[:,0]
            d=cancer_controlled_discordance(wp,rp,cancer[te&valid]); rows += [{"method":method,"task":"immune_discordance","metric":k,"value":float(v)} for k,v in d.items() if np.isscalar(v)]
        surv=master.reindex(ids); ok=tr & surv.survival_time.notna().to_numpy() & surv.survival_event.notna().to_numpy(); ev=te & surv.survival_time.notna().to_numpy() & surv.survival_event.notna().to_numpy()
        if ok.sum()>10 and ev.sum()>10:
            risk=train_ridge_predict(f,surv.survival_event.fillna(0).to_numpy()[:,None],ok,ev)[:,0]
            rows.append({"method":method,"task":"survival","metric":"harrell_cindex","value":harrell_cindex(surv.survival_time.to_numpy()[ev],surv.survival_event.to_numpy(bool)[ev],risk)})
    pd.DataFrame(rows).to_csv(out/"remaining_task_suite.csv",index=False); (out/"remaining_task_suite.json").write_text(json.dumps(rows,indent=2)); print(out/"remaining_task_suite.csv")
if __name__=="__main__": main()
