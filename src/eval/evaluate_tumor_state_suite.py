"""Leakage-safe shared discovery-support evaluation for frozen and V2 representations."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.metrics import roc_auc_score
from morpheus.src.eval.retrieval_metrics import paired_retrieval_metrics

GROUPS={
 "hallmark_all":None,
 "immune_tme":["ALLOGRAFT_REJECTION","COMPLEMENT","INFLAMMATORY_RESPONSE","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","IL2_STAT5_SIGNALING","IL6_JAK_STAT3_SIGNALING"],
 "cell_state":["APOPTOSIS","E2F_TARGETS","G2M_CHECKPOINT","EPITHELIAL_MESENCHYMAL_TRANSITION","KRAS_SIGNALING_UP","MYC_TARGETS_V1"],
 "metabolism_hypoxia":["HYPOXIA","GLYCOLYSIS","OXIDATIVE_PHOSPHORYLATION","FATTY_ACID_METABOLISM","CHOLESTEROL_HOMEOSTASIS"],
 "mechanobiology":["TGF_BETA_SIGNALING","EPITHELIAL_MESENCHYMAL_TRANSITION","ANGIOGENESIS","APICAL_JUNCTION","COAGULATION"],
}
def _zfit(x,tr):
 m=x[tr].mean(0,keepdims=True); s=x[tr].std(0,keepdims=True); return (x-m)/np.maximum(s,1e-6)
def _corr(a,b):
 vals=[np.corrcoef(a[:,j],b[:,j])[0,1] for j in range(a.shape[1]) if np.std(a[:,j])>1e-8 and np.std(b[:,j])>1e-8]
 return float(np.nanmean(vals)) if vals else float("nan")
def _bootstrap(y,p,fn,n=300,seed=42):
 r=np.random.default_rng(seed); vals=[]
 for _ in range(n):
  ix=r.integers(0,len(y),len(y)); vals.append(fn(y[ix],p[ix]))
 return [float(np.nanquantile(vals,.025)),float(np.nanquantile(vals,.975))]
def _ridge(x,y,tr,ev):
 x=_zfit(x,tr); ymean=y[tr].mean(0,keepdims=True); ystd=np.maximum(y[tr].std(0,keepdims=True),1e-6)
 model=RidgeCV(alphas=np.logspace(-3,4,12)).fit(x[tr],(y[tr]-ymean)/ystd)
 return model.predict(x[ev])*ystd+ymean
def _choose_eval(split, valid=None):
 mask=np.ones(len(split),dtype=bool) if valid is None else valid
 for name in ("test","val"):
  ix=np.where((split==name)&mask)[0]
  if len(ix)>1: return ix,name
 ix=np.where((split!="train")&mask)[0]
 return ix,("non_train" if len(ix)>1 else "none")
def _read_hallmark(root):
 d=pd.read_parquet(root/"data/processed/genesets/msigdb_hallmark_scores.parquet")
 cols=[c for c in d if c!="patient_id"]; d=d.set_index("patient_id"); return d,cols
def _targets(root,ids):
 snv=pd.read_parquet(root/"meta-intersurv/data/omics/processed/snv_features.parquet").rename(columns={"Patient ID":"patient_id"}).set_index("patient_id")
 cnv=pd.read_parquet(root/"meta-intersurv/data/omics/processed/cnv_features.parquet").rename(columns={"Patient ID":"patient_id"}).set_index("patient_id")
 snv=snv.reindex(ids).fillna(0); cnv=cnv.reindex(ids).fillna(0)
 return np.c_[snv.sum(1).to_numpy(),(snv!=0).sum(1).to_numpy(),cnv.abs().mean(1).to_numpy(),(cnv.abs()>=1).mean(1).to_numpy()]
def evaluate(npz_path,root,out):
 z=np.load(npz_path,allow_pickle=True); ids=z["patient_ids"].astype(str); split=z["split"].astype(str); cancers=z["cancers"].astype(str)
 tr=np.where(split=="train")[0]; ev,ev_name=_choose_eval(split)
 x=z["wsi_biology"] if "wsi_biology" in z.files else z["wsi_identity"]
 rows=[]; method=Path(npz_path).parent.name if Path(npz_path).stem == 'representations' else Path(npz_path).stem
 if "rna_identity" in z.files:
  r=paired_retrieval_metrics(z["wsi_identity"][ev],z["rna_identity"][ev],(1,5,10),cancers[ev].tolist(),cancers[ev].tolist())
  rows += [{"method":method,"task":"paired_retrieval","metric":k,"value":float(v)} for k,v in r.items() if isinstance(v,(float,int,np.floating))]
  by=[]
  for c in np.unique(cancers[ev]):
   ix=ev[cancers[ev]==c]
   if len(ix)>1: by.append(paired_retrieval_metrics(z["wsi_identity"][ix],z["rna_identity"][ix],(10,),[c]*len(ix),[c]*len(ix))["recall_at_10"])
  rows.append({"method":method,"task":"within_cancer_retrieval","metric":"macro_r10","value":float(np.mean(by))})
 h,cols=_read_hallmark(root); y=h.reindex(ids)[cols].to_numpy(float); valid=np.isfinite(y).all(1); trh=tr[valid[tr]]; evh,evh_name=_choose_eval(split,valid)
 if len(trh)<2 or len(evh)<2:
  rows.append({"method":method,"task":"hallmark_all","metric":"mean_pearson","value":float("nan"),"n":len(evh),"eval_split":evh_name,"note":"insufficient_target_overlap"})
  return pd.DataFrame(rows)
 pred=_ridge(x,y,trh,evh)
 for group,names in GROUPS.items():
  keep=np.arange(len(cols)) if names is None else np.array([i for i,c in enumerate(cols) if c.replace("HALLMARK_","") in names])
  yy=y[evh][:,keep]; pp=pred[:,keep]; val=_corr(yy,pp)
  rows.append({"method":method,"task":group,"metric":"mean_pearson","value":val,"ci_low":_bootstrap(yy,pp,lambda a,b:_corr(a,b))[0],"ci_high":_bootstrap(yy,pp,lambda a,b:_corr(a,b))[1],"n":len(evh),"eval_split":evh_name})
 g=_targets(root,ids); gp=_ridge(x,g,tr,ev)
 rows.append({"method":method,"task":"genomic_support","metric":"mean_pearson","value":_corr(g[ev],gp),"n":len(ev),"eval_split":ev_name})
 master=pd.read_parquet(root/"morpheus/data/processed/master_patient_table.parquet").set_index("patient_id").reindex(ids)
 if master["survival_event"].notna().sum()>20:
  # Risk score probe; C-index is intentionally deferred to lifelines availability.
  ysur=master["survival_event"].fillna(0).to_numpy(float)[:,None]; sp=_ridge(x,ysur,tr,ev)[:,0]
  rows.append({"method":method,"task":"clinical_support","metric":"event_auc","value":float(roc_auc_score(ysur[ev,0],sp)) if len(np.unique(ysur[ev,0]))>1 else float("nan"),"n":len(ev),"eval_split":ev_name})
 return pd.DataFrame(rows)
def main():
 p=argparse.ArgumentParser(); p.add_argument("--root",default="."); p.add_argument("--inputs",nargs="+",required=True); p.add_argument("--output",required=True); a=p.parse_args()
 root=Path(a.root); frames=[evaluate(Path(x),root,Path(a.output)) for x in a.inputs]; out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
 d=pd.concat(frames,ignore_index=True); d.to_csv(out/"task_suite.csv",index=False); (out/"task_suite.json").write_text(d.to_json(orient="records",indent=2)); print(out/"task_suite.csv")
if __name__=="__main__": main()
