## 2026-08-03 12:10 UTC — D2 decided: the PBS deficit survives on 40 untrained targets, and the contamination worry ran the other way

**Logged:** 2026-08-03 12:10 UTC. **How obtained:** Lambda A100 `150.136.45.194`. Full D2 re-run
`d2_v3` — three concurrent `phase_d d2` streams (`--seeds 42|43|44 --restrict-to-split
--pbs-components 128 --analysis-role primary --token-budget 8192`), 40 epochs/arm, 6/6
`TRAIN_SUCCESS`, CALIBRA G4 `gates_pass: true` on all three seeds. Readouts via
`morpheus.v2.research.rebase.d2_compare`, 2,000 repeats, one invocation per seed-pair per target
group, `OMP_NUM_THREADS=1`. All outputs on persistent NFS under `~/e0_run/d2_v3/`.

### Technical
Δ = PBS − Hallmark, top canonical correlation at 16 components on residualised held-out
`wsi_biology`, n_test = 2,766, cancer+TSS residualisation (84 sites).

| target set | n | seed 42 | seed 43 | seed 44 |
|---|---:|---:|---:|---:|
| unrestricted (all non-control) | 90 | −0.1359 | −0.1077 | −0.1192 |
| **untrained (heldout_pathway+immune_tme+tumour_state)** | **40** | **−0.1325** | **−0.1089** | **−0.1226** |
| `hallmark_in_training` only | 50 | −0.1091 | −0.0787 | −0.1149 |
| **`random_control` (negative control)** | 90 | **−0.0099** | **−0.0280** | **−0.0268** |

Untrained-40 CIs, patient / cancer:
s42 [−0.1605,−0.0993] / [−0.1792,−0.0632]; s43 [−0.1460,−0.0749] / [−0.1623,−0.0118];
s44 [−0.1502,−0.0866] / [−0.1653,−0.0411]. **Both exclude zero in 3/3.**
Random-control cancer CIs cover zero in 3/3; patient CIs cover zero in 2/3.

Three results matter beyond the headline:

1. **Stratifying barely moves the number.** Untrained-40 minus unrestricted is +0.0034 / −0.0012 /
   −0.0034. The 50 `hallmark_in_training` targets were not carrying the effect.
2. **The contamination ran the *opposite* way to the worry.** Scored on Hallmark's own supervision
   targets alone the gap is *smaller* (−0.1091/−0.0787/−0.1149) than on the untrained 40. Putting
   H's homework on the exam **understated** PBS's deficit.
3. **Effective rank contradicts the capacity explanation.** Residualised `wsi_biology` effective
   rank, H vs PBS: 23.39/14.87 (s42), 28.77/**34.12** (s43), 9.14/9.11 (s44). In seed 43 PBS has
   *higher* rank and still loses by −0.1089; in seed 44 they are equal to two decimals and PBS still
   loses by −0.1226.

Negative control, stated honestly: the arm gap on random targets is 4–13× smaller than on real
targets and mostly indistinguishable from zero, so the instrument is not manufacturing arm
differences — but it is not exactly zero and points the same way, so ~10–20% of the headline may be
generic representation quality rather than supervision content. The *absolute* level on random
controls (0.44–0.47) is not evidence of a broken instrument: a 200-draw permutation null puts chance
at **0.140** for every group, and random *gene sets* still track global expression covariance.

**A number that contradicts the earlier run.** Re-exporting the surviving seed-42 Hallmark
checkpoint reproduces the lost point estimate to five significant figures (**0.58612** vs recorded
0.5861), so the readout path is deterministic. But **retraining seed 42 with identical configuration
gives 0.6214**, with effective rank 23.39 against the recovered checkpoint's 8.68. Training is not
seed-reproducible on this stack. The paired contrast is unaffected (it is within-run), but no
individual D2 level should ever be quoted as reproducible from the seed.

### In plain terms
The fear was that the earlier verdict against the perturbation-based method was rigged, because more
than half the exam consisted of the rival method's own homework. We rebuilt the whole experiment and
re-marked it on forty questions neither method had ever seen. The gap is the same size, in all three
repeats. Marking both on deliberately meaningless questions leaves them nearly tied, so the marking
scheme is not simply biased toward one method. And the obvious alternative explanation — that the
winner merely had more room to express itself — fails, because in one repeat the loser had more room
and still lost. Separately: training the same model twice with the same seed does not give the same
model, so only differences measured within one run can be trusted.

### Meaning for the claim
**D2 is settled and quotable.** Perturbation-basis supervision is genuinely worse than Hallmark
supervision for held-out molecular structure, ~0.11–0.13 top-CCA in 3/3 seeds on targets neither arm
trained on. P3's headline hypothesis is refuted by its own predeclared test. The previous
`D2_RESULT.md` verdict of "pending" is discharged, and it resolves toward the first branch it named:
the gap was not measuring which arm's training targets were on the exam.

Two caveats must travel with the number: the negative control is small but non-zero, and point
estimates are not seed-reproducible — quote the paired difference, never the level.

### Files / commits
- `v2/research/rebase/nature/D2_RESULT.md` (updated; copy at `~/e0_run/d2_v3/D2_RESULT.md`)
- `~/e0_run/d2_v3/bootstrap/D2_{unrestricted,untrained40,random_control,hallmark_in_training}_seed{42,43,44}.json`
- `~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json` — points, permutation nulls, effective rank
- `~/e0_run/d2_v3/d2_v3_s{42,43,44}/` — checkpoints, artifacts, CALIBRA gates
- `~/e0_run/d2_v3/recovered_artifacts/` — re-exported seed-42 H (epoch 39) and partial PBS (epoch 13)
