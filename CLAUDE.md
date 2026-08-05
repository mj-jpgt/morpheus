# Read `PROJECT_GUIDE.md` in full before doing anything else in this repository.

It has the project context, the current per-paper build plan, and sixteen methodology rules — each one earned from a real, documented failure on this project (not style preferences). The rules that bite hardest if skipped:

- **Predeclare before you measure**, in `NOTEBOOK_ENTRIES/PREDECLARED_*.md`, committed before the code runs.
- **Report bad news first**, always.
- **Push favourable results until they break** — more repeats, harder statistics, out-of-sample tests. Do not bank a result you haven't tried to break.
- **Reuse `v2/calibra` canonical functions. Never reimplement a statistic inline** — an AST-scan test will fail you, and it will be right to.
- **Every agent logs everything to `NOTEBOOK_ENTRIES/`**, never edits `NOTEBOOK.md` directly.
- **Shared-tree git discipline**: `git add` only your own changed paths, never `-A`/`.`; if `git pull --rebase` fails on another agent's *uncommitted* files, leave them and just `git push` if you're only ahead.

Branch: `research/rebase-vision`. Four papers, kept deliberately separate: P1 (`paper/P1_CALIBRA_DRAFT.md`, the instrument, strongest result), P2 (`paper/P2_RANK_DRAFT.md`, a negative methods result), P3 (perturbation-basis → causal attribution), P4 (the promptable multimodal system, earliest stage).

If your work changes what's open or done in `PROJECT_GUIDE.md` §3, **update that section before you finish** — it is a living document, not a snapshot.
