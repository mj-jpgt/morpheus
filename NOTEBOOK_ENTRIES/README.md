# Notebook entries staging

Concurrent agents must not edit `NOTEBOOK.md` directly — simultaneous writes silently lose entries.
Write one file per finding here, named `<topic>_<UTC timestamp>.md`, in the running-log format:

```
## <YYYY-MM-DD HH:MM UTC> — <one-line headline>

**Logged:** <date time UTC>. **How obtained:** <command / script / machine>

### Technical
### In plain terms
### Meaning for the claim
### Files / commits
```

The main session merges these into `NOTEBOOK.md` and deletes the merged file.
