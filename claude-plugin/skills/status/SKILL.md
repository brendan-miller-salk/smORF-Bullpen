---
name: status
description: Check on a smORF-Bullpen run and diagnose failures. Use when the user asks whether their run is still going, how far along it is, why it failed or stopped, or how to resume it.
argument-hint: "[run directory]"
---

# Check and diagnose a smORF-Bullpen run

Work out what state the run is in, and if something failed, find the first real error, explain it
in plain terms and propose the fix. Read logs; do not change the run's outputs.

## 1. Locate the run

Use `$ARGUMENTS` as the run directory if given, otherwise the current directory. A run directory
has `config/bullpen.yaml` and a `logs/` folder. Read the config for `outdir` and `index_dir`
(both relative to the run directory unless absolute).

## 2. Establish the state

```bash
ls -t logs/bullpen_head_*.out | head -3       # head job logs, newest first
squeue -u $USER                               # is the head job (bullpen_head) still queued or running?
sacct -j <jobid> --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode
```

Read the end of the newest head log (`tail -n 80`). Snakemake reports progress as
`N of M steps (X%) done`.

| What the log shows | State |
|---|---|
| Head job in `squeue`, progress lines advancing | Running; report the latest `N of M steps` and which rules are active |
| `Nothing to be done` or the final `100%` line, head job gone | Finished; suggest the `results` skill |
| `Error in rule <name>` | A step failed; go to step 3 |
| `WorkflowError` before any job is submitted | Config problem; see [failures.md](failures.md#config-errors) |
| `LockException` / `Directory cannot be locked` | A previous head job was killed; see step 4 |
| Head log stops with no error, head job gone | Head job hit its time limit or was killed; check `sacct` State |

## 3. Find the failing step

For each `Error in rule <name>` block, the log names a `log:` file and a Slurm job ID. Read, in
order:

1. The rule's own log under `<outdir>/logs/` (or `<index_dir>/logs/` for reference-building
   steps). This is where the tool's error message is.
2. The Slurm log at `.snakemake/slurm_logs/rule_<name>/` in the run directory.
3. `sacct -j <slurm job id> --format=JobID,State,Elapsed,MaxRSS,ReqMem,ExitCode` to tell an
   error from `OUT_OF_MEMORY` or `TIMEOUT`.

Focus on the first failure; later ones are often consequences. Match the message against
[failures.md](failures.md), which lists every deliberate stop in the pipeline, what it means and
how to relax or fix it.

## 4. Propose the fix and the resume

- Explain the cause in one or two sentences, quoting the key log line.
- Give the specific change (a config key, a samplesheet fix, a `set-resources` override) and
  show the edit before making it. Never relax a QC gate (`strand_check: warn`,
  `psite_check: warn`) without telling the user what the gate protects against.
- Resuming re-runs the same launch command; finished steps are skipped. Show it, for example
  `sbatch -A <account> -p <partition> $BULLPEN_REPO/profiles/slurm/run_head.sbatch config/bullpen.yaml`,
  and submit only after the user approves.
- For a stale lock, confirm with `squeue` that no head job for this run is alive, then propose
  running Snakemake once with `--unlock` (append it to the launch command).
- A script edit is not noticed by Snakemake; add `--forcerun <rule>` to redo that step.

## 5. Report

State, then cause, then fix, then the exact command. Keep log excerpts short. If the cause is
unclear, say what you checked and what to look at next rather than guessing.
