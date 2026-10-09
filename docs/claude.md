# Using smORF-Bullpen with Claude

smORF-Bullpen includes a [Claude Code](https://code.claude.com) plugin that helps you set up
runs, check on them and interpret the results by asking Claude in plain language. The plugin
gives Claude the pipeline's configuration options, failure modes and output columns, so it can
work from your files rather than from general knowledge.

## Contents

- [Install](#install)
- [What the plugin does](#what-the-plugin-does)
- [Example session](#example-session)
- [Working safely on a cluster](#working-safely-on-a-cluster)
- [Without the plugin](#without-the-plugin)
- [Updating and removing](#updating-and-removing)

## Install

From inside a Claude Code session:

```
/plugin install smorf-bullpen --marketplace brendan-miller-salk/smORF-Bullpen
```

Or in a terminal:

```bash
claude plugin install smorf-bullpen --marketplace brendan-miller-salk/smORF-Bullpen
```

On Claude Code versions older than 2.1.275, add the marketplace first:

```bash
claude plugin marketplace add brendan-miller-salk/smORF-Bullpen
claude plugin install smorf-bullpen@smorf-bullpen
```

The plugin is available in every directory once installed, including run directories outside
the repository. Start a new session, or run `/reload-plugins`, for it to take effect.

## What the plugin does

The plugin adds three skills. Claude uses them automatically when your request matches, or you
can call one directly.

| Skill | Use it to | Example requests |
|---|---|---|
| `/smorf-bullpen:setup` | Go from FASTQs (or an ORF GTF) to a launched run: run directory, samplesheet, config, module choice, dry run and submission | "Set up a smORF-Bullpen run for the FASTQs in /data/brain_riboseq" · "I have RNA-seq too, should I use assembly?" |
| `/smorf-bullpen:status` | Check progress, find why a run stopped, and get the fix and the command to resume | "Is my run still going?" · "My run failed, what happened?" |
| `/smorf-bullpen:results` | Summarise a finished run and query the results | "How many novel smORFs passed?" · "Which uORFs are conserved in human?" · "Why was ORF_1234 excluded?" |

Each skill accepts a path, for example `/smorf-bullpen:status ~/runs/brain`. Without one, it
uses the current directory.

## Example session

```
$ cd ~/runs
$ claude
> Set up a smORF-Bullpen run for the Ribo-seq FASTQs in /data/brain_riboseq.
  I also have matched RNA-seq in /data/brain_rnaseq.
```

Claude checks your Snakemake and Apptainer setup, asks where your references are (or offers to
build them), proposes sample names from the FASTQ files for you to confirm, asks which modules
you want, writes `config/bullpen.yaml`, and shows the dry-run command for your approval. After
the dry run looks right, it shows the launch command, again for approval.

Later, from the run directory:

```
> /smorf-bullpen:status
> /smorf-bullpen:results How many Primary smORFs are novel, and how many of those are conserved in human?
```

## Working safely on a cluster

- **You approve every submission.** The skills show each `sbatch` command and wait for you
  before submitting. Claude Code also asks permission before running commands, depending on
  your permission settings.
- **Heavy work runs as jobs.** The skills submit the pipeline through Slurm and do not run
  alignment or other heavy tools on a login node.
- **Your files stay where they are.** The skills write only to the run directory (and, if you
  ask for one, a site profile). Result files are read, never modified.
- **Check the account and partition.** Claude looks up the accounts available to you with
  `sacctmgr` and asks which to use; the head job requests 2 days.

## Without the plugin

Claude can also help without the plugin if you point it at the documentation. Start Claude in
your run directory and paste one of these:

```
Read https://github.com/brendan-miller-salk/smORF-Bullpen/blob/main/docs/usage.md and help me
set up a smORF-Bullpen run for the FASTQs in <folder>. Show me every file before writing it and
every sbatch command before submitting it.
```

```
My smORF-Bullpen run in this directory stopped. Read the newest logs/bullpen_head_*.out, find
the first "Error in rule", read that rule's log under the output directory, and explain the
cause and fix. Failure modes are described in
https://github.com/brendan-miller-salk/smORF-Bullpen/blob/main/docs/output.md#quality-control-reports
```

```
Using https://github.com/brendan-miller-salk/smORF-Bullpen/blob/main/docs/output.md for column
definitions, summarise the results in <outdir>/orfs.csv and summary.csv: Primary and
Comprehensive smORF counts, types, and any QC warnings.
```

The plugin is more reliable, because it carries the exact failure messages, column definitions
and module dependencies.

## Updating and removing

```bash
claude plugin marketplace update smorf-bullpen     # fetch the latest catalogue
claude plugin update smorf-bullpen@smorf-bullpen   # update the plugin
claude plugin uninstall smorf-bullpen@smorf-bullpen
```

To update automatically, open `/plugin`, go to **Marketplaces** and enable auto-update for
`smorf-bullpen`.
