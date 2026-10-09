# Contributing to smORF-Bullpen

This guide covers setting up a machine to work on the pipeline, committing and pushing changes,
and keeping several machines in sync.

## Contents

- [Setting up a new machine](#setting-up-a-new-machine)
- [Making a change](#making-a-change)
- [Working from more than one machine](#working-from-more-than-one-machine)
- [Larger changes](#larger-changes)
- [What stays out of the repository](#what-stays-out-of-the-repository)
- [Commit messages](#commit-messages)

## Setting up a new machine

Do this once on each computer or cluster you work from.

**1. Set your identity.** Use the same name and email everywhere, so that every commit is
attributed to the same GitHub account:

```bash
git config --global user.name "Brendan Miller"
git config --global user.email "you@example.org"
```

**2. Create an SSH key for this machine** and add it to GitHub. Each machine gets its own key;
never copy a private key between machines.

```bash
ssh-keygen -t ed25519 -f ~/.ssh/github -C "you@example.org"
cat ~/.ssh/github.pub
```

Paste the printed line at [github.com/settings/keys](https://github.com/settings/keys)
(**New SSH key**, type **Authentication Key**). Then tell SSH to use it for GitHub by adding
these lines to `~/.ssh/config`:

```
Host github.com
    IdentityFile ~/.ssh/github
    IdentitiesOnly yes
```

Check the connection. GitHub should greet you by username:

```bash
ssh -T git@github.com
```

**3. Clone the repository:**

```bash
git clone git@github.com:brendan-miller-salk/smORF-Bullpen.git
cd smORF-Bullpen
```

**4. Restore anything machine-specific.** Site profiles and run configs are not in the
repository (see [below](#what-stays-out-of-the-repository)), so copy them over from your backup
or set them up again from `profiles/example/`.

## Making a change

**1. Start from the latest version:**

```bash
git pull
```

**2. Edit, then run the tests:**

```bash
for t in tests/test_*.py; do python "$t" > /dev/null && echo "ok   $t" || echo "FAIL $t"; done
```

For changes to rules or the Snakefile, a dry run from a run directory also catches syntax and
wiring errors:

```bash
snakemake -s $BULLPEN_REPO/workflow/Snakefile --configfile config/bullpen.yaml -n
```

**3. Review what changed:**

```bash
git status          # which files changed
git diff            # the changes themselves
```

**4. Commit and push:**

```bash
git add <files>     # or `git add -A` for everything listed by git status
git commit -m "Short description of the change"
git push
```

## Working from more than one machine

Every machine is a full copy of the repository, and GitHub holds the shared version. As long as
you **pull before you start and push when you finish**, all machines stay in step.

If you forget, `git push` is rejected with a message that the remote contains work you do not
have. Bring that work in, then push again:

```bash
git pull --rebase
git push
```

`--rebase` replays your new commits on top of the ones from the other machine, which keeps the
history linear. If both machines changed the same lines, git stops and marks the conflict in
the file:

```
<<<<<<< HEAD
the version from GitHub
=======
your version
>>>>>>> 1a2b3c4 (your commit message)
```

Edit the file to keep the right version, delete the markers, then:

```bash
git add <file>
git rebase --continue
git push
```

To abandon the attempt and return to where you were, run `git rebase --abort`.

> [!TIP]
> Uncommitted edits block `git pull`. Either commit them first, or set them aside with
> `git stash`, pull, and restore them with `git stash pop`.

## Larger changes

For work that spans several sessions or may not be kept, use a branch so that `main` always
holds a working pipeline:

```bash
git switch -c new-annotator        # create and switch to a branch
# ...commit as usual...
git push -u origin new-annotator   # first push of the branch
```

When the work is finished, open a pull request on GitHub and merge it into `main`, or merge
locally:

```bash
git switch main
git pull
git merge new-annotator
git push
```

## What stays out of the repository

`.gitignore` keeps machine-specific and large files out of the repository. These are never
committed, so back them up separately:

| Path | Contents |
|---|---|
| `config/*` (except `*.template.*`) | Run configs and samplesheets |
| `profiles/*` (except `slurm/`, `local/`, `example/`) | Site profiles, with account names and cluster paths |
| `refs/*` (except the fetch scripts) | Reference genomes, annotations and databases |
| `results*/`, `logs/`, `.snakemake/` | Pipeline outputs and Snakemake state |
| `*.sif`, `*.img` | Container images |

Before committing, check `git status` for anything that should not be public: account names,
cluster paths, email addresses, credentials or data files. Removing a file from a later commit
does not remove it from the history on GitHub.

## Commit messages

Write a short summary line in the imperative, under about 70 characters, describing what the
commit does. Add a blank line and a longer explanation when the reason for the change is not
obvious from the diff.

```
Fit counting flanks per dataset

A fixed 15-nt flank under-subtracted CDS spill-over on high-depth
libraries; the flank is now fitted from reads past canonical start
and stop codons.
```
