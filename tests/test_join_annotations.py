"""join_annotations.py: exclude_types moves ORFs of a type out of Comprehensive and
Primary. Prints PASS/FAIL, exits 1 on failure."""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workflow", "scripts"))
from join_annotations import (GENERAL_OF, exclude_chroms, exclude_no_unshared,  # noqa: E402
                              exclude_types, exclude_unique, general_smorf_type)

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


df = pd.DataFrame(dict(orf_id=["a", "b", "c"], smorf_type=["Trunc-Iso", "N-Iso", "Trunc-Iso"],
                       pass_all=[True, True, False],
                       exclude_reason=["", "", "fewer than 5 in-frame P-sites outside annotated CDS"]))
n = exclude_types(df, ["Trunc-Iso"])
check("two Trunc-Isos moved", n, 2)
check("pass_all cleared", df["pass_all"].tolist(), [False, True, False])
check("reason for Trunc-Iso", df.at[0, "exclude_reason"].startswith("Trunc-Iso: no own sequence"), True)
check("earlier reason kept", df.at[2, "exclude_reason"], "fewer than 5 in-frame P-sites outside annotated CDS")
check("other types untouched", df.at[1, "exclude_reason"], "")
check("no types: no-op", exclude_types(df.copy(), []), 0)
# psORFs need unique-read evidence; other types don't
u = pd.DataFrame(dict(orf_id=["p1", "p2", "l1"], smorf_type=["psORF", "psORF", "lncRNA"],
                      pass_all=[True, True, True], exclude_reason=["", "", ""],
                      unshared_inframe_psites_unique=[0.0, 9.0, 0.0]))
check("unique rule: one psORF moved", exclude_unique(u, ["psORF"], 5), 1)
check("unique rule: pass_all", u["pass_all"].tolist(), [False, True, True])
check("unique rule: reason", u.at[0, "exclude_reason"].startswith("psORF: fewer than 5 in-frame P-sites from uniquely"), True)
check("no unique column: no-op", exclude_unique(u.drop(columns=["unshared_inframe_psites_unique"]), ["psORF"], 5), 0)
# a pseudogene ORF the UniProt relabel renamed is still tested on the Annotator's type
v = pd.DataFrame(dict(orf_id=["p3", "p4"], smorf_type=["annotated_CDS", "N-Iso"],
                      smorf_type_annotator=["psORF", "uORF"], pass_all=[True, True],
                      exclude_reason=["", ""], unshared_inframe_psites_unique=[0.0, 0.0]))
check("unique rule: relabeled psORF moved", exclude_unique(v, ["psORF"], 5), 1)
check("unique rule: reason names the Annotator type", v.at[0, "exclude_reason"].startswith("psORF:"), True)
check("unique rule: non-psORF untouched", v.at[1, "exclude_reason"], "")
# no unshared sequence: excluded unless exactly a UniProt sequence; earlier reason kept
o = pd.DataFrame(dict(orf_id=["c", "a", "t", "n", "k"],
                      smorf_type=["C-Iso", "annotated_CDS", "TrEMBL", "uORF", "Trunc-Iso"],
                      blastoff_class=["sp_partial", "sp_identical", "tr_identical", "novel",
                                      "sp_N_truncated"],
                      unshared_nt=[0, 0, 0, 30, 0], pass_all=[True] * 5,
                      exclude_reason=["", "", "", "", "Trunc-Iso: no own sequence"]))
check("no-unshared: moved", exclude_no_unshared(o), 3)
check("no-unshared: pass_all", o["pass_all"].tolist(), [False, True, False, True, False])
check("no-unshared: tr_identical not exempt", o.at[2, "exclude_reason"], "TrEMBL: no sequence outside annotated CDS")
check("no-unshared: reason", o.at[0, "exclude_reason"], "C-Iso: no sequence outside annotated CDS")
check("no-unshared: sp_identical exempt", o.at[1, "exclude_reason"], "")
check("no-unshared: earlier reason kept", o.at[4, "exclude_reason"], "Trunc-Iso: no own sequence")
check("no-unshared: missing column no-op", exclude_no_unshared(o.drop(columns=["unshared_nt"])), 0)
# trembl_identical_exempt: tr_identical joins sp_identical as exempt
o2 = o.assign(pass_all=True, exclude_reason=["", "", "", "", "Trunc-Iso: no own sequence"])
check("no-unshared exempt tr: moved", exclude_no_unshared(o2, {"sp_identical", "tr_identical"}), 2)
check("no-unshared exempt tr: tr_identical kept", (o2.at[2, "exclude_reason"], o2.at[2, "pass_all"]), ("", True))
# chrM ORFs move to excluded; an earlier reason is kept
m = pd.DataFrame(dict(orf_id=["m1", "m2", "a1"], chrom=["chrM", "chrM", "chr1"],
                      pass_all=[True, False, True], exclude_reason=["", "earlier", ""]))
check("chroms: moved", exclude_chroms(m, ["chrM"]), 2)
check("chroms: pass_all", m["pass_all"].tolist(), [False, False, True])
check("chroms: reason", m["exclude_reason"].tolist(),
      ["chrM: called with the standard genetic code", "earlier", ""])
check("chroms: none given no-op", exclude_chroms(m.copy(), []), 0)
# general_smorf_type: architecture groups; blank for a missing type; unknown label stops
check("general: groups", general_smorf_type(["uORF", "uaoORF", "riORF", "C-Iso", "daoORF", "aORF",
                                              "annotated_CDS", "TrEMBL", float("nan")]),
      ["Upstream", "Upstream", "Internal", "Isoform", "Downstream", "Non-coding", "Reviewed",
       "TrEMBL", ""])
try:
    general_smorf_type(["uORF", "newORF"])
    check("general: unknown label stops", "no error", "SystemExit")
except SystemExit as e:
    check("general: unknown label stops", "newORF" in str(e), True)
# every label the Annotator, the UniProt relabel and annotator_table can emit has a group
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workflow", "scripts", "annotator"))
emitted = {"uORF", "uaORF", "dORF", "daORF", "iORF", "uoORF", "uaoORF", "doORF", "daoORF", "udORF",
           "psORF", "aORF", "riORF", "lncRNA", "eORF", "UA", "intergenic", "annotated_CDS", "TrEMBL",
           "Trunc-Iso", "N-Iso", "C-Iso", "NC-Iso", "Ins-Iso", "Del-Iso", "not_annotated_scaffold",
           "not_annotated_unknown"}
check("general: every emitted label grouped", sorted(emitted - set(GENERAL_OF)), [])

sys.exit(1 if failed else 0)
