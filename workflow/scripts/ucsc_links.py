"""UCSC genome-browser links for ORFs, from BED (0-based half-open) coordinates.

Ported from the panda project's build_master_table.py, generalised past mouse:
any numbered chromosome, X, Y and M/MT map; scaffolds get no link.
"""
from __future__ import annotations

import re
from urllib.parse import quote, unquote, urlencode

UCSC_BASE = "https://genome.ucsc.edu/cgi-bin/hgTracks"
SHORT_SESSION = re.compile(r"^https?://[^/]*genome[^/]*\.ucsc\.edu/s/([^/?#]+)/([^/?#]+)", re.I)


def ucsc_chrom(chrom: str):
    """Seqname -> UCSC name, or None. Accepts 'chr1' or bare '1'.

    Scaffolds are left unmapped: UCSC renames them (GL456233.2 is
    chrUn_GL456233v2 on mm39), which needs the assembly's chromAlias table, and
    a link to a seqname UCSC does not know lands on an error page.
    """
    c = chrom[3:] if chrom.startswith("chr") else chrom
    if re.fullmatch(r"\d+|X|Y", c):
        return "chr" + c
    if c in ("M", "MT"):
        return "chrM"
    return None


def _short_link(user: str, sess: str) -> str:
    return "https://genome.ucsc.edu/s/%s/%s" % (quote(user, safe=""), quote(sess, safe=""))


def parse_session(spec, genome: str):
    """(base_url, extra query params) for a session spec, or for `genome` if none.

    spec: 'user/sessionname' or a short link https://genome.ucsc.edu/s/<user>/<name>.
    UCSC expands the short link server-side and keeps position=. No db= is sent
    with a session: the session defines the assembly, so it must be saved on the
    same build as these coordinates. A live-cart hgsid= URL is refused: it
    expires, and every reader would drive the same cart.
    """
    if not spec:
        return UCSC_BASE, {"db": genome}
    if "://" in spec:
        m = SHORT_SESSION.match(spec)
        if not m:
            raise ValueError("expected https://genome.ucsc.edu/s/<user>/<session>, got %r" % spec)
        return _short_link(unquote(m.group(1)), unquote(m.group(2))), {}
    user, _, sess = spec.partition("/")
    if not user or not sess:
        raise ValueError("expected user/sessionname, got %r" % spec)
    return _short_link(user, sess), {}


def ucsc_link(chrom, start, end, session, pad: int = 50, fmt: str = "hyperlink") -> str:
    """Link for one ORF, or "" when the seqname has no UCSC name.

    fmt="hyperlink": =HYPERLINK("url","chr1:123-456"). Excel and Google Sheets
    show the position and click through; a parser reading the column gets the
    formula text. Comma argument separator (en-US Excel).
    fmt="url": the bare URL.
    """
    uc = ucsc_chrom(str(chrom or ""))
    if uc is None:
        return ""
    try:
        s, e = int(start), int(end)
    except (TypeError, ValueError):
        return ""
    pos = "%s:%d-%d" % (uc, max(1, s + 1 - pad), e + pad)   # BED -> 1-based inclusive
    base, extra = session
    url = base + "?" + urlencode({**extra, "position": pos}, quote_via=quote, safe=":")
    return '=HYPERLINK("%s","%s")' % (url, pos) if fmt == "hyperlink" else url
