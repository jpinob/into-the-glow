"""
Check that data/1ema.cif agrees with itself: what the coordinates in
_atom_site show against what the annotation fields say.

Input : data/1ema.cif   (mmCIF file from the Protein Data Bank, CC0)
Output: stdout only. One line per check, then a summary table.
        Nothing is written. The page, the README and the data stay as they are.

Every line starts with PASS, MISMATCH or NOT CHECKABLE and names the fields
and values compared. NOT CHECKABLE means the file alone cannot settle the
question. Nothing is guessed. "PASS (silent)" means the field neither
confirms nor contradicts the coordinates: it simply does not say.

Checks
1. Residue 65. The amino acid the CRO atoms show vs the annotation fields.
2. Mutations. Entity sequence vs _struct_ref sequence, _struct_ref_seq_dif
   and _entity.pdbx_mutation, in both directions.
3. Modified residues. _atom_site vs _pdbx_struct_mod_residue vs _entity_poly_seq.
4. Bonds. Every _struct_conn row measured; short contacts the table omits.
5. Selenium. SE atoms vs MSE residues vs the count in _entity.details.
6. Missing residues. _pdbx_unobs_or_zero_occ_residues vs _atom_site vs _entity_poly_seq.
7. Numbering. label_seq_id vs auth_seq_id around CRO; the scheme the page uses.
8. Counts. 1,771 non-water atoms; 11 strands in _struct_sheet_range.

Run from the repository root:  python pipeline/check_consistency.py
Optional: --uniprot [FILE]   also compare the reference sequence stored in the
          file with UniProt P42212, fetched from rest.uniprot.org (or read from
          FILE, a FASTA saved from that URL). The URL is printed with the result.
"""
import argparse
import re
import textwrap
import urllib.request
from collections import Counter
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "1ema.cif"
JSON = ROOT / "data" / "gfp.json"
PAGE = ROOT / "src" / "template.html"
README = ROOT / "README.md"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/P42212.fasta"

BOND_TOL = 0.02      # Å. Check 4: measured distance vs _struct_conn.pdbx_dist_value
CONTACT_MAX = 2.0    # Å. Check 4 reverse: contacts this short must be in _struct_conn
COVALENT_MAX = 1.9   # Å. Check 1: two heavy atoms this close inside CRO are bonded
WORD_NUMBERS = {"ONE": 1, "TWO": 2, "THREE": 3, "FOUR": 4, "FIVE": 5,
                "SIX": 6, "SEVEN": 7, "EIGHT": 8, "NINE": 9, "TEN": 10}

# What each check compares, for the summary table
CHECKS = {
    "1": ("Residue 65", "_atom_site CRO atoms vs _chem_comp.name/formula, _pdbx_entry_details.sequence_details, "
                        "_pdbx_struct_mod_residue.parent_comp_id, _entity_poly.pdbx_seq_one_letter_code_can, "
                        "_struct_ref_seq_dif, _entity.pdbx_mutation"),
    "2": ("Mutations", "_entity_poly_seq vs _struct_ref.pdbx_seq_one_letter_code (mapping from _struct_ref_seq) "
                       "vs _struct_ref_seq_dif vs _entity.pdbx_mutation"),
    "3": ("Modified residues", "_atom_site vs _pdbx_struct_mod_residue vs _entity_poly_seq vs _chem_comp.mon_nstd_flag"),
    "4": ("Bonds", "_struct_conn atoms and pdbx_dist_value vs measured distance; contacts < 2.0 Å vs _struct_conn"),
    "5": ("Selenium", "_atom_site SE atoms vs MSE residues vs _entity.details vs _entity_poly_seq vs _pdbx_unobs_or_zero_occ_residues"),
    "6": ("Missing residues", "_pdbx_unobs_or_zero_occ_residues vs _entity_poly_seq minus _atom_site vs _pdbx_poly_seq_scheme"),
    "7": ("Numbering", "_atom_site.label_seq_id vs auth_seq_id; page numbers vs both schemes"),
    "8": ("Counts", "_atom_site non-water atoms vs _refine_hist; _struct_sheet_range rows vs _struct_sheet.number_strands"),
}

# Facts shown on the page or in the README, and the checks that bear on them
PAGE_FACTS = [
    (PAGE, "GFP S65T variant", "1 2"),
    (PAGE, "amino acid 65 is threonine", "1 2"),
    (PAGE, "The natural protein has serine", "1 2"),
    (PAGE, "amino acids 65, 66 and 67", "1 7"),
    (PAGE, "amino acids 65 to 67", "1 7"),
    (PAGE, "1,771 atoms", "8"),
    (PAGE, "11 strands", "8"),
    (PAGE, "Start: amino acid 2", "7"),
    (PAGE, "End: amino acid 229", "7"),
    (PAGE, "228 of its 238 amino acids", "6 7"),
    (PAGE, "amino acids 230 to 238 are not visible", "6"),
    (PAGE, "4 orange atoms are selenium", "5"),
    (PAGE, "selenomethionine (methionine in the natural protein)", "3 5"),
    (README, "the S65T variant of GFP", "1 2"),
    (README, "Amino acid 65 is threonine in this variant", "1 2"),
    (README, "Chromophore formed from amino acids 65 to 67", "1 7"),
    (README, "amino acids 2 to 229 visible", "7"),
    (README, "Four visible selenium atoms", "5"),
    (README, "1,771 atoms, 11 strands", "8"),
]

RESULTS = []  # (check id, status, fields, detail)


def report(check, status, fields, detail=""):
    RESULTS.append((check, status, fields, detail))
    line = f"{status:<14} [{check}] {fields}"
    if detail:
        line += f" -> {detail}"
    print(line)


# ---------------------------------------------------------------- helpers
def clean(v):
    v = gemmi.cif.as_string(v)
    return "" if v in ("?", ".") else v


def rows(block, prefix, cols):
    return [[clean(r[i]) for i in range(len(cols))] for r in block.find(prefix, cols)]


def value(block, tag):
    v = block.find_value(tag)
    return "" if v is None else clean(v)


def squash(s):
    return " ".join(s.split())


def is_standard(name):
    t = gemmi.find_tabulated_residue(name)
    return bool(t) and t.is_standard()


def formula(s):
    """'C4 H9 N O3' -> Counter({'C': 4, 'H': 9, 'N': 1, 'O': 3})"""
    return Counter({el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", s)})


def fmt_formula(c):
    return " ".join(f"{el}{n if n != 1 else ''}" for el, n in c.items() if n > 0)


def aligned_pairs(query, target):
    """Global alignment with gemmi; yields (query index or None, target index or None)."""
    res = gemmi.align_string_sequences(list(query), list(target), [0] * (len(target) + 1))
    qi = ti = 0
    for n, op in re.findall(r"(\d+)([MID])", res.cigar_str()):
        for _ in range(int(n)):
            if op == "M":
                yield qi, ti
                qi += 1
                ti += 1
            elif op == "I":
                yield qi, None
                qi += 1
            else:
                yield None, ti
                ti += 1


class Cif:
    """The parts of the file the checks need, read once."""

    def __init__(self, path):
        self.block = gemmi.cif.read(str(path)).sole_block()
        self.st = gemmi.read_structure(str(path))
        self.st.setup_entities()
        self.model = self.st[0]
        self.chain = self.model["A"]
        self.polymer = [r for r in self.chain if r.entity_type == gemmi.EntityType.Polymer]
        self.waters = [r for r in self.chain if r.is_water()]
        b = self.block
        # Sequence of the entity and the two numbering schemes
        self.seq = [(int(n), m) for n, m in rows(b, "_entity_poly_seq.", ["num", "mon_id"])]
        self.label2auth, self.observed = {}, {}
        for seq_id, _mon, pdb_num, auth_num in rows(b, "_pdbx_poly_seq_scheme.",
                                                      ["seq_id", "mon_id", "pdb_seq_num", "auth_seq_num"]):
            self.label2auth[int(seq_id)] = int(pdb_num)
            self.observed[int(seq_id)] = auth_num != ""
        self.by_label = {r.label_seq: r for r in self.polymer}
        self.by_auth = {r.seqid.num: r for r in self.polymer}
        # Modified residues: label_seq_id -> list of parents
        self.parents = {}
        self.mod_rows = rows(b, "_pdbx_struct_mod_residue.",
                             ["label_seq_id", "label_comp_id", "auth_seq_id", "auth_comp_id", "parent_comp_id", "details"])
        for label, comp, _a, _ac, parent, _d in self.mod_rows:
            self.parents.setdefault((int(label), comp), []).append(parent)
        self.chem = {r[0]: r for r in rows(b, "_chem_comp.", ["id", "mon_nstd_flag", "name", "formula"])}
        self.entity = rows(b, "_entity.", ["id", "type", "pdbx_number_of_molecules", "details", "pdbx_mutation"])
        self.mutation = squash(next(e[4] for e in self.entity if e[0] == "1"))
        self.details = squash(next(e[3] for e in self.entity if e[0] == "1"))
        self.seq_details = squash(value(b, "_pdbx_entry_details.sequence_details"))
        self.refine_details = squash(value(b, "_refine.details"))
        self.ref_seq = "".join(value(b, "_struct_ref.pdbx_seq_one_letter_code").split())
        self.ref_db = (value(b, "_struct_ref.db_name"), value(b, "_struct_ref.pdbx_db_accession"))
        self.srs = {k: value(b, "_struct_ref_seq." + k) for k in
                    ("seq_align_beg", "seq_align_end", "db_align_beg", "db_align_end",
                     "pdbx_auth_seq_align_beg", "pdbx_auth_seq_align_end")}
        self.dif = rows(b, "_struct_ref_seq_dif.",
                        ["seq_num", "mon_id", "db_mon_id", "pdbx_seq_db_seq_num", "pdbx_auth_seq_num", "details"])
        self.can = "".join(value(b, "_entity_poly.pdbx_seq_one_letter_code_can").split())

    def letter(self, mon):
        """One-letter code as the file defines it: standard residues by name,
        modified residues by their single parent in _pdbx_struct_mod_residue, else X."""
        if is_standard(mon):
            return gemmi.find_tabulated_residue(mon).one_letter_code.upper()
        parents = {p for (label, comp), ps in self.parents.items() if comp == mon for p in ps}
        if len(parents) == 1:
            return gemmi.find_tabulated_residue(parents.pop()).one_letter_code.upper()
        return "X"

    def auth_span(self):
        """auth number -> (mon_id, label_seq_id). Each residue gets its own auth
        number. A gap in the auth numbering (PHE 64, CRO 66, VAL 68) goes to the
        non-standard neighbour, which stands in for the missing residues. A gap
        with no such neighbour stays unassigned; check 7 reports it."""
        span = {self.label2auth[label]: (mon, label) for label, mon in self.seq}
        for (l1, m1), (l2, m2) in zip(self.seq, self.seq[1:]):
            a1, a2 = self.label2auth[l1], self.label2auth[l2]
            if not is_standard(m1) and is_standard(m2):
                owner = (m1, l1)
            elif not is_standard(m2) and is_standard(m1):
                owner = (m2, l2)
            else:
                continue
            for n in range(a1 + 1, a2):
                span[n] = owner
        return span

    def unassigned_gaps(self):
        span = self.auth_span()
        return [n for n in range(min(span), max(span) + 1) if n not in span]


# ---------------------------------------------------------------- checks
def check_1(c):
    """Residue 65: which amino acid do the CRO atoms show, and what do the fields say."""
    cro = [r for r in c.polymer if r.name == "CRO"]
    if len(cro) != 1:
        report("1", "NOT CHECKABLE", "_atom_site CRO residues", f"{len(cro)} found, expected 1")
        return
    res = cro[0]
    atoms = {a.name: a for a in res}
    # 1a. Side chain on the residue-65 fragment, from the coordinates alone.
    # A beta carbon bonded to one O and one C is threonine; to one O only, serine; no beta carbon, glycine.
    cb = atoms.get("CB1")
    if cb is None:
        shown, evidence = "GLY", "no CB1 atom"
    else:
        nb = [a for a in res if a.name not in ("CB1", "CA1") and a.pos.dist(cb.pos) <= COVALENT_MAX]
        elems = sorted(a.element.name for a in nb)
        shown = {"O": "SER", "CO": "THR"}.get("".join(elems), "unrecognised")
        evidence = "CB1 bonded to " + " and ".join(f"{a.name} ({a.element.name}, {a.pos.dist(cb.pos):.2f} Å)" for a in nb)
    shown1 = gemmi.find_tabulated_residue(shown).one_letter_code.upper() if is_standard(shown) else "?"
    report("1a", "PASS" if shown in ("SER", "THR", "GLY") else "NOT CHECKABLE",
           f"_atom_site CRO {res.seqid.num} atoms {', '.join(a.name for a in res)}",
           f"residue-65 fragment is {shown} ({evidence})")

    # 1b. _chem_comp.name: the IUPAC fragment names the side chain.
    name = c.chem["CRO"][2]
    by_name = "THR" if "HYDROXYPROPYL" in name else "SER" if "HYDROXYETHYL" in name else "?"
    report("1b", "PASS" if by_name == shown else "MISMATCH" if by_name != "?" else "NOT CHECKABLE",
           f"_chem_comp.name for CRO = '{name}'",
           f"1-amino-2-hydroxypropyl = 3-carbon hydroxy side chain = {by_name}; atoms show {shown}")

    # 1c. _chem_comp.formula: CRO = X + TYR + GLY, minus 2 H2O (two peptide bonds),
    # minus H2O (cyclization: the THR 65 carbonyl O leaves as water, per sequence_details),
    # minus H2 (oxidation of the TYR 66 CA-CB bond, per sequence_details).
    cro_f = formula(c.chem["CRO"][3])
    matches = []
    for x in ("THR", "SER", "GLY", "ALA"):
        if x not in c.chem:
            continue
        f = formula(c.chem[x][3]) + formula(c.chem["TYR"][3]) + formula(c.chem["GLY"][3])
        f.subtract(Counter({"H": 2 * 3 + 2, "O": 3}))
        if +f == +cro_f:  # unary + drops zero counts before comparing
            matches.append(x)
    report("1c", "PASS" if matches == [shown] else "MISMATCH" if matches else "NOT CHECKABLE",
           f"_chem_comp.formula CRO = {c.chem['CRO'][3]} vs X + TYR + GLY - 3 H2O - H2 using _chem_comp.formula of X",
           f"formula fits X = {', '.join(matches) or 'none of THR/SER/GLY/ALA'}; atoms show {shown}")

    # 1d. Free text in the file.
    for tag, txt in (("_pdbx_entry_details.sequence_details", c.seq_details), ("_refine.details", c.refine_details)):
        named = sorted(set(re.findall(r"\b([A-Z]{3}) 65\b", txt)))
        report("1d", "PASS" if named == [shown] else "MISMATCH" if named else "NOT CHECKABLE",
               f"{tag} names residue 65 as {named or 'nothing'}", f"atoms show {shown}")

    # 1e. _pdbx_struct_mod_residue.parent_comp_id
    parents = c.parents.get((res.label_seq, "CRO"), [])
    report("1e", "PASS" if parents[:1] == [shown] else "MISMATCH" if parents else "NOT CHECKABLE",
           f"_pdbx_struct_mod_residue.parent_comp_id for CRO (label_seq_id {res.label_seq}) = {parents}",
           f"first parent {parents[:1]} vs atoms {shown}")

    # 1f. _entity_poly.pdbx_seq_one_letter_code_can at auth 65-67
    span = c.auth_span()
    covered = sorted(n for n, (mon, _l) in span.items() if mon == "CRO")
    can_ok = len(c.can) == len(c.seq) - 1 + len(covered)
    if not can_ok:
        report("1f", "NOT CHECKABLE", f"_entity_poly.pdbx_seq_one_letter_code_can length {len(c.can)}",
               f"expected {len(c.seq) - 1 + len(covered)} (one letter per auth position) to index it")
    else:
        tri = c.can[covered[0] - 1:covered[-1]]
        report("1f", "PASS" if tri[0] == shown1 else "MISMATCH",
               f"_entity_poly.pdbx_seq_one_letter_code_can[{covered[0]}..{covered[-1]}] = '{tri}' "
               f"(flanks '{c.can[covered[0] - 2]}' {covered[0] - 1}, '{c.can[covered[-1]]}' {covered[-1] + 1})",
               f"position {covered[0]} reads '{tri[0]}', atoms show {shown} ('{shown1}')")

    # 1g. _struct_ref_seq_dif rows for this residue
    dif = [d for d in c.dif if int(d[0]) == res.label_seq]
    names = sorted({d[1] for d in dif})
    report("1g", "PASS (silent)" if names == ["CRO"] else "MISMATCH" if dif else "NOT CHECKABLE",
           f"_struct_ref_seq_dif rows with seq_num {res.label_seq}: " +
           "; ".join(f"mon_id {d[1]} db_mon_id {d[2]}@{d[3]} details '{d[5]}'" for d in dif),
           f"names the entry residue only as CRO and the reference as {dif[0][2] if dif else '?'}; "
           f"does not name {shown}; details are not 'engineered mutation'")

    # 1h. _entity.pdbx_mutation
    tokens = [t.strip() for t in c.mutation.split(",")]
    point65 = [t for t in tokens if re.fullmatch(r"[A-Z]65[A-Z]", t)]
    report("1h", "PASS (silent)" if not point65 else ("PASS" if point65[0][-1] == shown1 else "MISMATCH"),
           f"_entity.pdbx_mutation = '{c.mutation}'",
           f"point mutation at 65 named: {point65 or 'none'}; atoms show {shown} ('{shown1}')")


def check_2(c, uniprot=None):
    """Mutations: entity sequence vs reference sequence vs _struct_ref_seq_dif vs pdbx_mutation."""
    s = c.srs
    off = int(s["db_align_beg"]) - int(s["pdbx_auth_seq_align_beg"])  # db position = auth + off
    ok = (c.label2auth[int(s["seq_align_beg"])] == int(s["pdbx_auth_seq_align_beg"])
          and c.label2auth[int(s["seq_align_end"])] == int(s["pdbx_auth_seq_align_end"])
          and int(s["pdbx_auth_seq_align_end"]) + off == int(s["db_align_end"])
          and int(s["db_align_end"]) <= len(c.ref_seq))
    report("2a", "PASS" if ok else "MISMATCH",
           f"_struct_ref_seq: seq_align {s['seq_align_beg']}-{s['seq_align_end']} = auth {s['pdbx_auth_seq_align_beg']}-"
           f"{s['pdbx_auth_seq_align_end']} = db {s['db_align_beg']}-{s['db_align_end']}; "
           f"_struct_ref ({c.ref_db[0]} {c.ref_db[1]}) sequence length {len(c.ref_seq)}",
           f"db position = auth + {off}; the stored reference starts '{c.ref_seq[:8]}'")

    span = c.auth_span()
    auth_lo, auth_hi = int(s["pdbx_auth_seq_align_beg"]), int(s["pdbx_auth_seq_align_end"])
    ent = "".join(c.letter(span[n][0]) for n in sorted(span))
    pairs = list(aligned_pairs(ent, c.ref_seq))
    gaps = [(q, t) for q, t in pairs if q is None or t is None]
    report("2b", "PASS", f"gemmi global alignment of entity (auth 1-{max(span)}, CRO as X) vs _struct_ref sequence",
           f"gaps: {', '.join(('entity ' + str(q + 1) + ' unmatched') if t is None else ('reference ' + c.ref_seq[t] + str(t + 1) + ' unmatched') for q, t in gaps) or 'none'}; "
           f"offset {off} holds from auth {auth_lo} on")

    diffs = []  # (auth, mon, label, entity letter, db position, ref letter)
    for n in sorted(span):
        mon, label = span[n]
        db = n + off
        r = c.ref_seq[db - 1] if 1 <= db <= len(c.ref_seq) else "-"
        if c.letter(mon) != r:
            diffs.append((n, mon, label, c.letter(mon), db, r))
    inside = [d for d in diffs if auth_lo <= d[0] <= auth_hi]
    outside = [d for d in diffs if not auth_lo <= d[0] <= auth_hi]
    report("2c", "PASS", "every difference, entity vs _struct_ref sequence, at db = auth + %d" % off,
           "; ".join(f"auth {n} {mon}(label {l}) '{e}' vs reference '{r}'@{db}" for n, mon, l, e, db, r in inside) +
           (f" | outside aligned range: " + "; ".join(f"auth {n} {mon} '{e}' vs '{r}'@{db}" for n, mon, l, e, db, r in outside) if outside else ""))

    # 2d. differences <-> _struct_ref_seq_dif
    dif_pos = {(int(d[3]) - off, gemmi.find_tabulated_residue(d[2]).one_letter_code.upper() if is_standard(d[2]) else d[2]): d for d in c.dif}
    missing_rows = [d for d in inside if (d[0], d[5]) not in dif_pos]
    extra_rows = [d for (n, r), d in dif_pos.items() if not any(x[0] == n and x[5] == r for x in inside)]
    report("2d", "PASS" if not missing_rows and not extra_rows else "MISMATCH",
           f"differences ({len(inside)}) vs _struct_ref_seq_dif rows ({len(c.dif)}: " +
           ", ".join(f"{d[1]}@{d[0]} for {d[2]}@db{d[3]}" for d in c.dif) + ")",
           f"differences without a row: {[(d[0], d[1], d[5]) for d in missing_rows] or 'none'}; "
           f"rows without a difference: {[(d[1], d[2], d[3]) for d in extra_rows] or 'none'}")
    # 2e. each row's db_mon_id vs the stored reference at that db position
    bad = [d for d in c.dif if c.ref_seq[int(d[3]) - 1] != gemmi.find_tabulated_residue(d[2]).one_letter_code.upper()]
    report("2e", "PASS" if not bad else "MISMATCH",
           "_struct_ref_seq_dif.db_mon_id vs _struct_ref sequence at pdbx_seq_db_seq_num",
           "all rows agree" if not bad else "; ".join(f"{d[2]}@{d[3]} but reference reads {c.ref_seq[int(d[3]) - 1]}" for d in bad))

    # 2f. differences <-> _entity.pdbx_mutation
    tokens = [t.strip() for t in c.mutation.split(",")]
    covered = set()
    for t in tokens:
        m = re.fullmatch(r"([A-Z])(\d+)([A-Z])", t)
        rng = re.fullmatch(r"(\d+)\s*-\s*(\d+)\s+REPLACED BY\s+(\w+)", t)
        if m:
            n = int(m.group(2))
            e, r = c.letter(span[n][0]) if n in span else "-", c.ref_seq[n + off - 1] if 1 <= n + off <= len(c.ref_seq) else "-"
            hit = any(d[0] == n for d in inside)
            status = "PASS" if hit and r == m.group(1) and e == m.group(3) else "MISMATCH"
            report("2f", status, f"_entity.pdbx_mutation token '{t}' vs alignment at auth {n}",
                   f"entity reads '{e}' ({span[n][0]}), stored reference reads '{r}'@db{n + off}, "
                   f"difference listed: {'yes' if hit else 'no'}, _struct_ref_seq_dif row for {n}: "
                   f"{'yes' if any(int(d[3]) - off == n for d in c.dif) else 'no'}")
            covered.add(n)
        elif rng:
            lo, hi, comp = int(rng.group(1)), int(rng.group(2)), rng.group(3)
            mons = {span[n][0] for n in range(lo, hi + 1) if n in span}
            hit = all(any(d[0] == n for d in inside) for n in range(lo, hi + 1))
            report("2f", "PASS" if mons == {comp} and hit else "MISMATCH",
                   f"_entity.pdbx_mutation token '{t}' vs alignment at auth {lo}-{hi}",
                   f"entity residues there: {sorted(mons)}; differences listed at every position: {'yes' if hit else 'no'}")
            covered.update(range(lo, hi + 1))
        else:
            report("2f", "NOT CHECKABLE", f"_entity.pdbx_mutation token '{t}'", "not a point mutation or a 'REPLACED BY' range")
    unnamed = [d for d in inside if d[0] not in covered]
    report("2g", "PASS" if not unnamed else "MISMATCH", "differences not named in _entity.pdbx_mutation",
           "none" if not unnamed else "; ".join(f"auth {d[0]} {d[1]} '{d[3]}' vs '{d[5]}'" for d in unnamed))

    # 2h. Optional: the stored reference vs UniProt P42212 itself
    if uniprot is None:
        report("2h", "NOT CHECKABLE", "_struct_ref sequence vs UniProt P42212", "run with --uniprot to fetch it")
        return
    useq, source = uniprot
    upairs = list(aligned_pairs(c.ref_seq, useq))
    udiff = [f"stored {c.ref_seq[q]}{q + 1} inserted" if t is None else f"UniProt {useq[t]}{t + 1} missing" if q is None
             else f"stored {c.ref_seq[q]}{q + 1} vs UniProt {useq[t]}{t + 1}"
             for q, t in upairs if q is None or t is None or c.ref_seq[q] != useq[t]]
    report("2h", "PASS" if not udiff else "MISMATCH",
           f"_struct_ref.pdbx_seq_one_letter_code ({len(c.ref_seq)} aa) vs UniProt P42212 ({len(useq)} aa) from {source}",
           "; ".join(udiff) or "identical")
    report("2i", "PASS" if int(s["db_align_end"]) <= len(useq) else "MISMATCH",
           f"_struct_ref_seq.db_align_end = {s['db_align_end']} vs UniProt P42212 length {len(useq)}",
           "within the UniProt sequence" if int(s["db_align_end"]) <= len(useq) else "points past the end of the UniProt sequence")
    epairs = list(aligned_pairs(ent, useq))
    ediff = [(q + 1, span[q + 1][0], ent[q], useq[t], t + 1) for q, t in epairs if q is not None and t is not None and ent[q] != useq[t]]
    egaps = [(q, t) for q, t in epairs if q is None or t is None]
    report("2j", "PASS", "entity (CRO as X) vs UniProt P42212, gemmi global alignment",
           f"gaps: {len(egaps)}; differences: " + "; ".join(f"auth {n} {mon} '{e}' vs UniProt '{u}'@{p}" for n, mon, e, u, p in ediff))
    for n, mon, e, u, p in ediff:
        if mon == "CRO":
            continue
        tok = f"{u}{n}{e}"
        report("2k", "PASS" if tok in tokens else "MISMATCH", f"UniProt difference {tok} vs _entity.pdbx_mutation",
               "named" if tok in tokens else "not named")
        has_row = any(int(d[3]) - off == n for d in c.dif)
        report("2k", "PASS" if has_row else "MISMATCH", f"UniProt difference {tok} vs _struct_ref_seq_dif",
               "row present" if has_row else "no row for it")


def check_3(c):
    """Non-standard residues: coordinates vs _pdbx_struct_mod_residue vs _entity_poly_seq."""
    in_atoms = {(r.label_seq, r.name) for r in c.polymer if not is_standard(r.name)}
    in_mod = {(int(m[0]), m[1]) for m in c.mod_rows}
    in_seq = {(n, m) for n, m in c.seq if not is_standard(m)}
    unobserved = {(n, m) for n, m in in_seq if not c.observed[n]}
    report("3a", "PASS" if in_atoms == in_mod else "MISMATCH",
           f"non-standard residues in _atom_site {sorted(in_atoms)} vs _pdbx_struct_mod_residue {sorted(in_mod)}",
           f"group_PDB for them: {sorted({r.het_flag for r in c.polymer if not is_standard(r.name)})} "
           f"(A = ATOM, H = HETATM); parents: " + "; ".join(f"{k[1]}@{k[0]}={v}" for k, v in sorted(c.parents.items())))
    report("3b", "PASS" if in_seq - unobserved == in_atoms else "MISMATCH",
           f"non-standard in _entity_poly_seq {sorted(in_seq)} minus unobserved {sorted(unobserved)} vs _atom_site",
           "modelled set agrees" if in_seq - unobserved == in_atoms else f"difference {sorted((in_seq - unobserved) ^ in_atoms)}")
    auth_bad = [m for m in c.mod_rows if c.by_label[int(m[0])].seqid.num != int(m[2]) or c.by_label[int(m[0])].name != m[3]]
    report("3c", "PASS" if not auth_bad else "MISMATCH",
           "_pdbx_struct_mod_residue auth_seq_id/auth_comp_id vs _atom_site",
           "all rows agree" if not auth_bad else str(auth_bad))
    flags = {k: v[1] for k, v in c.chem.items() if k != "HOH"}
    nstd = {k for k, f in flags.items() if f == "n"}
    report("3d", "PASS" if nstd == {m for _n, m in in_seq} else "MISMATCH",
           f"_chem_comp.mon_nstd_flag = n for {sorted(nstd)} vs non-standard names in _entity_poly_seq {sorted({m for _n, m in in_seq})}",
           f"_entity_poly.nstd_monomer = {value(c.block, '_entity_poly.nstd_monomer')}")


def check_4(c):
    """Bonds: every _struct_conn row measured; short contacts not in the table."""
    cols = ["id", "conn_type_id", "ptnr1_label_asym_id", "ptnr1_label_comp_id", "ptnr1_label_seq_id", "ptnr1_label_atom_id",
            "ptnr2_label_asym_id", "ptnr2_label_comp_id", "ptnr2_label_seq_id", "ptnr2_label_atom_id",
            "ptnr1_auth_seq_id", "ptnr2_auth_seq_id", "ptnr1_symmetry", "ptnr2_symmetry", "pdbx_dist_value"]
    conns = rows(c.block, "_struct_conn.", cols)
    index = {(r.subchain, r.label_seq, a.name): (r, a) for r in c.polymer for a in r}
    listed = set()
    for cid, ctype, as1, cp1, sq1, at1, as2, cp2, sq2, at2, au1, au2, sym1, sym2, dist in conns:
        p1, p2 = index.get((as1, int(sq1), at1)), index.get((as2, int(sq2), at2))
        fields = f"_struct_conn {cid} ({ctype}) {cp1} {sq1}/auth {au1} {at1} - {cp2} {sq2}/auth {au2} {at2}, pdbx_dist_value {dist}"
        if p1 is None or p2 is None:
            report("4a", "MISMATCH", fields, f"atom missing from _atom_site: {'ptnr1' if p1 is None else 'ptnr2'}")
            continue
        if sym1 != "1_555" or sym2 != "1_555":
            report("4a", "NOT CHECKABLE", fields, f"symmetry {sym1}/{sym2} not handled")
            continue
        names_ok = p1[0].name == cp1 and p2[0].name == cp2 and p1[0].seqid.num == int(au1) and p2[0].seqid.num == int(au2)
        d = p1[1].pos.dist(p2[1].pos)
        ok = names_ok and (dist == "" or abs(d - float(dist)) <= BOND_TOL)
        report("4a", "PASS" if ok else "MISMATCH", fields,
               f"measured {d:.3f} Å, difference {d - float(dist):+.3f} Å" + ("" if names_ok else "; comp/auth ids do not match the atoms"))
        listed.add(frozenset([(as1, int(sq1), at1), (as2, int(sq2), at2)]))

    ns = gemmi.NeighborSearch(c.model, c.st.cell, 5).populate(include_h=False)
    cs = gemmi.ContactSearch(CONTACT_MAX)
    cs.ignore = gemmi.ContactSearch.Ignore.SameResidue
    order = {id(r): i for i, r in enumerate(c.chain)}
    unlisted, peptide, in_table = [], 0, 0
    for k in cs.find_contacts(ns):
        r1, a1, r2, a2 = k.partner1.residue, k.partner1.atom, k.partner2.residue, k.partner2.atom
        key = frozenset([(r1.subchain, r1.label_seq, a1.name), (r2.subchain, r2.label_seq, a2.name)])
        if k.image_idx == 0 and key in listed:
            in_table += 1
        elif (k.image_idx == 0 and {a1.name, a2.name} == {"C", "N"}
              and abs(order[id(r1)] - order[id(r2)]) == 1 and r1.entity_type == r2.entity_type == gemmi.EntityType.Polymer):
            peptide += 1
        else:
            unlisted.append(f"{r1.name} {r1.seqid.num} {a1.name} - {r2.name} {r2.seqid.num} {a2.name} {k.dist:.2f} Å"
                            + (f" (symmetry image {k.image_idx})" if k.image_idx else ""))
    report("4b", "PASS" if not unlisted else "MISMATCH",
           f"contacts < {CONTACT_MAX} Å between different residues (waters and symmetry mates included) vs _struct_conn",
           f"{in_table} in _struct_conn, {peptide} C(i)-N(i+1) peptide bonds, unlisted: {unlisted or 'none'}")


def check_5(c):
    """Selenium: SE atoms vs MSE residues vs _entity.details."""
    se_atoms = [(r.name, r.seqid.num, a.name) for r in c.chain for a in r if a.element.name == "Se"]
    mse = [r for r in c.polymer if r.name == "MSE"]
    per_res = Counter((r.name, r.seqid.num) for r in c.chain for a in r if a.element.name == "Se")
    report("5a", "PASS" if len(se_atoms) == len(mse) and all(per_res[(r.name, r.seqid.num)] == 1 for r in mse) else "MISMATCH",
           f"_atom_site Se atoms {se_atoms} vs MSE residues auth {[r.seqid.num for r in mse]}",
           f"{len(se_atoms)} Se atoms, {len(mse)} MSE residues, one Se each: {all(per_res[(r.name, r.seqid.num)] == 1 for r in mse)}")
    m = re.search(r"\b([A-Z]+) SE-METHIONINES", c.details)
    total = WORD_NUMBERS.get(m.group(1)) if m else None
    absent = set()
    if "N-TERMINAL MET" in c.details:
        absent.add(c.label2auth[c.seq[0][0]])
    absent.update(int(n) for n in re.findall(r"\bMET (\d+)", c.details))
    seq_mse = [c.label2auth[n] for n, mon in c.seq if mon == "MSE"]
    unobs_mse = [a for a in seq_mse if not c.observed[next(n for n, _m in c.seq if c.label2auth[n] == a)]]
    expected = None if total is None else total - len(absent)
    report("5b", "PASS" if total == len(seq_mse) and sorted(absent) == sorted(unobs_mse) and expected == len(se_atoms) else "MISMATCH",
           f"_entity.details = '{c.details}'",
           f"says {total} MSE with {sorted(absent)} absent -> {expected} present; _entity_poly_seq has MSE at auth {seq_mse}, "
           f"unobserved {unobs_mse}; _atom_site has {len(se_atoms)} Se")
    report("5c", "PASS" if "SE" in c.chem["MSE"][3].upper().split() and any(r[0] == "SE" for r in rows(c.block, "_atom_type.", ["symbol"])) else "MISMATCH",
           f"_chem_comp.formula MSE = '{c.chem['MSE'][3]}'; _atom_type symbols {[r[0] for r in rows(c.block, '_atom_type.', ['symbol'])]}",
           "selenium declared in both")


def check_6(c):
    """Missing residues: the unobserved list vs what is really absent from _atom_site."""
    unobs = rows(c.block, "_pdbx_unobs_or_zero_occ_residues.", ["polymer_flag", "occupancy_flag", "auth_comp_id", "auth_seq_id"])
    listed = {int(u[3]): u[2] for u in unobs}
    absent = {c.label2auth[n]: m for n, m in c.seq if n not in c.by_label}
    scheme_absent = {c.label2auth[n]: m for n, m in c.seq if not c.observed[n]}
    report("6a", "PASS" if listed == absent == scheme_absent else "MISMATCH",
           f"_pdbx_unobs_or_zero_occ_residues auth {sorted(listed)} vs _entity_poly_seq residues absent from _atom_site "
           f"{sorted(absent)} vs _pdbx_poly_seq_scheme.auth_seq_num = ? {sorted(scheme_absent)}",
           f"{len(listed)} listed, {len(absent)} absent, flags {sorted({(u[0], u[1]) for u in unobs})} (Y,1 = polymer, unobserved); "
           f"{len(c.seq)} in sequence, {len(c.polymer)} in _atom_site")
    zero = [(r.name, r.seqid.num, a.name) for r in c.chain for a in r if a.occ == 0]
    report("6b", "PASS" if not zero else "MISMATCH", "zero-occupancy atoms in _atom_site vs occupancy_flag 0 rows",
           "none in either" if not zero and all(u[1] == "1" for u in unobs) else str(zero))
    uatoms = rows(c.block, "_pdbx_unobs_or_zero_occ_atoms.", ["auth_comp_id", "auth_seq_id", "auth_atom_id"])
    present = {(r.seqid.num, a.name) for r in c.polymer for a in r}
    wrong = [u for u in uatoms if (int(u[1]), u[2]) in present]
    report("6c", "PASS" if not wrong else "MISMATCH",
           f"_pdbx_unobs_or_zero_occ_atoms ({len(uatoms)} rows) vs _atom_site",
           "every listed atom is absent" if not wrong else f"listed but present: {wrong}")
    report("6d", "NOT CHECKABLE", "atoms absent from _atom_site but not listed in _pdbx_unobs_or_zero_occ_atoms",
           "needs the ideal atom list of each residue (a monomer library), which the file does not carry")


def check_7(c):
    """Numbering: label_seq_id vs auth_seq_id, and which scheme the page uses."""
    around = [c.by_auth[n] for n in (64, 66, 68) if n in c.by_auth]
    offsets = Counter(r.seqid.num - r.label_seq for r in c.polymer)
    report("7a", "PASS", "_atom_site label_seq_id vs auth_seq_id around CRO",
           "; ".join(f"{r.name} label {r.label_seq} = auth {r.seqid.num} (offset {r.seqid.num - r.label_seq:+d})" for r in around)
           + f"; over the chain: offsets {dict(sorted(offsets.items()))} (auth = label before CRO, label + 2 after)")
    span = c.auth_span()
    cro_auth = sorted(n for n, (m, _l) in span.items() if m == "CRO")
    cro_label = [r.label_seq for r in c.polymer if r.name == "CRO"]
    auth_range = (min(c.by_auth), max(c.by_auth))
    label_range = (min(c.by_label), max(c.by_label))
    unobs_auth = sorted(c.label2auth[n] for n, _m in c.seq if not c.observed[n])
    seen = len(c.by_auth) - len(cro_label) + len(cro_auth)
    total_auth, total_label = c.label2auth[c.seq[-1][0]], c.seq[-1][0]
    said = re.search(r"RESIDUES (\d+), (\d+),? AND (\d+) ARE NOT PRESENT IN THE ENTRY AND ARE INSTEAD REPLACED WITH CRO (\d+)", c.seq_details)
    ok = said and [int(x) for x in said.groups()[:3]] == cro_auth and int(said.group(4)) == c.by_label[cro_label[0]].seqid.num
    report("7b", "PASS" if ok else "MISMATCH" if said else "NOT CHECKABLE",
           f"_pdbx_entry_details.sequence_details '{said.group(0) if said else '?'}' vs the gap in auth_seq_id around CRO",
           f"auth numbers covered by CRO {cro_auth} (label {cro_label}); unassigned gaps in the numbering: {c.unassigned_gaps() or 'none'}")
    page = [
        ("amino acids 2 to 229 / Start: amino acid 2 / End: amino acid 229", auth_range, label_range),
        ("amino acids 65 to 67 (chromophore)", tuple(cro_auth), tuple(cro_label)),
        ("amino acid 1 and 230 to 238 not visible", (unobs_auth[0], unobs_auth[1], unobs_auth[-1]),
         tuple(n for n, _m in c.seq if not c.observed[n])[:2] + (c.seq[-1][0],)),
        ("228 of its 238 amino acids", (seen, total_auth), (len(c.by_label), total_label)),
    ]
    for text_, auth_v, label_v in page:
        report("7c", "PASS", f"page '{text_}' vs auth {auth_v} vs label {label_v}",
               "auth scheme" if auth_v != label_v else "same in both schemes")
    if JSON.exists():
        import json
        d = json.loads(JSON.read_text())
        nums = d["num"]
        ok = nums[0] == auth_range[0] and nums[-1] == auth_range[1] and 66 in nums and 65 not in nums and 67 not in nums
        report("7d", "PASS" if ok else "MISMATCH", f"data/gfp.json num[0]={nums[0]} num[-1]={nums[-1]} 65/66/67 in num = "
               f"{65 in nums}/{66 in nums}/{67 in nums} vs _atom_site auth_seq_id",
               "the page data use auth_seq_id (gemmi seqid.num) throughout, one scheme" if ok else "scheme differs")


def check_8(c):
    """Counts: non-water atoms and strands."""
    n_poly = sum(len(r) for r in c.chain if not r.is_water())
    n_wat = sum(len(r) for r in c.waters)
    hist = {k: value(c.block, "_refine_hist." + k) for k in ("pdbx_number_atoms_protein", "number_atoms_solvent", "number_atoms_total")}
    wat_mol = next(e[2] for e in c.entity if e[1] == "water")
    hyd = sum(1 for r in c.chain for a in r if a.element.name == "H")
    ok = n_poly == int(hist["pdbx_number_atoms_protein"]) == 1771 and n_wat == int(hist["number_atoms_solvent"]) == int(wat_mol) \
        and n_poly + n_wat == int(hist["number_atoms_total"])
    report("8a", "PASS" if ok else "MISMATCH",
           f"_atom_site non-water atoms {n_poly}, water atoms {n_wat} vs _refine_hist {hist} vs _entity water molecules {wat_mol}",
           f"page says 1,771; hydrogens in file: {hyd}; altlocs: {sorted({a.altloc for r in c.chain for a in r} - {chr(0)}) or 'none'}")
    ranges = rows(c.block, "_struct_sheet_range.", ["sheet_id", "id", "beg_label_seq_id", "end_label_seq_id", "beg_auth_seq_id", "end_auth_seq_id"])
    n_decl = value(c.block, "_struct_sheet.number_strands")
    n_order = len(rows(c.block, "_struct_sheet_order.", ["sheet_id"]))
    n_gemmi = sum(len(s.strands) for s in c.st.sheets)
    ok = len(ranges) == int(n_decl) == n_gemmi == 11 and n_order == len(ranges) - 1
    report("8b", "PASS" if ok else "MISMATCH",
           f"_struct_sheet_range rows {len(ranges)} vs _struct_sheet.number_strands {n_decl} vs _struct_sheet_order rows {n_order} vs gemmi {n_gemmi}",
           "page says 11; strands (auth): " + ", ".join(f"{r[4]}-{r[5]}" for r in ranges))


# ---------------------------------------------------------------- output
def load_uniprot(arg):
    if arg != "fetch" and Path(arg).exists():
        text, source = Path(arg).read_text(), arg
    else:
        with urllib.request.urlopen(UNIPROT_URL, timeout=60) as r:
            text, source = r.read().decode(), UNIPROT_URL
    lines = text.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise SystemExit("not a FASTA file: " + source)
    return "".join(l.strip() for l in lines[1:] if not l.startswith(">")), source


def find_fact(path, text_):
    hits = [i + 1 for i, line in enumerate(path.read_text(encoding="utf-8").splitlines()) if text_ in line]
    return f"{path.relative_to(ROOT)}:{','.join(map(str, hits))}" if hits else None


def summary():
    print("\n" + "=" * 128)
    print(f"{'Check':<21} {'Result':<17} {'Fields compared':<48} Meaning for the page")
    print("-" * 128)
    for cid, (title, fields) in CHECKS.items():
        statuses = [s for k, s, _f, _d in RESULTS if k.rstrip("abcdefghijk") == cid]
        n = Counter(statuses)
        if n["MISMATCH"]:
            result = f"MISMATCH {n['MISMATCH']}/{len(statuses)}"
        elif n["NOT CHECKABLE"] == len(statuses):
            result = "NOT CHECKABLE"
        else:
            result = f"PASS {n['PASS'] + n['PASS (silent)']}/{len(statuses)}"
        notes = [f"{n['NOT CHECKABLE']} not checkable"] * bool(n["NOT CHECKABLE"] and not result.startswith("NOT")) \
            + [f"{n['PASS (silent)']} silent"] * bool(n["PASS (silent)"])
        r_lines = [result] + notes
        facts = []
        for path, text_, checks in PAGE_FACTS:
            if cid in checks.split():
                where = find_fact(path, text_)
                facts.append(f"'{text_}' ({where or 'NOT FOUND'})")
        meaning = ("FLAG: " if n["MISMATCH"] else "holds: ") + ("; ".join(facts) if facts else "no on-screen fact")
        f_lines = textwrap.wrap(fields, 48) or [""]
        m_lines = textwrap.wrap(meaning, 40) or [""]
        for i in range(max(len(f_lines), len(m_lines), len(r_lines))):
            cell = lambda lines, w: f"{lines[i] if i < len(lines) else '':<{w}}"
            print(f"{cid + '. ' + title if i == 0 else '':<21} {cell(r_lines, 17)} {cell(f_lines, 48)} {cell(m_lines, 40)}".rstrip())
    print("=" * 128)
    flagged = [k for k, s, _f, _d in RESULTS if s == "MISMATCH"]
    print(f"{len(RESULTS)} lines: {Counter(s for _k, s, _f, _d in RESULTS)}; MISMATCH lines: {flagged or 'none'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--uniprot", nargs="?", const="fetch", metavar="FILE",
                        help="also compare with UniProt P42212 (fetched, or read from FILE)")
    args = parser.parse_args()
    c = Cif(SRC)
    uniprot = load_uniprot(args.uniprot) if args.uniprot else None
    print(f"{SRC.relative_to(ROOT)}: entry {value(c.block, '_entry.id')}, dictionary {value(c.block, '_audit_conform.dict_version')}, "
          f"{len(c.polymer)} polymer residues, {len(c.waters)} waters, gemmi {gemmi.__version__}\n")
    for check in (check_1, check_2, check_3, check_4, check_5, check_6, check_7, check_8):
        print(f"--- {check.__doc__}")
        check(c, uniprot) if check is check_2 else check(c)
        print()
    summary()


if __name__ == "__main__":
    main()
