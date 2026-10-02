"""
Turn the PDB entry 1EMA (green fluorescent protein) into the compact data
the page needs.

Input : data/1ema.cif   (mmCIF file from the Protein Data Bank, CC0)
Output: data/gfp.json   (atoms, backbone trace, secondary structure)

What it does
1. Reads the structure with gemmi and drops water molecules.
2. Collects the helices and the 11 beta strands annotated in the file.
3. Finds the barrel axis as the average direction of the 11 strands.
4. Rotates every coordinate so the barrel axis points up (Y) and the
   chromophore (residue CRO, amino acids 65-67) sits at the origin.
5. Writes atoms, the C-alpha trace and the chromophore atoms as JSON.

Run from the repository root:  python pipeline/prep.py
"""
import json
from pathlib import Path

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "1ema.cif"
OUT = ROOT / "data" / "gfp.json"

# Chromophore atoms that form the conjugated rings (the part that absorbs and emits light)
CONJUGATED = ["C1", "N2", "N3", "C2", "O2", "CA2", "CB2", "CG2", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]
ELEMENT_INDEX = {"C": 0, "N": 1, "O": 2, "S": 3, "Se": 4}


def main():
    st = gemmi.read_structure(str(SRC))
    st.remove_waters()
    st.remove_hydrogens()
    chain = st[0]["A"]

    # Secondary structure as annotated in the file
    ss = {}
    for h in st.helices:
        for n in range(h.start.res_id.seqid.num, h.end.res_id.seqid.num + 1):
            ss[n] = "H"
    strands = []
    for sheet in st.sheets:
        for x in sheet.strands:
            a, b = x.start.res_id.seqid.num, x.end.res_id.seqid.num
            strands.append((a, b))
            for n in range(a, b + 1):
                ss[n] = "E"

    # Barrel axis: average direction of the 11 strands (antiparallel ones flipped)
    by_num = {r.seqid.num: r for r in chain}
    directions = []
    for a, b in strands:
        pa, pb = by_num[a]["CA"][0].pos, by_num[b]["CA"][0].pos
        v = np.array([pb.x - pa.x, pb.y - pa.y, pb.z - pa.z])
        v /= np.linalg.norm(v)
        if directions and np.dot(v, directions[0]) < 0:
            v = -v
        directions.append(v)
    axis = np.mean(directions, axis=0)
    axis /= np.linalg.norm(axis)

    # The chromophore centre becomes the origin
    cro = next(r for r in chain if r.name == "CRO")
    centre = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in cro]).mean(axis=0)

    # Rotation that maps the barrel axis onto +Y
    y = axis
    x = np.cross([0, 0, 1], y)
    x /= np.linalg.norm(x)
    z = np.cross(x, y)
    rot = np.vstack([x, y, z])

    def transform(p):
        return rot @ (np.array(p) - centre)

    atoms, ca, cro_index = [], [], []
    for i, res in enumerate(chain):
        for atom in res:
            p = transform([atom.pos.x, atom.pos.y, atom.pos.z])
            if res.name == "CRO":
                cro_index.append(len(atoms))
            atoms.append([round(float(p[0]), 1), round(float(p[1]), 1), round(float(p[2]), 1),
                          ELEMENT_INDEX.get(atom.element.name, 0), i])
        ca_name = "CA2" if res.name == "CRO" else "CA"
        a = res[ca_name][0].pos
        p = transform([a.x, a.y, a.z])
        code = "X" if res.name == "CRO" else ss.get(res.seqid.num, "C")
        ca.append([round(float(p[0]), 2), round(float(p[1]), 2), round(float(p[2]), 2), code, res.seqid.num])

    data = {
        "rn": [r.name for r in chain],
        "croConj": [i for i, a in enumerate(cro) if a.name in CONJUGATED],
        "a": [v for atom in atoms for v in atom],
        "ca": [v for c in ca for v in c[:3]],
        "ss": "".join(c[3] for c in ca),
        "num": [c[4] for c in ca],
        "cro": [cro_index[0], cro_index[-1] + 1],
        "strands": strands,
    }
    OUT.write_text(json.dumps(data, separators=(",", ":")))
    print(f"{len(atoms)} atoms, {len(ca)} residues, {len(strands)} strands -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
