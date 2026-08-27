#!/usr/bin/env python3
"""Generate the ligand_docking scene - GROOVE-SPANNING SUBSTRATE edition.

STANDALONE COPY: identical to the in-engine script except for three path
adaptations - the PDB is read from ../data/1SFB.pdb (auto-downloaded from
RCSB if missing) and the scene JSON is written to ../scene/. No other changes.

Both molecules come from the same crystal structure, PDB 1SFB (wild-type hen
egg-white lysozyme complexed with penta-N-acetylchitopentaose, (NAG)5 - a
room-temperature powder-diffraction study; the pentasaccharide occupies
binding sub-sites A-B-C-D-E exactly, spanning the groove with no overhang):

  Receptor: chain-A C-alpha trace, one bead per residue, harmonically
            restrained to the crystal coordinates (restrained MD) with small
            deterministic thermal velocities.
  Ligand:   penta-NAG coarse-grained at three beads per sugar (pyranose ring,
            N-acetyl arm, hydroxymethyl group) = 15 beads. Intra-sugar bead
            triangles are stiff; glycosidic ring-ring links and soft 1-3
            braces let the chain FLEX at its linkages - real oligosaccharide
            flexibility - while resisting jackknife folding.

The ligand starts displaced outward from its crystallographic pose and
RE-BINDS into the groove - pose reproduction. The binding_reaction_coordinate
validator is calibrated against the crystal-pose value printed by this script.
"""
import json, math, os

SIGMA = 1.0
D = 2 ** (1 / 6) * SIGMA
CA_SPACING_A = 3.8
SCALE = D / CA_SPACING_A

PDB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "1SFB.pdb")
if not os.path.exists(PDB):
    import urllib.request
    os.makedirs(os.path.dirname(PDB), exist_ok=True)
    print("downloading 1SFB from RCSB...")
    urllib.request.urlretrieve("https://files.rcsb.org/download/1SFB.pdb", PDB)

# Camera sits on -Z (yaw 0, pitch 0.3; +x renders screen-LEFT). The protein is
# rotated so the cleft points along U: screen-right, up, toward the camera.
_u = (-0.62, 0.28, -0.76)
_n = math.sqrt(sum(c * c for c in _u))
U = tuple(c / _n for c in _u)

SITE_RESIDUES = [35, 52, 62, 63, 101, 108]   # cleft sub-sites A-F lining, Glu35/Asp52
SITE_PAINT_RADIUS = 2.6
TETHER_K = 100.0
JIGGLE_V = 0.030        # receptor thermal speed
LIG_JIGGLE_V = 0.025    # ligand internal thermal speed
K_INTRA = 250.0         # intra-sugar bead triangle (rigid unit)
K_LINK = 150.0          # glycosidic ring-ring links
# <ring-backbone elastic network: progressively softer with separation.
#  Bonded ring pairs are LJ-excluded, which removes the self-attraction that
#  would otherwise collapse the chain (poor-solvent curl), while soft long
#  braces still allow visible bending/flexing. Arms stay unbonded across
#  sugars: free to wiggle.>
K_BRACES = {2: 80.0, 3: 50.0, 4: 35.0, 5: 25.0}
LIGAND_SPEED = 0.045
DISPLACE = 6.0

RING_ATOMS = {"C1", "C2", "C3", "C4", "C5", "O5", "O1", "O3", "O4"}
ACETYL_ATOMS = {"N2", "C7", "C8", "O7"}
HYDROXYMETHYL_ATOMS = {"C6", "O6"}

def cross(p, q): return (p[1]*q[2]-p[2]*q[1], p[2]*q[0]-p[0]*q[2], p[0]*q[1]-p[1]*q[0])
def dot(p, q): return sum(x*y for x, y in zip(p, q))
def dist(p, q): return math.dist(p, q)
def centroid(pts):
    n = len(pts)
    return tuple(sum(p[i] for p in pts) / n for i in range(3))

# --- parse chain-A C-alpha trace + chain-B NAG heavy atoms ---
ca, site_idx_by_resseq = [], {}
nag = {}   # resseq -> {"ring": [...], "acetyl": [...], "hm": [...]}
with open(PDB) as f:
    for line in f:
        if line.startswith("ATOM") and line[12:16] == " CA " and line[21] == "A":
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            site_idx_by_resseq[int(line[22:26])] = len(ca)
            ca.append((x, y, z))
        elif line.startswith("HETATM") and line[17:20] == "NAG" and line[21] == "B":
            name = line[12:16].strip()
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            r = nag.setdefault(int(line[22:26]), {"ring": [], "acetyl": [], "hm": []})
            if name in ACETYL_ATOMS:
                r["acetyl"].append((x, y, z))
            elif name in HYDROXYMETHYL_ATOMS:
                r["hm"].append((x, y, z))
            else:
                # <ring atoms + any leftover hydroxyls fold into the ring bead>
                r["ring"].append((x, y, z))
assert len(ca) == 129, len(ca)
assert len(nag) == 5, f"expected penta-NAG (5 residues), got {len(nag)}"

# ligand beads in crystal frame: list of (position, kind, sugar_index)
# per sugar: [ring, acetyl, hm] in that order (missing groups skipped)
lig_raw = []
ring_bead_of_sugar = {}
for si, resseq in enumerate(sorted(nag)):
    for kind in ("ring", "acetyl", "hm"):
        pts = nag[resseq][kind]
        if pts:
            if kind == "ring":
                ring_bead_of_sugar[si] = len(lig_raw)
            lig_raw.append((centroid(pts), kind, si))

# --- shared transform: center on protein centroid, scale, rotate cleft -> U ---
cx = sum(p[0] for p in ca) / len(ca)
cy = sum(p[1] for p in ca) / len(ca)
cz = sum(p[2] for p in ca) / len(ca)
def xform0(p): return ((p[0]-cx)*SCALE, (p[1]-cy)*SCALE, (p[2]-cz)*SCALE)
ca_s = [xform0(p) for p in ca]
lig_s = [(xform0(p), kind, si) for p, kind, si in lig_raw]

site_pts = [ca_s[site_idx_by_resseq[r]] for r in SITE_RESIDUES]
sc = centroid(site_pts)
sn = math.sqrt(dot(sc, sc))
a = tuple(c / sn for c in sc)
v = cross(a, U); c_ = dot(a, U); s = math.sqrt(dot(v, v))
def rotate(p):
    if s < 1e-12:
        return p if c_ > 0 else (-p[0], -p[1], -p[2])
    k = tuple(x / s for x in v)
    kxp = cross(k, p); kdp = dot(k, p)
    return tuple(p[i]*c_ + kxp[i]*s + k[i]*kdp*(1-c_) for i in range(3))

ca_f = [rotate(p) for p in ca_s]
site_center = rotate(sc)
lig_crystal = [(rotate(p), kind, si) for p, kind, si in lig_s]

# --- crystal-pose reaction coordinate (validator calibration) ---
lig_com = centroid([p for p, _, _ in lig_crystal])
site_bead_idx = [i for i, p in enumerate(ca_f) if dist(p, site_center) < SITE_PAINT_RADIUS]
site_bead_centroid = centroid([ca_f[i] for i in site_bead_idx])
crystal_rc = dist(lig_com, site_bead_centroid)

# --- displaced start pose along the ligand's own outward radial direction ---
rn = math.sqrt(dot(lig_com, lig_com))
u_out = tuple(c / rn for c in lig_com)
disp = DISPLACE
while True:
    start = [(tuple(p[i] + disp * u_out[i] for i in range(3)), kind, si) for p, kind, si in lig_crystal]
    gap = min(dist(p, q) for q in ca_f for p, _, _ in start)
    if gap > 2.6:
        break
    disp += 0.5
drift = tuple(-LIGAND_SPEED * c for c in u_out)

# --- flexible bond topology from crystal geometry (indices offset 129) ---
bonds = []
def add_bond(a_, b_, k):
    r0 = dist(lig_crystal[a_][0], lig_crystal[b_][0])
    bonds.append({"i": 129 + a_, "j": 129 + b_, "r0": round(r0, 4), "k": k})

by_sugar = {}
for idx, (_, kind, si) in enumerate(lig_crystal):
    by_sugar.setdefault(si, []).append(idx)
# intra-sugar triangles (stiff)
for si, idxs in by_sugar.items():
    for x in range(len(idxs)):
        for y in range(x + 1, len(idxs)):
            add_bond(idxs[x], idxs[y], K_INTRA)
# glycosidic ring-ring links (moderate) + ring-backbone braces (soft, graded)
sugars = sorted(by_sugar)
for si in sugars[:-1]:
    add_bond(ring_bead_of_sugar[si], ring_bead_of_sugar[si + 1], K_LINK)
for sep, kb in K_BRACES.items():
    for si in sugars[:-sep]:
        add_bond(ring_bead_of_sugar[si], ring_bead_of_sugar[si + sep], kb)

ext = max(dist(a_[0], b_[0]) for a_ in lig_crystal for b_ in lig_crystal)
print(f"protein beads: {len(ca_f)}  ligand beads: {len(lig_crystal)}  bonds: {len(bonds)}  ligand extent: {ext:.2f} sigma")
print(f"crystal-pose reaction coordinate: {crystal_rc:.2f} sigma  (validator bound-basin reference)")
print(f"start displacement: {disp:.1f} sigma  gap to protein: {gap:.2f} (want > 2.5)")

# --- deterministic thermal velocities (net zero) ---
def h3(i):
    def h(n):
        n = (n ^ 0x9E3779B9) * 2654435761 % (2**32)
        n = (n ^ (n >> 16)) * 2246822519 % (2**32)
        return ((n ^ (n >> 13)) % 10007) / 10007.0 * 2.0 - 1.0
    return (h(3*i+1), h(3*i+2), h(3*i+3))

def net_zero(vs):
    out = vs
    for k in range(3):
        m = sum(vv[k] for vv in out) / len(out)
        out = [tuple(vv[j] - m if j == k else vv[j] for j in range(3)) for vv in out]
    return out

rec_vels = [tuple(JIGGLE_V * c for c in vv) for vv in net_zero([h3(i) for i in range(len(ca_f))])]
# ligand: thermal component (net zero) + shared drift => COM moves at drift exactly
lig_therm = [tuple(LIG_JIGGLE_V * c for c in vv) for vv in net_zero([h3(1000 + i) for i in range(len(lig_crystal))])]
lig_vels = [tuple(drift[k] + t[k] for k in range(3)) for t in lig_therm]

# --- emit entities ---
def ent(name, p, vel=None, group=None, tether=None):
    e = {"particle": name, "position": [round(c, 4) for c in p]}
    if vel:
        e["velocity"] = [round(c, 5) for c in vel]
    if group is not None:
        e["group"] = group
    if tether is not None:
        e["tether_k"] = tether
    return e

R_VARIANTS = ["receptor_bead", "receptor_bead_md", "receptor_bead_dk"]
S_VARIANTS = ["receptor_site", "receptor_site_dk"]

entities = []
for i, p in enumerate(ca_f):
    is_site = i in site_bead_idx
    name = (S_VARIANTS[(i * 2654435761) % 2] if is_site
            else R_VARIANTS[(i * 2654435761) % 3])
    entities.append(ent(name, p, rec_vels[i], group=1, tether=TETHER_K))
print(f"site-highlighted beads: {len(site_bead_idx)}")

for i, (p, kind, si) in enumerate(start):
    if kind == "ring":
        name = "ligand_bead" if si % 2 == 0 else "ligand_bead_dk"
    elif kind == "acetyl":
        name = "ligand_bead_o"
    else:
        name = "ligand_bead"
    entities.append(ent(name, p, lig_vels[i], group=2))

start_com = centroid([p for p, _, _ in start])
anchor_r = round(max(12.0, max(math.sqrt(dot(p, p)) for p, _, _ in start) + 3.0), 1)
entities += [ent("frame_anchor", (anchor_r, 0.0, 0.0)), ent("frame_anchor", (-anchor_r, 0.0, 0.0))]

scene = {
    "name": "Lysozyme + penta-NAG (1SFB pose rebinding)",
    "description": "Hen egg-white lysozyme re-binding penta-N-acetylchitopentaose, (NAG)5, both coarse-grained from the same crystal structure (PDB 1SFB). The receptor is the chain-A C-alpha trace under restrained thermal dynamics; the ligand is the pentasaccharide at three beads per sugar with flexible glycosidic linkages, starting displaced from its crystallographic pose and re-binding across binding sub-sites A-B-C-D-E - the groove's five occupied sites, spanned exactly, with no overhanging sugar.",
    "preset": "standard_model/molecular_dynamics/ligand_docking",
    "environment": "none",
    "lagrangian": {"interactions": {"qdesic": {"enabled": False}}},
    "_physics_note": "Pose-rebinding of the groove-spanning substrate: penta-NAG starts displaced outward from its 1SFB crystallographic pose and re-binds under Lennard-Jones attraction + shape complementarity. Ligand topology: stiff intra-sugar bead triangles, moderate glycosidic ring-ring links, soft 1-3 braces - the chain flexes at its linkages like a real oligosaccharide while resisting jackknife folding. Bonded pairs are excluded from nonbonded LJ (standard 1-2 exclusion). Registration note: the chain reaches the crystallographic reaction-coordinate value (~2.3 sigma) and then slides ~1.5 sigma deeper along the groove to the contact-maximizing position - with a uniform LJ potential there are no site-specific hydrogen bonds to pin the exact crystallographic registration. Known coarse-graining limitation, not an integrator artifact. Receptor beads are position-restrained to crystal coordinates (external forces: momentum conservation intentionally broken, see expected_fail). Bond, restraint, and LJ potentials are all included in the Hamiltonian for the energy validator.",
    "entities": entities,
    "bonds": bonds,
    "solver": {"dt": 0.005, "steps_per_frame": 40},
    "expected_pass": [
        "charge_conservation",
        "subluminal_velocity",
        "energy_stability",
        "receptor_rmsd",
        "binding_reaction_coordinate"
    ],
    "_validator_notes": [
        "energy_stability: E = KE + V_LJ + V_restraint + V_bond conserved (all Hamiltonian terms bookkept).",
        "receptor_rmsd: restrained C-alpha beads hold ~0.4-0.5 sigma RMSD vs the 1SFB crystal reference.",
        f"binding_reaction_coordinate: |COM(ligand) - centroid(site)| re-enters the bound basin; crystal-pose value {crystal_rc:.2f} sigma. FAILS during approach by construction, flips to PASS at rebinding."
    ],
    "expected_fail": ["momentum_conservation"],
    "_momentum_note": "Position restraints anchor the receptor to the lab frame: external forces break translational symmetry (Noether), so momentum conservation is expected to fail once the ligand's momentum is absorbed by the restraints. Harmonic bonds are internal (Newton's third law) and momentum-safe. This is correct physics for restrained dynamics, not an integrator defect.",
    "terminate_when": {"steps_above": 60000}
}

out = os.path.join(os.path.dirname(os.path.abspath(PDB)), "..", "scene", "ligand_docking.json")
with open(out, "w") as f:
    json.dump(scene, f, indent=2)
print(f"wrote {os.path.normpath(out)} with {len(entities)} entities, {len(bonds)} bonds")
