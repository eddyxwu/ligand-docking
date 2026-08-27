# Ligand Docking — Lysozyme + Penta-NAG (Coarse-Grained, PDB 1SFB)

A coarse-grained molecular-recognition simulation: **hen egg-white lysozyme re-binding
penta-N-acetylchitopentaose, (NAG)₅** — a five-ring fragment of its natural bacterial
cell-wall substrate. Both molecules are built from the same crystal structure
(**PDB 1SFB**), and the docking is **emergent, not animated**: the ligand starts
displaced from its crystallographic pose and re-enters the enzyme's binding groove
under nothing but a Lennard-Jones potential and shape complementarity.

Built as an extension to Lagra, a Rust physics engine.

## Scope — read this first

This repository contains **everything I wrote for this scene**, and only that.
The engine itself is private: the patches in `engine-patches/` document the physics
I added and **do not compile standalone**. The scene generator, however, **runs
standalone** — it fetches the crystal structure and emits the exact scene the
engine consumes:

```bash
python3 generator/gen_ligand_docking_scene.py
# downloads PDB 1SFB (if needed) -> writes scene/ligand_docking.json
```

## The model

| Component | Model |
|---|---|
| **Receptor** | Lysozyme's 129-residue Cα trace (one bead per residue, real 1SFB coordinates), each bead harmonically **position-restrained** to its crystal position with small thermal velocities — restrained MD. The protein visibly fluctuates (~0.4 σ RMSD, the coarse analog of crystallographic B-factors) but its fold cannot degrade. |
| **Ligand** | Penta-NAG at three beads per sugar (pyranose ring, N-acetyl arm, hydroxymethyl group = 15 beads), its crystal geometry held by **harmonic bonds**: stiff intra-sugar triangles, moderate glycosidic links, and soft graded backbone braces — the chain flexes at its linkages like a real oligosaccharide. Bonded pairs are excluded from nonbonded LJ (the standard 1-2 exclusion). |
| **Interaction** | A single truncated Lennard-Jones potential (ε, σ, 2.5 σ cutoff from config). Intermolecular pairs use ε × 0.38 (cohesion > adhesion — the coarse-grained proxy for covalent bonding vs. van der Waals adhesion). No electrostatics, no solvent, no hydrogen bonds — deliberately minimal. |
| **Experiment** | **Pose rebinding**: the ligand starts 6 σ outside its crystallographic pose and re-binds. The binding reaction coordinate — \|COM(ligand) − centroid(site residues)\| — settles into the bound basin around the crystal-pose value (2.07 σ). |

## Runtime validation (measured, full 60k-step run)

| Validator | Result | Meaning |
|---|---|---|
| Energy conservation | ✅ 1.43 % drift | E = KE + V_LJ + V_restraint + V_bond — every Hamiltonian term bookkept |
| Receptor RMSD vs crystal | ✅ 0.39 σ plateau | fold holds at thermal fluctuation level |
| Binding reaction coordinate | ✅ 1.66 σ (bound) / ❌ 7.19 σ (unbound check) | certifies the docking event; fails during approach *by construction* and flips on binding |
| Charge conservation | ✅ exact | |
| Subluminal velocities | ✅ max 0.31 c | integrator health |
| Momentum conservation | ❌ **expected failure, documented** | position restraints are external forces — broken translational symmetry (Noether), measured drift matches the ligand momentum absorbed by the restraints |

The deliberate momentum failure is the point, not a bug: restrained dynamics
genuinely does not conserve momentum, and the validator measures exactly that.

## What I built (see `engine-patches/CONTRIBUTIONS.md` for the guided tour)

- **Lennard-Jones interaction** with config-driven ε/σ/cutoff and intermolecular ε scaling
- **Harmonic bonds** with standard 1-2 bonded exclusions (internal forces — momentum-safe)
- **Position restraints** (restrained MD) with correct Hamiltonian bookkeeping
- **Two validators**: RMSD-vs-reference and a binding reaction-coordinate outcome check
- **PDB-driven scene generator**: parsing, coarse-graining, bond-topology construction, staging, validator calibration
- **Shaded-sphere rendering** for coarse-grained particles (visual only)

## Repository map

```
generator/        the standalone PDB -> scene generator (start here)
scene/            its output: the scene the engine consumes (146 entities, 25 bonds)
configs/          clean excerpts: LJ preset, bead definitions, validator entries
data/             PDB 1SFB (public crystallographic data)
engine-patches/   my engine additions as numbered patches + reading guide
```

## Data & citations

- Crystal structure: **PDB 1SFB** — Von Dreele, R.B., *Binding of N-acetylglucosamine
  oligosaccharides to hen egg-white lysozyme: a powder diffraction study*,
  Acta Crystallogr. D (2005). Data from [RCSB PDB](https://www.rcsb.org/structure/1SFB).
- Lysozyme background: the first enzyme structure ever solved (Phillips, 1965);
  (NAG)₅ spans binding sub-sites A–E of its substrate-binding cleft.

## License

MIT — see [LICENSE](LICENSE).
