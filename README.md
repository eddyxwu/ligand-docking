# ligand-docking

A coarse-grained simulation of **lysozyme catching its sugar substrate**. Both
molecules come from a real crystal structure (PDB 1SFB) — the protein as one bead
per residue, the five-ring sugar (penta-NAG) as a flexible bonded chain. Pull the
sugar out of the binding groove, let physics run, and it finds its way back to
where the X-ray says it belongs. Nothing is scripted — the docking emerges from a
Lennard-Jones potential and shape complementarity alone.

Built as an extension to Lagra, a Rust physics engine.

## Try it

```bash
python3 generator/gen_ligand_docking_scene.py
```

Downloads the crystal structure and builds the complete scene (146 particles,
25 bonds) from scratch. No dependencies beyond Python 3.

## What's here

```
generator/       PDB -> scene builder: parsing, coarse-graining, bond topology, staging
scene/           the generated scene the engine runs
configs/         the physics knobs: LJ potential, bead types, validators
data/            PDB 1SFB (public crystallographic data)
engine-patches/  the ~520 lines I added to the engine, with a reading guide
```

The engine itself is private — the patches show my additions (Lennard-Jones
interaction, harmonic bonds with 1-2 exclusions, position restraints, and two
MD-style validators) but won't compile on their own. The generator runs anywhere.

## Does it work?

The simulation checks itself while running:

- energy conserved to **1.4%** over the full run (every potential term bookkept)
- protein holds **0.39σ RMSD** vs the crystal reference (it jiggles, it doesn't unfold)
- the binding coordinate drops from 7σ to the crystallographic bound basin (**~2σ**)
- one *deliberate* failure: momentum — position restraints are external forces,
  so it shouldn't be conserved, and the validator confirms exactly that

## Credit

Crystal structure: PDB [1SFB](https://www.rcsb.org/structure/1SFB) (Von Dreele,
powder diffraction). Lysozyme is the first enzyme ever solved (1965); penta-NAG
is a fragment of the bacterial cell-wall sugar it naturally cuts.

MIT license.
