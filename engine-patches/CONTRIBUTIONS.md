# Edits overview

Each patch is a `git diff` of one area of the Lagra Rust engine. They are numbered
in physics-first reading order. Every force term follows the engine's existing
convention: the functional form lives in Rust with a derivation comment; every
parameter comes from JSON config.

## `1_lennard_jones.patch` (`src/lagrangians.rs`, ~100 lines)

The Lennard-Jones interaction as a first-class Lagrangian term:
**V = 4ε[(σ/r)¹² − (σ/r)⁶]**, force **F = 24ε[2(σ/r)¹² − (σ/r)⁶]/r**, truncated at
`cutoff_factor × σ`. All parameters (ε, σ, cutoff, inter-group factor) parse from JSON.

Two design points worth noting:

- **Intermolecular ε scaling** (`inter_group_factor`): pairs in *different* molecules
  attract at 0.38 × the intramolecular strength. Discovered the hard way — with
  symmetric ε the ligand *wets* the receptor surface like a droplet on glass
  (per-bead adhesion beats cohesion). Real ligands stay intact because covalent
  bonds beat van der Waals adhesion ~100×; the scaling is the coarse-grained proxy,
  conceptually the same as intra/inter nonbonded scaling in molecular force fields.
- The `bonded` flag on the force context implements the exclusion described in patch 2.

## `2_bonds_restraints_energy.patch` (`src/bin/renderer/world.rs`, ~140 lines)

Three mechanics:

- **Harmonic bonds** `V = ½k(r − r₀)²`, applied equal-and-opposite (Newton's third
  law → internal forces → total momentum unaffected). Used to hold the ligand's
  crystal geometry: stiff intra-sugar triangles, glycosidic links, soft graded braces.
- **1-2 bonded exclusions**: bonded pairs interact through their bond *only*, never
  through nonbonded LJ. This was forced by physics, not taste — coarse-grained bond
  lengths (0.5–0.9 σ) sit deep inside the LJ repulsive wall, and the first unexcluded
  run blew the ligand apart at 40 σ/frame. Every real force field applies this same
  exclusion for the same reason.
- **Position restraints** (restrained MD): each receptor bead is tethered
  harmonically to its crystal coordinate, `F = −k(p − p₀)`. These are *external*
  forces: they anchor the protein to the lab frame and deliberately break momentum
  conservation (Noether) — measured and documented by the validator rather than hidden.
- **Hamiltonian bookkeeping**: both the restraint PE and bond PE are added to
  `total_potential_energy`, so the energy validator measures true drift (1.43 %)
  instead of reporting the restraint energy as spurious drift.

## `3_validators.patch` (`src/bin/renderer/runtime_validator.rs`, ~53 lines)

Two new validator types, chosen to mirror how the MD field actually validates:

- **`rmsd_from_reference`** — RMSD of restrained particles vs. their crystal
  coordinates: the canonical structural-stability metric (measured plateau: 0.39 σ).
- **`binding_site_distance`** — the binding reaction coordinate:
  |COM(ligand) − centroid(site residues)|. An *outcome* validator: it fails during
  the approach by construction and flips to pass when the ligand enters the bound
  basin (tolerance calibrated against the 1SFB crystal pose: 2.07 σ + thermal margin).

## `4_scene_plumbing.patch` (`config.rs` + `scene_loader.rs`, ~45 lines)

Scene-format extensions: per-entity `group` (molecule membership, drives the
inter-group ε scaling), per-entity `tether_k` (restraints), and a scene-level
`bonds` list `(i, j, r₀, k)`. The generator emits all three; the engine stays
free of any scene-specific knowledge.

## `5_renderer.patch` (render/, ~150 lines — visual only, zero physics)

Shaded-sphere rendering for coarse-grained particles (Lambert diffuse + Blinn
specular + antialiased edge + soft halo) replacing flat circles — the molecular
space-fill look. Plus an optional scene-label overlay (title + 3D-anchored
annotations), currently unused by the scene.

## `6_capture_plumbing.patch` (`main.rs`, ~26 lines)

Headless capture mode now honors `--render-profile` and `-f` (force-visualization
mode), so validated frame sequences can be rendered without the interactive window.
