# AF interface descriptors

These descriptors are reported by `alphafold-pdockq2-patches.py`. They do not estimate binding probability, affinity, or interaction energy. The former standalone comparison command has been retired.

## Run

```sh
sieve-py scripts/alphafold-pdockq2-patches.py predictions.zip
```

CIF input with `--full-data` and directories of ZIPs are also supported. Summaries go to stdout; `--output patches.tsv` optionally writes the existing patch TSV. See [pDockQ2.md](pDockQ2.md#structural-descriptors-in-the-patch-stdout-summary) for arguments and output definitions.

DSSP is found on PATH or at `tmp/dssp/bin/mkdssp`; use `--dssp PATH` to override or `--no-dssp` to omit secondary structure. Missing DSSP is explicitly reported as unassigned; a located executable that fails stops the run. Context defaults are `--sequence-flank 5`, `--spatial-radius 8`, and `--sasa-points 240`.

Other numerical thresholds are defined by `DescriptionOptions`: salt 4 Å, nonpolar 4.5 Å, van der Waals overlap 0.6 Å. These are exploratory geometric conventions, not learned cutoffs or validated binding classifiers.

## Data and confidence

Each ZIP's matching CIF and full confidence JSON are read with the existing AF input parser. CIF coordinates supply residue identity, chain ID, numbering, Cα and side-chain geometry. CIF Cα B-factor values supply pLDDT. The confidence JSON supplies directional PAE, preserving original token indexing. The request JSON and summary ipTM are not used in the descriptors. Sequence grouping uses the modeled CIF sequence, not unmodeled full-length sequence. Only the first coordinate model in each CIF is used.

Selected chains must exist and differ. Canonical amino acids are required; insertion codes and duplicate numbers are rejected. Coordinates and PAE must be finite, PAE nonnegative, and pLDDT within 0–100. Missing side chains are not reconstructed; completeness should be checked before interpreting atom-level counts. Other chains are excluded from the pairwise surface-burial calculation, while DSSP runs on the supplied complex.

Each patch retains the existing directional pDockQ2, PAE-only confidence, contact distances, and graph metadata. Chemistry and context are reported separately from prediction confidence; no alignment is performed.

## Interface descriptors

1. **Contacts and patches.** A residue pair contacts when its Cα distance is ≤ the chosen cutoff. The existing recursive graph algorithm supplies the full interface and all intermediate/leaf patch nodes, with the requested bridge and minimum-size settings. Descriptions use the patch's contacting residue sets. Full nodes can exist even when no leaf qualifies.
2. **Secondary structure.** DSSP assigns residue codes from coordinates. Contiguous, consecutively numbered runs of identical codes are enumerated from N to C independently per chain: H → alpha1, alpha2; E → beta1, beta2; G → 3_10_1; I → pi; B → beta_bridge; P → polyproline. Other assigned codes are grouped into separately numbered loop runs. A beta run is a strand, not an entire sheet; sheet topology is not inferred. Labels are local to a model and are not cross-species correspondence. A patch can intersect multiple runs. [DSSP implementation](https://github.com/PDB-REDO/dssp), [Biopython DSSP wrapper](https://biopython.org/docs/1.81/api/Bio.PDB.DSSP.html).
3. **Hydropathy.** Arithmetic mean of the existing Kyte–Doolittle scale over unique residues. Also report the fraction in the explicitly chosen set AVILMFWYPC. This set is a pipeline convention; it is not a binding score. Positive mean hydropathy means more hydrophobic on this scale. [Kyte and Doolittle, 1982](https://pubmed.ncbi.nlm.nih.gov/7108955/).
4. **Charge proxy.** Count K/R as +1 and D/E as −1; other residues, including H, as zero. Report H separately. This is nominal side-chain charge, not protonation modeling or electrostatic potential. Termini, pKa shifts, solvent, salt, and screening are omitted.
5. **Context.** Upstream/downstream windows surround every patch residue, including discontinuous patches. Remove patch residues and deduplicate within each set. Their union is the sequence shell; upstream and downstream sets may overlap. The spatial shell contains other residues of that chain with any heavy atom within the radius of a patch heavy atom. Shells exclude the patch itself. Each reports composition, hydropathy, and charge.
6. **Buried area.** Biopython Shrake–Rupley SASA with a 1.4 Å probe: for each atom, ΔSASA = SASA(chain alone) − SASA(selected two-chain complex), at unchanged coordinates. Sum per residue and per chain, in Å². No chain is relaxed or cropped to a patch. Patch area sums contributions of its selected residues and may omit other buried residues. Values are reported separately for each side; they are not halved or called a single interface area. Carbon/sulfur contributions are also recorded. Sample-point approximation affects precision; repeat at higher point counts when small differences matter.
7. **Salt candidates.** Interchain Lys NZ/Arg NE,NH1,NH2 to Asp OD1,OD2/Glu OE1,OE2 at ≤4 Å. Each atom pair is retained; residue-pair summaries deduplicate by interaction type. Histidine and protonation-dependent assignments are omitted.
8. **Nonpolar candidates.** Interchain side-chain C/S atoms (excluding N,CA,C,O,OXT names) in residues AVILMFWYPC at ≤4.5 Å. This deliberately restricted definition is not a complete description of hydrophobic packing, aromatic stacking, or cation–π geometry.
9. **Steric-overlap candidates.** For C,N,O,S use radii 1.70,1.55,1.52,1.80 Å; flag radius sum minus distance >0.6 Å. SG–SG distances 1.8–2.2 Å are excluded as possible disulfides. This is a geometric screen, not MolProbity clashscore. Hydrogen bonds, protonation, other covalent links, and atom-specific radii are not resolved. Do not count overlaps as favorable interaction evidence.

Atom-level features are calculated for the selected chain pair. Patch features retain atom pairs whose two residue endpoints belong to the patch; they need not independently satisfy the Cα cutoff. Thus atom-contact counts and Cα-contact counts describe different geometries.

## Interpretation

The weakest links in contact chemistry are atom completeness, side-chain geometry, and protonation. These are geometric candidates, not confirmed bonds. Check their consistency across predicted models and, when available, against experimental structures. Nominal charge is not electrostatic potential, and hydropathy is not binding energy.

## Local DSSP installation used for the example

DSSP 4.6.1 was built from the official PDB-REDO/dssp source with CMake, FastFloat 8.0.2, and fetched libcifpp/Eigen dependencies. The executable is `tmp/dssp/bin/mkdssp`; this ignored local installation is not bundled with the project. The runtime dictionary directory contains the installed PDBx and DSSP dictionaries, the official ModelCIF dictionary, and libcifpp's `ccd-subset.cif` copied as `components.cif` (canonical amino acids sufficient for these inputs). Install a complete supported DSSP distribution for broader chemical-component coverage. Every residue in all ten supplied models received a DSSP assignment in the example run.

For single-input analysis, use the expanded stdout summary in `alphafold-pdockq2-patches.py`; see [pDockQ2.md](pDockQ2.md#structural-descriptors-in-the-patch-stdout-summary). It does not perform alignment or require a reference ensemble.
