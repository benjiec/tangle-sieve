# AF interface description and comparison, version 1

This exploratory pipeline describes two AlphaFold ensembles and compares their interfaces under an explicit chain correspondence. It does not estimate binding probability, affinity, or interaction energy. Prediction confidence and structural similarity are separate outputs. Existing pDockQ2 scoring is unchanged; its equations and patch graph algorithm are documented in [pDockQ2.md](pDockQ2.md).

## Run

From the `sieve` directory, using the sieve Python environment:

```sh
sieve-py scripts/alphafold-interface-compare.py \
  data/fold_hs_myd88_dd_1x_hs_irak4_dd_1x.zip \
  data/fold_nv_myd88_dd_irak_xp_048575654_dd.zip \
  --reference-chains A B --query-chains A B \
  --output-dir tmp/myd88-comparison
```

The positional inputs are the reference (known-protein AF) ZIP and query (novel-protein AF) ZIP. Use a new output directory. All models are analyzed; five models on each side produce 25 comparisons. The first selected chain anchors the structural superposition. `--reference-chains A B --query-chains X Y` explicitly means A corresponds to X and B to Y. It does not mean residue numbers match. For multiple copies, choose the biological pair explicitly; symmetry-related chain assignments are not searched.

Arguments:

| Argument | Default | Meaning |
|---|---|---|
| `--reference-chains`, `--query-chains` | required | Two ordered chain IDs in each ZIP |
| `--output-dir` | required | New report directory; existing paths refused |
| `--distance-cutoff` | 8 Å | Inclusive interchain Cα contact and patch graph threshold |
| `--sequence-flank` | 5 | Residues before/after each patch residue in coordinate sequence order |
| `--spatial-radius` | 8 Å | Same-chain neighborhood, minimum heavy-atom distance to patch |
| `--sasa-points` | 240 | Surface sample points per atom |
| `--dssp` | PATH or `tmp/dssp/bin/mkdssp` | DSSP executable |
| `--no-dssp` | false | Explicitly omit secondary structure |

DSSP errors stop the run. The pipeline never silently substitutes a secondary-structure prediction. Other numerical thresholds are defined by `DescriptionOptions`: salt 4 Å, nonpolar 4.5 Å, van der Waals overlap 0.6 Å. These are exploratory geometric conventions, not learned cutoffs or validated binding classifiers.

## Data and confidence

Each ZIP's matching CIF and full confidence JSON are read with the existing AF input parser. CIF coordinates supply residue identity, chain ID, numbering, Cα and side-chain geometry. CIF Cα B-factor values supply pLDDT. The confidence JSON supplies directional PAE, preserving original token indexing. The request JSON and summary ipTM are not used in this comparison. Sequence correspondence uses the modeled CIF sequence, not unmodeled full-length sequence. Only the first coordinate model in each CIF is used.

Selected chains must exist and differ. Sequences and residue numbering must remain consistent across ensemble models. Canonical amino acids are required; insertion codes and duplicate numbers are rejected. Coordinates and PAE must be finite, PAE nonnegative, and pLDDT within 0–100. Missing backbone atoms are listed per residue. Missing side chains are not reconstructed; completeness should be checked before interpreting atom-level counts. Other chains are excluded from the pairwise surface-burial calculation, while DSSP runs on the supplied complex.

Each patch retains the existing directional pDockQ2, PAE-only confidence, contact distances, and graph metadata. The comparison F1 and coordinate fits do **not** use pLDDT or PAE. Confidence is not used to select a favorable alignment or discard poorly matching residues.

## Interface descriptors

1. **Contacts and patches.** A residue pair contacts when its Cα distance is ≤ the chosen cutoff. The existing recursive graph algorithm supplies the full interface and all intermediate/leaf patch nodes, with its default bridge and minimum-size settings. Descriptions use the patch's contacting residue sets. Full nodes can exist even when no leaf qualifies.
2. **Secondary structure.** DSSP assigns residue codes from coordinates. Contiguous, consecutively numbered runs of identical codes are enumerated from N to C independently per chain: H → alpha1, alpha2; E → beta1, beta2; G → 3_10_1; I → pi; B → beta_bridge; P → polyproline. Other assigned codes are grouped into separately numbered loop runs. A beta run is a strand, not an entire sheet; sheet topology is not inferred. Labels are local to a model and are not cross-species correspondence. A patch can intersect multiple runs. [DSSP implementation](https://github.com/PDB-REDO/dssp), [Biopython DSSP wrapper](https://biopython.org/docs/1.81/api/Bio.PDB.DSSP.html).
3. **Hydropathy.** Arithmetic mean of the existing Kyte–Doolittle scale over unique residues. Also report the fraction in the explicitly chosen set AVILMFWYPC. This set is a pipeline convention; it is not a binding score. Positive mean hydropathy means more hydrophobic on this scale. [Kyte and Doolittle, 1982](https://pubmed.ncbi.nlm.nih.gov/7108955/).
4. **Charge proxy.** Count K/R as +1 and D/E as −1; other residues, including H, as zero. Report H separately. This is nominal side-chain charge, not protonation modeling or electrostatic potential. Termini, pKa shifts, solvent, salt, and screening are omitted.
5. **Context.** Upstream/downstream windows surround every patch residue, including discontinuous patches. Remove patch residues and deduplicate within each set. Their union is the sequence shell; upstream and downstream sets may overlap. The spatial shell contains other residues of that chain with any heavy atom within the radius of a patch heavy atom. Shells exclude the patch itself. Each reports composition, hydropathy, and charge.
6. **Buried area.** Biopython Shrake–Rupley SASA with a 1.4 Å probe: for each atom, ΔSASA = SASA(chain alone) − SASA(selected two-chain complex), at unchanged coordinates. Sum per residue and per chain, in Å². No chain is relaxed or cropped to a patch. Patch area sums contributions of its selected residues and may omit other buried residues. Values are reported separately for each side; they are not halved or called a single interface area. Carbon/sulfur contributions are also recorded. Sample-point approximation affects precision; repeat at higher point counts when small differences matter.
7. **Salt candidates.** Interchain Lys NZ/Arg NE,NH1,NH2 to Asp OD1,OD2/Glu OE1,OE2 at ≤4 Å. Each atom pair is retained; residue-pair summaries deduplicate by interaction type. Histidine and protonation-dependent assignments are omitted.
8. **Nonpolar candidates.** Interchain side-chain C/S atoms (excluding N,CA,C,O,OXT names) in residues AVILMFWYPC at ≤4.5 Å. This deliberately restricted definition is not a complete description of hydrophobic packing, aromatic stacking, or cation–π geometry.
9. **Steric-overlap candidates.** For C,N,O,S use radii 1.70,1.55,1.52,1.80 Å; flag radius sum minus distance >0.6 Å. SG–SG distances 1.8–2.2 Å are excluded as possible disulfides. This is a geometric screen, not MolProbity clashscore. Hydrogen bonds, protonation, other covalent links, and atom-specific radii are not resolved. Do not count overlaps as favorable interaction evidence.

Atom-level features are calculated for the selected chain pair. Patch features retain atom pairs whose two residue endpoints belong to the patch; they need not independently satisfy the Cα cutoff. Thus atom-contact counts and Cα-contact counts describe different geometries.

## Residue mapping and quantitative comparison

For each explicit chain pair, globally align CIF sequences with Biopython `PairwiseAligner`, BLOSUM62, gap-open −10 and gap-extension −0.5, including terminal gaps. Nongap columns define reference↔query residue correspondence. Report identities among mapped columns, both sequence coverages, and whether another equally optimal alignment exists. Version 1 chooses the first optimal alignment; it does not enumerate all tied alignments. Repeat the entire comparison with gap-open penalties −8 and −12. These checks probe sensitivity but do not establish biological correctness. [Biopython alignment documentation](https://biopython.org/docs/1.84/Tutorial/chapter_pairwise.html).

For each model pair, express query contacts in reference numbering. Exclude contacts with either endpoint unmapped and report their counts separately. Let R be eligible reference contacts, Q the mapped query contacts, and M = R ∩ Q:

- Precision = |M| / |Q|: fraction of query contacts reproduced in the reference.
- Recall = |M| / |R|: fraction of reference contacts reproduced in the query.
- Contact F1 = 2|M| / (|R| + |Q|).

If either set is empty, F1 is NA rather than an apparent perfect match or definitive mismatch. Precision/recall are NA only when their own denominator is zero. A nonempty disjoint comparison has F1=0. Contact F1 weights every pair equally; it is not a contact-quality or energy score. This is an exploratory overlap metric, not the published DockQ score.

Repeat set comparisons separately for salt candidates, nonpolar candidates, and steric-overlap candidates. Do not combine these with confidence or collapse favorable contacts and clashes into one score. The typed results are in `comparison.json` under each model pair's `typed_contacts`.

Rigidly fit all mapped Cα coordinates of the query's first selected chain onto the reference's first chain using an SVD least-squares proper rotation. Report anchor RMSD and partner RMSD after applying that same transform. Independently fit the partner to report its own RMSD; this helps expose fold/mapping differences that can inflate displacement. Three noncollinear mapped points are required for a defined fit. No iterative trimming or confidence weighting is used. Partner RMSD is in Å, not a pure rotation angle or binding score.

Within each ZIP, contact recurrence is the number of models containing the same numbered residue pair divided by model count. Repeated rows count once per model. These models are not independent experimental replicates, so recurrence is not a statistical confidence interval.

## Outputs and interpretation

- `report.md`: human-readable full/leaf patch descriptions, all model-pair comparisons, mapping sensitivity and limitations.
- `reference.json`, `query.json`: all residue descriptors, DSSP segments/status/version, full atom interactions, contacts/PAE, every patch node, parameters and source names.
- `comparison.json`: primary and alternative mappings, 25 comparisons per mapping when inputs have five models, typed comparisons, recurrence, input SHA-256 hashes and software/source versions.
- `mapping.tsv`: primary nongap residue correspondence for inspection.
- `model-comparisons.tsv`: scalar contact and coordinate metrics.
- `*-interactions.tsv`: individual atom-level geometric candidates.
- `*-recurrence.tsv`: per-contact recurrence.

The weakest link for divergent/truncated constructs is residue correspondence. Inspect mapping coverage, secondary structure, tied alignments, independent domain RMSDs, and sensitivity. A low F1 alone does not show loss of binding; mapping errors or different predicted poses can cause it. Validate with fuller constructs and curated sequence/structural correspondence, then compare the known-protein AF prediction against its experimental structure. A PDB validation stage and curated mapping import are not implemented here. Long-term calibration would require multiple known positive and negative interface pairs, not these two ensembles.

For the supplied MyD88/IRAK example, human A ↔ Nv A and human B ↔ Nv B are the user-provided biological assignments. The Nv IRAK construct boundary is provisional. Do not infer conservation or loss of interaction from this first report alone.

## Local DSSP installation used for the example

DSSP 4.6.1 was built from the official PDB-REDO/dssp source with CMake, FastFloat 8.0.2, and fetched libcifpp/Eigen dependencies. The executable is `tmp/dssp/bin/mkdssp`; this ignored local installation is not bundled with the project. The runtime dictionary directory contains the installed PDBx and DSSP dictionaries, the official ModelCIF dictionary, and libcifpp's `ccd-subset.cif` copied as `components.cif` (canonical amino acids sufficient for these inputs). Install a complete supported DSSP distribution for broader chemical-component coverage. Every residue in all ten supplied models received a DSSP assignment in the example run.
