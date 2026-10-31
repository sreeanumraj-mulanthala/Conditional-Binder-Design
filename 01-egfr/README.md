# pH-Specific Cross-Reactive EGFR Binder
This design aims to bind to the extracellular domain of EGFR found in tumor cells in acidic environments (not at pH 7.4). It is also designed to bind to both mouse and human EGFR. It was submitted to Challenge 01 of the Anthropic x Adaptyv Protein Design Competition (Track 03), 2026.

## Research Question
Cancers of epithelial origin, including head and neck cancers and colorectal cancers, overexpress EGFR. Antibodies directed against EGFR, such as cetuximab and panitumumab, are useful in the treatment of these cancers. However, binding of these antibodies to EGFR on normal epithelial cells causes side effects. It was found that normal epithelial cells (i.e. healthy cells) have an extracellular pH of around 7.4, while solid tumor cells have an extracellular pH of nearly 6.5. Therefore, the designed binder should bind to EGFR at a pH of 6.5 and not at a pH of 7.4. This would help to act on the tumor cells and not on the healthy cells. 
The requirement of the challenge is to design a binder that recognizes EGFR at pH 6.5 and does not recognize it at pH 7.4. The binder should recognize mouse EGFR.

## Method
BindCraft was utilized in order to design the final binding protein.

## Epitope Selection
- Epitope recognition of cetuximab was determined using the crystal structure (PDB 6ARU) and 24 residues of EGFR were determined to be within 4.5 Angstrom of the Fab.
- 7 of the 24 residue contacts were different for mouse (Q01279) and human (P00533) EGFR.
- Four out of seven present at single-clusters spanning residues 465 to 473.
- Pairwise CB distances across the conserved 17 contacts, separating them into two lobes with about 15 to 20 Å apart.
- One lobe contains residues 382, 408, 409, 411, 412, 415 and 438.
- This lobe was fully conserved in the mouse variant which had a pairwise distance of less than 11 Å.
- It is also at least 12 Å apart from the closest divergent residue.
- The solvent exposed residues A408, A409, A411, and A412 were also defined as hotspots.

## Reason for selection of binders
- H409 is protonated at pH 6.5, and becomes deprotonated at pH 7.4.
- H409 is part of a conserved region of EGFR.
- At pH 6.5, H409 is protonated and solvent exposed. Its solvent accessible surface area is 147.8 Å². Of all residues in this region, H409 is the most exposed.
- There are no charged side chains within 10 Å of H409. So, we can expect the pKa of - H409 to be similar to that of histidine side chains, which is between 6.0 to 6.5.
- Thus, the goal was to choose binders that had a carboxylate, oriented to salt bridge the protonated H409.
- At pH 6.5, the protonated imidazole and carboxylate groups form a charged hydrogen bond, while at pH 7.4, the imidazole is mostly neutral and the hydrogen bond with the carboxylate is lost.
- Therefore, the interaction of the binder at the salt bridge would depend upon the system pH.
- Binding at lower pH as compared to normal system pH was hypothesized to be one of the binding mechanisms.

## Selection of pH-Sensitive Binders
- This trimmed target was used for 49 BindCraft trajectories. Twenty-nine reached the - ProteinMPNN redesign stage, generating structures. In total, 374 complexes were analysed, with 187 different designs.
- Of the 187 unique designs, 62 had a binder carboxylate group within 4.0 Angstroms of H409 across 15 of the 29 productive trajectories.
- The carboxylate group was not designed to bind, was able to obtain a carboxylate within 4 Å of the imidazole in 62 of 187 designs by setting H409 as a hotspot.
- The designs were ranked based on the following criteria:
  Distance between the carboxylate and H409.
  Closeness of neighboring negative charges to H409.
- Observation of binding contacts with human-mouse divergent EGFR residues. This was weighted using the BLOSUM62 score to define severity of amino acid substitution and changes in physical-chemical properties of the substitution.

## Submitted Designs
- Of the 14 of the 29 productive BindCraft trajectories, 20 designs were chosen.
- Designs ranged between 63-96 residues.
- The maximum pairwise sequence identity in the set is 84.6%.
- 18 designs of the 20 had a predicted salt bridge to H409, ranging between 2.52-3.93 Å.
- The 2 negative control designs, which were mechanism-negative, had predicted binding quality of i_pTM of 0.82 and 0.85, respectively. The negative controls had a carboxylate group 7.5 and 10.8 Å from H409, respectively.
- The controls were introduced to determine if the proposed binding was pH dependent.

## BindCraft Acceptance Analysis
- Of the 187 unique designs, 6 of the 62 designs (9.7%) with the proposed salt bridge, passed the BindCraft acceptance filters.
- 10 of the 48 designs (20.8%) without a carboxylate within 7 Å of H409 also passed the acceptance filters.
- While the acceptance of the proposed design was low, the Fisher's Exact test showed this was not statistically significant (p = 0.11).
- Thus there is no clear correlation that could be extrapolated between salt-bridge formation and BindCraft.

## Limitations
- This pH sensitive interaction was selected on the basis of the predicted geometry.
- In order to assess the protonation behavior of H409, no explicit pKa calculations or conditions of constant pH were conducted.
- One salt bridge contributes to ~ 1-2 kcal/mol, or 5-30 fold affinity shift.
- As described in the previous section, it is probable that the influence of one salt-bridge will be beneficial at acidic pH, but less so at neutral pH.
- All 20 designs made contact with at least one of the human-mouse divergent residues (418, 467 or 468). This may affect the cross-reactivity. The two substitutions (G471A and N473K) that alter the surface chemistry most are not contacted by any of the designs.

## Residue Numbering
- EGFR residue numbering is based on the mature sequence. No precursor is included (which contains a signal peptide of 24 residues).
- Sequence alignment confirmed that PDB 6ARU uses mature numbering (607/609 positions aligned).
- BindCraft has renamed the target chain from 1, and therefore residue positions were mapped to the ordinal sequence, to ensure consistency.

## Repository Structure

```
analysis/
  align_egfr.py                   # Human–mouse sequence alignment
  analyse_6aru.py                 # Cetuximab epitope from 6ARU coordinates
  epitope_analysis.py             # Epitope, SASA and H409 analysis
  make_target_pdb.py              # EGFR target preparation
  rank_designs.py                 # Design ranking and conservation
  select_submission2.py           # Design selection
  build_submission.py             # Submission CSV assembly
target/
   6aru.cif                       # Source crystal structure
   EGFR_domainIII.pdb             # Trimmed design target
   *.fasta                        # Reference sequences
design/                           # 20 predicted complexes
results/                          # Per-design analysis results
```

## Reproducibility
- Requirements: Python 3, NumPy and Biopython.

- The analysis can be reproduced using:
```
python3 analysis/analyse_6aru.py       # Cetuximab epitope identification from PDB 6ARU
python3 analysis/align_egfr.py         # Human–mouse EGFR sequence conservation
python3 analysis/epitope_analysis.py   # SASA, H409 environment and epitope geometry
python3 analysis/make_target_pdb.py    # Preparation of the EGFR design target
```
- Three other scripts (rank_designs.py, select_submission2.py and build_submission.py) are included but require the original BindCraft trajectory files (377 MB), which are not included, due to their size.
- Processed outputs are available in the results/ directory. In particular, all_designs_h409_geometry.csv, contains the geometry analysis of all 374 predicted complexes generated by rank_designs.py.
- Designs were generated on Google Colab using BindCraft. Reproducibility of the design generation may vary slightly based on the GPU hardware and GPU driver versions.

## Software and AI Assistance
- Design generation: BindCraft using AlphaFold2 and ProteinMPNN with the four-stage multimer protocol and default advanced/filter settings.
- Data analysis: Python, Biopython and NumPy.
- AI assistance: The script and analysis were done using Python and Claude (Opus). Methodology was polished by Claude.
- The author made all research choices including selection of epitope, hypothesis on pH-sensitive binding and criteria for design selection.

## References
1. Pacesa, M., Nickel, L., Schellhaas, C. et al. (2025). One-shot design of functional protein binders with BindCraft. Nature, 646, 483–492. https://doi.org/10.1038/s41586-025-09429-6
2. Jumper, J., Evans, R., Pritzel, A. et al. (2021). Highly accurate protein structure prediction with AlphaFold. Nature, 596, 583–589. https://doi.org/10.1038/s41586-021-03819-2
3. Dauparas, J., Anishchenko, I., Bennett, N. et al. (2022). Robust deep learning–based protein sequence design using ProteinMPNN. Science, 378(6615), 49–56. https://doi.org/10.1126/science.add2187
4. Cotet, T.-S., Krawczuk, I., Stocco, F. et al. (2025). Crowdsourced protein design: Lessons from the Adaptyv EGFR binder competition. bioRxiv, 2025.04.17.648362. https://doi.org/10.1101/2025.04.17.648362
5. Tien, M. Z., Meyer, A. G., Sydykova, D. K., Spielman, S. J. & Wilke, C. O. (2013). Maximum allowed solvent accessibilites of residues in proteins. PLOS ONE, 8(11), e80635. https://doi.org/10.1371/journal.pone.0080635
6. Henikoff, S. & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks. Proceedings of the National Academy of Sciences, 89(22), 10915–10919. https://doi.org/10.1073/pnas.89.22.10915
7. Christie, M. & Christ, D. (2018). Structure of cetuximab Fab mutant in complex with EGFR extracellular domain (PDB ID: 6ARU). Protein Data Bank. https://doi.org/10.2210/pdb6ARU/pdb












