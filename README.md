# Supplemental materials for manuscript "Isolating the Impact of Post-Translational Modification on MHC Peptide Binding and TCR Engagement"
Created on: 3/22/2023
Created by: Nathaniel Bloodworth
This file serves as a key for PTM .params files and rotamer libraries generated for the above manuscript.

## PTM descriptions, 3-letter codes, and generation method
ROSETTA patch file descriptions or the name of the custom ROSETTA .params files are provided.
If a .params files were used to model the PTM, it is provided in the supplementary materials along with the PDB rotamer library

|     PTM  Description        | PTM  |        Method          |
| ----------------------------| ---- | ---------------------- |
|N-terminal acetylated GLU    |  NGL | Patch: N_ACETYLATION   |
|Citrullinated ARG			  |  CAR | CAR.params             |    
|Hydroxylated PRO             |  HPR | Patch: HYDROXYLATION1  |
|Phosphorylated SER           |  SEP | Patch: PHOSPHORYLATION |
|1-methyl LYS                 |  1MK | Patch: METHYLATION     |         
|1,2-methyl LYS               |  2MK | Patch: DIMETHYLATION   |
|1,2,3-methyl LYS             |  3MK | Patch: TRIMETHYLATION  |
|Succinated LYS               |  SLY | SLY.params             |
|Acetylated LYS               |  ALY | Patch: ACETYLATION     |
|Biotinylated LYS             |  BIK | BIK.params             |

## Peptides sequences tested:

|          Peptide Description          |    Name   |     Sequence   | PDB Template |
| ------------------------------------- | --------- | -------------- | ------------ |
|COVID-19 Spike Protein derived peptide | sarsWT    | ESIVRFPNI      | 4PG9         |
|Spike protein: acetylated N-term       | sarsNAc   | [NGL]SIVRFPNI  | 4PG9         |
|Spike protein: R-5 citrullination      | sarsR5cit | ESIV[CAR]FPNI  | 4PG9         |
|Spike protein: P-7 hydroxylation       | sarsP7hy  | ESIVRF[HPR]NI  | 4PG9         |
|Ovalbumin derived peptide              | ovaWT     | SIINFEKL       | 3P9L         |
|Ovalbumin: S-1 phosphorylation         | ovaS1p    | [SEP]IINFEKL   | 3P9L         |
|Ovalbumin: K-7 methylation             | ovaK7m1   | SIINFE[1MK]L   | 3P9L         |
|Ovalbumin: K-7 dimethylation           | ovaK7m2   | SIINFE[2MK]L   | 3P9L         |
|Ovalbumin: K-7 trimethylation          | ovaK7m3   | SIINFE[3MK]L   | 3P9L         |
|Ovalbumin: K-7 scuccination            | ovaK7suc  | SIINFE[SLY]L   | 3P9L         |
|Ovalbumin: K-7 acetylation             | ovaK7Ac   | SIINFE[ALY]L   | 3P9L         |
|Ovalbumin: K-7 biotinylation           | ovaK7bio  | SIINFE[BIK]L   | 3P9L         |
|MBP peptide wild-type                  | mbpWT     | RTAHYGSL       | 3P4O         |
|MBP peptide: R-1 citrullination        | mbpR1cit  | [CAR]TAHYGSL   | 3P4O         |

## Protocol Capture

### Generating custom ROSETTA .params files and rotamer libraries
We use the biochemical library (BCL), developed by the Meiler lab and collaborators (1), to 
generate custom ROSETTA params file and rotamer libraries for specified PTMs above. The 
process is briefly described with example code below.

PATH_TO_BCL: path to the biochemical library executable
PTM_preoptimized.sdf: Initial PTM amino acid created in step 1
PTM.sdf: Optimized PTM amino acid created in step 2
PTM_rotamers.sdf: Rotamers file created in step 3
PATH_TO_ROSETTA: path to ROSETTA 3.13 installation /main folder

1. Starting with a 2D or 3D chemical structure in SDF format, use Pymol or a similar molecular
graphics software tool to add or remove additional molecules until the side chain structure
matches that of the desired PTM. Add additional N-terminus acetylated and C-teminus methylated
capping groups.
2. Generate an 3D structure of your PTM amino acid with optimized bond angles and steriochemistry using the BCL molecule:GenerateRosettaNCAAInstructions application
```
$ PATH_TO_BCL/bcl.exe molecule:GenerateRosettaNCAAInstructions \
    -input_filenames PTM_preoptimized.sdf \
    -generate_3D \
    -chirality L_AA \
    -explicit_aromaticity
```
3. Generate a rotamer library using the BCL molecule:ConformerGenerator application
```
$ PATH_TO_BCL/bcl.exe molecule:ConformerGenerator \
    -ensemble_filenames PTM.sdf \
    -conformation_comparer SymmetryRMSD 0.25 \
    -max_iterations 2000 \
    -top_models 200 \
    -cluster \
    -explicit_aromaticity \
    -conformers_single_file PTM_rotamers.sdf
```
4. Generate a ROSETTA .params file using molfile_to_params_polymer.py. NOTE: this script is located in PATH_TO_ROSETTA/demos/public/using_ncaas_protein_peptide_interface_design/HowToMakeResidueTypeParamFiles/scripts
```
$ python molfile_to_params_polymer.py \
    -n PTM \
    -p PTM \
    --polymer \
    --no-pdb \
    PTM_rotamers.sdf
```
5. Convert PTM_rotamers.sdf to PDB format using obabel or similar program
6. Rename heavy atoms in PTM_rotamers.pdb so they match ROSETTA atom names in PTM.params file
7. Add the following line to the end of the PTM.params file:
```
    PDB_ROTAMERS PTM_rotamers.pdb
```
8. You can now run the ROSETTA FlexPepDock application with the -extra_res and/or -extra_res_fa flags followed by the path to the PTM.params file. Ensure PTM_rotamers.pdb is located in the same path as PTM.params at runtime.

### Command line implementation of FlexPepDock

template.pdb: The MHC-I/peptide complex to use as a template
H-2Kb_alphaFold2.pdb: The idealized H-2Kb structure, generated by AlphaFold2, to serve as the receptor for the docking simulation
SIINFEKL: The peptide sequence of interest (single letter code, ova peptide used for this example)
PATH_TO_ROSETTA: path to ROSETTA 3.13 installation /main/source/bin folder

1. Select a template PDB peptide backbone structure based on sequence similarity
3. Using the ROSETTA Scripts mover SimpleThreadingMover, mutate the template residue to match the desired sequence. Add the PTM with the MutateResidue mover (for PTMs created using custom .params files) or the ModifyVariantType mover (for PTMs created with ROSETTA patch files). Add the H-2Kb MHC-I model alpha1 and alpha2 domains (produced by AlphaFold2, also available at https://alphafold.ebi.ac.uk/entry/P01901) with the AddChain mover. Run the scripted protocol (see Ref (2) for details on using ROSETTA Scripts).
```
    Example ROSETTA scripts file (rosettaScriptExample.xml)
    This example adds the custom PTM defined by the ROSETTA .params file PTM.params to the 7th 
    residue of the ova peptide.
    <ROSETTASCRIPTS>
        <RESIDUE_SELECTORS>
            <Index name="select_res_7" resnums="7B"/>
        </RESIDUE_SELECTORS>
        <MOVERS>
            <SimpleThreadingMover name="simpleThread" start_position="1B" skip_unknown_mutant="1" thread_sequence="SIINFEKL" />
            <MutateResidue name="addPTM" residue_selector="select_res_7" new_res="PTM" />
            <AddChain name="addAlphaFoldH2Kb" file_name="H-2Kb_alphaFold2.pdb" swap_chain_number="A"/>
            <FlexPepDock name="FlexPepPrepack" ppk_only="1"/>
        </MOVERS>
        <PROTOCOLS>
            <Add mover="simpleThread" />
            <Add mover="addPTM" />
            <Add mover="addAlphaFoldH2Kb"/>
            <Add mover="FlexPepPrepack"/>
        </PROTOCOLS>
    </ROSETTASCRIPTS>

    Command-line implementation:
    $ PATH_TO_ROSETTA/rosetta_scripts.linuxgccrelease \
        -in:file:s template.pdb \
        -parser:protocol rosettaScriptExample.xml \
        -extra_res_fa PTM.params \
        -ex1 -ex2aro -nstruct 1 \

    This command produces the prepacked structure named template_0001.pdb
```
4. Use FlexPepDock to generate 1000 structures, sorting the final decoys by reweighted_sc
```
    $ PATH_TO_ROSETTA/FlexPepDocking.linuxgccrelease \
        -in:file:s template_0001.pdb \
        -extra_res_fa PTM.params \
        -ex1 -ex2aro -nstruct 1000
        -out:file:silent example.silent

    This command generates the silent file example.silent. Extract the score data and store in the scorefile example.sc for further processing:

    $ grep '^SCORE' example.silent >example.sc
```

## References (also cited in main text)
1. Brown BP et al. Front Pharmacol. 2022;13: 833099.
2. Fleishman SJ et al. PLoS One. 2011;6: e20161.
