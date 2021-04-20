# Things we won't do

* Load a structure with missing atoms
* Check for clashes or structural errors
* Support automatic filling in missing atoms/loops 

# User contracts

* Input data should have all the hydrogens. Complete molecules.
* Molecules can only become invalid through inplace modifications. Methods would not return invalid molecules.

# Guided example(s)


```python
# Load structures and perceive info
# This PDB file must have explicity hydrogens
protein = Molecule.from_pdb_file('protein.pdb')

# Pick some peptide C=O bond. The bond orders have been filled in for a library of substructures.
protein.atoms[9]
> <Atom with name O with element O>

protein.bonds[10]
> <Bond with bond_order 2 and atom1_index 9 and atom2_index 10>

print(protein.residues)
> ALA ALA LYS ALA

dna = Molecule.from_file('dna.sdf')
dna.bonds[10]
> <Bond with bond_order 1 and atom1_index 7 and atom2_index 11>

# Note that the SDF format doesn't contain residue names or numbers
dna.residues
> AttributeError
print(dna.atoms[30].metadata)
> {}
# This method searches a dictionary of substructures to assign common residue names and numbers
dna.perceive_residues(nucleotides=True)
dna.residues
> A T C G
print(dna.atoms[30].metadata)
> {'residue_num':2, 'residue_name':'T', 'chain':'A'}

# Where they attach
protein_attachment = protein.mdtraj_sel('resname LYS and name HD2')[0]
dna_attachment = dna.chemical_environment_matches('c1ccn([H:1])c1')[0]


# Merge the two
new_molecule = Molecule.merge_molecules(protein, dna, protein_attachment, dna_attachment, 
                                        bond_order=1, keep_metadata=False)

# Clear out old hierarchy definition and re-assign residue numbers and names
new_molecule.perceive_residues(amino_acids=True, 
                               nucleotides=True, 
                               clear_existing=True,
                               scheme='default')
new_molecule.residues
> ALA ALA UNK ALA T C G
```

 What if we wanted to keep the old residue names (and avoid residue `UNK`?)


```python
# Merge the two
new_molecule = Molecule.merge_molecules(protein, dna, protein_attachment, dna_attachment, 
                                        bond_order=1, keep_metadata=True)


# The residue numbers will now collide because we specified keep_metadata=True
new_molecule.residues
> ALA A ALA T LYS C ALA G

new_molecule.residue(5)  # Selects atoms in 6th residue in iter, but this thinks it's residue number 3
> <HierarchyElement with type residues with residue_num 3 residue_name C with chain A>


# Increment the residue numbers before the merge to keep them distinct
n_protein_res = len(protein.residues)
for residue in dna.residues:
    residue.residue_num += n_protein_res + 1
    
new_molecule = Molecule.merge_molecules(protein, dna, protein_attachment, dna_attachment, 
                                        bond_order=1, keep_metadata=True)   
    

new_molecule.residues
> ALA ALA LYS ALA A T C G
```

Why do I keep seeing the word "Hierarchy"?


```python
print(protein.hierarchy_schemes)
> [<HierarchyScheme with method_name='residues' uniqueness_criteria=['chain', 'residue_num', 'residue_name']>,
   <HierarchyScheme with method_name='chains' uniqueness_criteria=['chain']>]

protein.add_hierarchy_scheme(...)


```

# Reading from file

Reading biopolymer molecules from files, either in `sdf` or in `pdb` formats.

* from `sdf`: No residue information.
* from `pdb`: No detailed chemical information. Bond order, formal charges, stereochemistry, etc.

**Multiple molecules in single file**, so far. Molecule expected to be complete, including hydrogens.

* Load from SDF, perceive residues
* Load from PDB, fill in details according to known residues?


```python
protein = Molecule.from_file('some_mol.sdf')
protein.perceive_residues(pattern=['SMARTS_PATTERN', 'SMARTS_PATTERN2'])

dna_helix = Molecule.from_file('some_dna.sdf')
dna_helix.perceive_residues(nucleotides=True)
# Do we want to differentiate deoxy with ribonucleotides?

protein = Molecule.from_pdb_file('some_mol.pdb')

print(protein)
> Molecule with 5000 atoms

protein = Molecule.from_pdb_file('some_mol_missing_hs.pdb')
> IncompletePDBError: The Open Force Field Toolkit requires all 
> protons to be assigned before reading from PDB format. 
> (possibly more info here)

# Perceiving residues - Should be automatic for PDB
protein.perceive_residues()

```

Or we could have a special atom-typed object


```python
# Or a special atom-typed molecule
protein = PDBMolecule.from_pdb_file('some_mol.pdb')
print(protein)
> <PDBMolecule ...>
```

Or we could have a special method for atom-typed molecules/files


```python
# Or special method for PDB (atom-typed)
protein = Molecule.from_pdbfile('some_mol.pdb')
print(protein)
> <Molecule with X atoms>
```

# Reading from sequence

Reading biopolymer molecules from monomer/aa sequence string.

> - What about terminal caps? 

**Not much of a use for this**

    How to come up with conformers? Protonation? 
    How to match conformers to sequence?


```python
# from 3-letter code
protein = Molecule.from_sequence("AlaValGly", code='3-letter')

# from 1-letter code (default?)
protein = Molecule.from_sequence("AVQ")

print(protein)
> Molecule with 5000 atoms

# from 1-letter code (default?)
protein = Molecule.from_sequence("AVQZ")
> UnrecognizedAAError: Unrecognized amino acid in position 3.
# Analogous for 3-letter code. Could simply a ValueError do it?
```

# Histidines Protonation states

If we know the bonds orders we should be able to know the protonation state at a given pH and vice versa.


```python

```

# Flexible bookkeeping and grouping

Some (many?) of these things can be accomplished using some selection algebra language, as many visualizers and packages do (e.g. `name CA and resn LYS` ). We could eventually have a `to_mdtraj` as a way to use its selection language capabilities.


## Atoms


```python
protein = Molecule.from_file('protein.sdf')
protein.atoms
> <generator object at ... >

```


## Residues

Residues in molecule should be iterable.

**NOTE:** Do we need a `Residue` type? Would the `Molecule` type be sufficient? 


```python
# From SDF
protein = Molecule.from_file('protein.sdf')
protein.atoms[10].metadata['residue_name']
> KeyError: Whats a residue?
protein.perceive_residues()
protein.atoms[10].metadata['residue_name']
> ALA
protein.atoms[10].metadata['residue_num']
> 2

# From PDB
protein_from_pdb = Molecule.from_pdb_file('protein.pdb')
protein_from_pdb.atoms[10].metadata['residue_name']
> ALA


```


```python
protein.perceive_hierarchy(['residues', 'chains'])
type(protein.residues[10])
> <HierarchyElement of type residues with id ('A', 10, 'ALA') ...>
dir(protein.residues[10])
> [atoms, chain, residue_num, residue_name]
protein.residues[10].atoms
[150, 151, 152... 171]
protein.atoms[150].hier_info['residues']
> <HierarchyElement of type residues with id ('A', 10, 'ALA') ...>
protein.atoms[150].hier_info['chains']
> <HierarchyElement of type chain with id ('A',) ...>

# The user is allowed to shoot themself in the foot, and use the public API to fix it
protein.residues[10].atoms.append(900)
protein.residues[10].atoms
[150, 151, 152... 171, 900]
protein.atoms[900].hier_info['residues'] = protein.residues[10]

# manually modifying a residue can also get the backreferences broken
protein.residues[10].atoms = protein.residues[10].atoms[1:]
protein.atoms[150].hier_info['residues']
> <HierarchyElement of type residues with id ('A', 10, 'ALA') ...>
# fix
del protein.atoms[150].hier_info['residues']

# What if there's some special hand-modifications already performed, but the user has added something 
# new and wants that to have only the new parts to have their hierarchy percieved?
protein.perceive_hierarchy()
```

* Support iteration over `Biopolymer.residues` when residues are separated by
    * resname/residue + resnum/resid


```python
protein.residues
> <generator object residues at 0x7f7268acecf0>

residues_list = []
for residue in protein.residues:
    residues_list.append(residue)

print(residues_list)
> [<Residue ABC, ID1>, ..., <Residue XYZ, IDN>]
```


```python
print(protein.hierarchy_schemes)
> [<HierarchyScheme with method_name='residues' uniqueness_criteria=['chain', 'residue_num', 'residue_name']>,
   <HierarchyScheme with method_name='chains' uniqueness_criteria=['chain']>]

protein.add_hierarchy_scheme(...)


```

We should expose only the setters defined by the hierarchy scheme


```python
# Using the hierarchy scheme above
protein.perceive_hierarchy(scheme=['residues'])
protein.residues[10]
> <HierarchyElement of type "residues" with id ('A','11','ALA') with residue_num 11, chain A, residue_name ALA>
protein.residues[('A','11','ALA')]
> <HierarchyElement of type "residues" with id ('A','11','ALA') with residue_num 11, chain A, residue_name ALA>



protein.residues[10].residue_name
> ALA
protein.residues[10].residue_num
> 2
protein.residues[10].chain
> 0
or 
> "A"
protein.residues[10].residue_name = "LYS"
protein.residues[10].residue_name
> LYS
protein.residues[10].residue_num = 5
protein.residues[10].residue_num
> 5


# Now Change the hierarchy scheme -- What does it do with the previous hierarchy metadata?
protein.perceive_hierarchy(scheme=['chains'])
[*protein.chains]
> SomeSortOfDictyList[<HierarchyElement of type "chains" with id "A" with chain A>]
protein.chains[0]
> <HierarchyElement of type "chains" with id "A" with chain A>
protein.chains['A']
> <HierarchyElement of type "chains" with id "A" with chain A>

protein.residues[10].chain = 'B'
protein.residues[10].chain
> "B"

# Option 1 -- Let's do this one!
protein.chains['B']
> KeyError
protein.perceive_hierarchy(scheme=['chains'])
protein.chains['B']
> <HierarchyElement of type "chains" with id "B" with chain B>

# Option 2 -- This is probably bad, because it means the hierarchies are silently being re-percieved in the
# background, which could screw with things that the user manually set
protein.chains['B']
> <HierarchyElement of type "chains" with id "B" with chain B>



protein.chains['A'].chain = 'B'
[*protein.chains]
> SomeSortOfDictyList[<HierarchyElement of type "chains" with id "B" with chain B>]


```


```python
protein.residues[10]
> <HierarchyElement of type "residues" with id ('A','11','ALA') with residue_num 11, chain A, residue_name ALA>

protein.residues[('A','11','ALA')]
> <HierarchyElement of type "residues" with id ('A','11','ALA') with residue_num 11, chain A, residue_name ALA>

protein.residues[10].chain
> A
protein.residues[10].chain='Z'
protein.residues[10].chain
> A
protein.residues[-1].chain
> Z

protein.residues[('A','11','ALA')]
> KeyError

protein.residues[10]
> <HierarchyElement of type "residues" with id ('A','12','HIS') with residue_num 12, chain A, residue_name HIS>

protein.residues[('Z','11','ALA')]
> <HierarchyElement of type "residues" with id ('Z','11','ALA') with residue_num 11, chain Z, residue_name ALA>
protein.residues[-1]
> <HierarchyElement of type "residues" with id ('Z','11','ALA') with residue_num 11, chain Z, residue_name ALA>

```

Drawbacks of the above approach
* It is dangerous to allow integer-based access to hierarchy-lists, since they can go out of date
    * Maybe only expose iterator?
    * Should we expose subscriptable access at all? (allow `protein.chains[('A',)]`?). Or only allow acces via things like `protein.select_residues`?
* It is dangerous to allow atom-level metadata assignment though HierarchyElements, since they can change the groupings
* 

Possible future ideas?
* Allow hierarchies within hierarchies? `protein.chains['A'].residues`?

### Nonstandard residue handling

* Detect unnatural AA and assign parameters, gracefully handle backbone interface w/ natural AAs
* Modified AA becomes a single residue
* Modified AA becomes several residues (maybe one for each monosaccharide)

What if the lipid is attached in more than one site? 


```python
protein = Molecule.from_file('protein.sdf')
lipid = Molecule.from_file('lipid.sdf')

print(protein.atoms[10].metadata['resname'])
> KeyError: 

protein.perceive_residues(skip_hierarchy=True)
print(protein.atoms[10].metadata['resname'])
> ALA

protein.residues
> AttributeError

protein.perceive_hierarchy()  # Exposes iterables
protein.residues
> ALA ALA LYS GLU

protein_attachment_h = protein.mdtraj_sel('resname LYS and name HD2')[0]
lipid_attachment_h = lipid.mdtraj.sel('name H1')[0]
# This would end up with "invalid" mols -- get rid of it
# protein.delete_atom(protein_attachment_h)
# lipid.delete_atom(lipid_attachment_h)
#add a new bond between the atoms that had their Hs removed
#protein.add_bond(protein_attachment_h-1, lipid_attachment_h-1) 

new_molecule = Molecule.merge_molecules(protein, lipid, protein_attachment_h, lipid_attachment_h, 
                                        bond_order=1, keep_metadata=False)

# Scheme is how to handle unknown regions of the molecule
# e.g. The modified region becomes a single unknown residue
new_molecule.perceive_hierarchy(scheme='default')


print([*new_molecule.residues])
> ALA, ALA, UNK, GLU

new_molecule.clear_hierarchy()
# Different scheme. The modified region tries to retain the original res label
new_molecule.perceive_hierarchy(scheme='separate_ptms')
print([*new_molecule.residues])
> ALA, ALA, LYS, LIP, GLU



# change atom properties/data
protein.atoms[10].metadata['resname'] = 'GLU'
protein.atoms[10].update_metadata(key='resname', value='GLU')
```

Merge_molecules should optionally return a mapping from the original atom indices to the new ones


```python
len(protein.atoms)
> 10000
len(dna.atoms)
> 5000
protein_attachment
> <Atom with name H with element H>
dna_attachment
new_molecule, mol1_to_new_mol_idxs, mol2_to_new_mol_idx2 = Molecule.merge_molecules(protein, 
                                                                                    lipid, 
                                                                                    protein_attachment_h, 
                                                                                    lipid_attachment_h, 
                                                                                    bond_order=1, 
                                                                                    keep_metadata=False, 
                                                                                    return_atom_maps=True)
mol1_to_new_mol_idxs
> {0: 0, 1:1, ..., 9999:9999}
mol2_to_new_mol_idx2
> {0:, 1:, ..., 4999:}
```

Get residues given IDs


```python
protein.select_residues(resid=23)
> [<Residue GLY, 23>]

# Do we want to be able to get groups of residues? (what about getting atoms from these?)
protein.select_residues(resid=[23, 72, 12])
> [<Residue GLY, 23>, <Residue GLY, 23>, <Residue GLY, 23>]  # Maybe this should be a generator
```

Get the atoms from a given residue ID. Then the `Atom` must have some `resid` attribute?


```python
protein.select_residues(resid=23).atoms
> [<Atom name='' atomic number='6' resid='23'>, <Atom name='' atomic number='7' resid='23'>, 
   ..., <Atom name='' atomic number='1' resid='23'>]

# Do we want to be able to get groups of residues? (what about getting atoms from these?)
protein.select_residues(resid=[23, 72, 12]).atoms
> [
    [<Atom name='' atomic number='6' resid='23'>, ..., <Atom name='' atomic number='1' resid='23'>],
    [<Atom name='' atomic number='6' resid='72'>, ..., <Atom name='' atomic number='1' resid='72'>],
    [<Atom name='' atomic number='6' resid='12'>, ..., <Atom name='' atomic number='1' resid='12'>]
]
```

Get all residues matching a name


```python
protein.select_residues(name="LYS")
> [<Residue LYS, 3>, <Residue LYS, 21>, ..., <Residue LYS, 344>]
```

Renaming/transforming residues. 

> Should be just a metadata change


```python
protein.select_residues(resid=23)[0].name
> "GLY"
protein.select_residues(resid=23)[0].name = "LYS"
protein.select_residues(resid=23)[0].name
> "LYS"
protein.select_residues(resid=23)
> [<Residue LYS, 23>]
protein.select_residues(resid=23)[0].name = "WAT"
> 
protein.select_residues(resid=23)[0].name 
> "WAT"

# Do we want to be able to rename groups of residues?
# Then we need to match shapes/lengths -- is it getting out of control?
protein.select_residues(resid=[23, 72, 12]).name
> ["GLY", "ARG", "TYR"]
protein.select_residues(resid=[23, 72, 12]).name = ["LYS", "VAL", "CYS"]
protein.select_residues(resid=[23, 72, 12]).name
> ["LYS", "VAL", "CYS"]
```

We should be able to see numbers/ids and atoms


```python
protein.select_residues(resid=23)[0].
```

## Chains/Segments

We want to be able to group or cluster chains.

> Do we want a `Chain` class?

Assuming that you have multiple chains in your input files.

Chains should be iterable.


```python
protein = Molecule.from_file('protein.sdf')
protein.chains
> AttributeError: ...
protein.perceive_chains()
protein.chains
> ['A', 'B', ...,  'Z', '0', ..., '9']
protein.group_by(key='chains', values=['A', 'B'])
> AtomGroup
```


```python
protein = Molecule.from_file('protein.sdf')
protein.atoms[10].metadata['chain']
> 'B'
```

## Example. Protein - DNA merging

We expect coordinates of both molecules to be in the same frame of reference.

* Attach a cofactor like heme, which connects to 4 other residues


```python
# Load structures and perceive info
protein = Molecule.from_file('protein.sdf')
dna = Molecule.from_file('dna.sdf')
protein.perceive_residues()
protein.residues
> ALA ALA LYS ALA
dna.perceive_residues(nucleotides=True)
dna.residues
> A T C G
# Where they attach
protein_attachment = protein.mdtraj_sel('resname LYS and name HD2')[0]
dna_attachment = dna.chemical_environment_matches('c1ccn([H:1])c1')[0]
# Merge the two
new_molecule = Molecule.merge_molecules(protein, dna, protein_attachment, dna_attachment, 
                                        bond_order=1, keep_metadata=True)
                                        # What metadata to keep then?
new_molecule.residues
> ALA ALA LYS ALA A T C G
new_molecule.residue(5)  # Selects atoms in 6th residue in iter
> T

new_molecule.perceive_residues(amino_acids=True, 
                               nucleotides=True, 
                               clear_existing=True,
                               scheme='default')
new_molecule.residues
> ALA ALA UNK ALA T C G
or
> ALA ALA UNK T ALA C G
```

Can a atom be in more than one residue/HierarchyElement?

Will an atom have a pointer to its residue/HierarchyElement?
* DD -- In MDAnalysis, to keepthings organized and valid, we made the hierarchies very strict. Eg, every atom was a member of exactly one residue, residues were part of a chain. Then there were just a few tables, and it was trivial to check whether they were consistent/valid. 
* JW -- This may be too inflexible for our use cases and the variety of users we'll have. 

What if we want no residue information on DNA


```python
# Load structures and perceive info - No residues for DNA
protein = Molecule.from_file('protein.sdf')
dna = Molecule.from_file('dna.sdf')
protein.perceive_residues()
protein.residues
> ALA ALA LYS ALA
dna.residues
> AttributeError
# Where they attach
protein_attachment = protein.mdtraj_sel('resname LYS and name HD2')[0]
dna_attachment = dna.chemical_environment_matches('c1ccn([H:1])c1')[0]
# Merge the two
new_molecule = Molecule.merge_molecules(protein, dna, protein_attachment, dna_attachment, 
                                        bond_order=1, keep_metadata=True)
                                        # What metadata to keep then?
new_molecule.residues
> ALA ALA LYS ALA   # Note the difference here
new_molecule.residue(5)  # Selects atoms in 6th residue in iter
> IndexError: ...
new_molecule.residue(3)
> ALA

new_molecule.perceive_residues(amino_acids=True, 
                               nucleotides=True, 
                               clear_existing=True,
                               scheme='default')
new_molecule.residues
> ALA ALA UNK ALA T C G
or
> ALA ALA UNK T ALA C G

```

If we have multiple attachment sites do we just run `merge_molecules` multiple times?

# Modifications

## Adding/removing atoms

* Change protonation state of a residue


```python
mol = Molecule.from_file('protein.sdf')
mol.perceive_residues()
mol.residue(6)
> Residue ALA / AtomGroup with 15 atoms with index 90-105
new_idx = mol.add_atom(atomic_number=99, 
                       formal_charge=0, 
                       sterochemistry=None)
mol.add_bond(100, new_idx, bond_order=1, stereochemistry=None)
# EITHER
mol.remove_bond(100, 101)
mol.remove_atom(101)
# OR
mol.remove_atom(101)

mol.residue(6)
> Residue ALA / AtomGroup with 14 atoms

mol.perceive_residues()
mol.residue(6)
> Residue UNK / AtomGroup with 15 atoms


```

## Adding removing/residues

# Handling bonds


```python

```


```python

```

# Handling `pydantic` errors


```python
from pydantic import BaseModel, ValidationError, validator


class UserModel(BaseModel):
    name: str
    username: str
    password1: str
    password2: str

    @validator('name')
    def name_must_contain_space(cls, v):
        if ' ' not in v:
            raise ValueError('must contain a space')
        return v.title()

    @validator('password2')
    def passwords_match(cls, v, values, **kwargs):
        if 'password1' in values and v != values['password1']:
            raise ValueError('passwords do not match')
        return v

    @validator('username')
    def username_alphanumeric(cls, v):
        assert v.isalnum(), 'must be alphanumeric'
        return v
    
    @classmethod
    def from_sequence(cls, seq):
        name = ' '.join(seq)
        if 'GLY' in seq:
            raise ValueError("GLY doesn't exist")
        return cls(name=name, username='blah', password1='foo', password2='bar')


user = UserModel(
    name='samuel colvin',
    username='scolvin',
    password1='zxcvbn',
    password2='zxcvbn',
)
print(user)
#> name='Samuel Colvin' username='scolvin' password1='zxcvbn' password2='zxcvbn'

try:
    protein = UserModel.from_sequence("AlaValGLY")

    #UserModel(
    #    name='samuel',
    #    username='scolvin',
    #    password1='zxcvbn',
    #    password2='zxcvbn2',
    #)
except ValidationError as e:
    print(e)
    """
    2 validation errors for UserModel
    name
      must contain a space (type=value_error)
    password2
      passwords do not match (type=value_error)
    """
```

    name='Samuel Colvin' username='scolvin' password1='zxcvbn' password2='zxcvbn'



    ---------------------------------------------------------------

    ValueError                    Traceback (most recent call last)

    <ipython-input-5-62c3b29bec81> in <module>
         43 
         44 try:
    ---> 45     protein = UserModel.from_sequence("AlaValGLY")
         46 
         47     #UserModel(


    <ipython-input-5-62c3b29bec81> in from_sequence(cls, seq)
         29         name = ' '.join(seq)
         30         if 'GLY' in seq:
    ---> 31             raise ValueError("GLY doesn't exist")
         32         return cls(name=name, username='blah', password1='foo', password2='bar')
         33 


    ValueError: GLY doesn't exist



```python
# Pre-parameterization MoSDeF molecule
at_mol = AtomTypedMol()
# Implementation: https://github.com/mosdef-hub/gmso/blob/3ff3829cb4bc492b41e5e520d26d35c09c5338a4/gmso/core/atom.py#L42-L60
for mosdef_atom in mosdef_mol.atoms:
    at_mol.add_particle/atom(element=modef_atom.element, #(optional, atoms only)
                             name, #(required)
                             type, # (optional, sometimes overridden by manually specified mass or charge)
                             position=mosdef_atom.position, # (required)
                             mass, #(sometimes contained in type)
                             charge #(required, float, sometimes contained in type)) 
                            )
                             
for mosdef_bond in mosdef_mol.bond:
    at_mol.add_bond(mosdef_bond.particle1, mosdef_bond.particle2)
    
# MoSDeF atoms also have angle and dihedral classes, which mostly
# provide iterators, but ALSO have types. How will we support this
# in AtomTypedTopology?

```


```python
# Post-parameterization MoSDeF molecule
at_mol = AtomTypedMol()
for mosdef_atom in mosdef_mol.atoms:
    at_mol.add_particle/atom(element=modef_atom.element, (optional, atoms only)
                             type, # (optional, sometimes overridden by manually specified mass or charge)
                             position=mosdef_atom.position,
                             mass,
                             charge (float)) 

# JBG -- The final post-parameterization Topology has links to types, which
# conain the actual physical parameter values

# What's the goal of storing an atom typed molecule?
# Pre-parameterization: Keep the door open for having atom-type based
# parameter assignment functionality

# Post-parameterization: This will largely be the grouping and
# identifiers that users expect when they do their trajectory analysis

# JBG -- mBuild compounds contain hierarchies, but we don't yet carry
# that forward into GMSO/simulation output formats
# Eg image on https://mbuild.mosdef.org/en/stable/
#In the bottom figure, you can descend through hierarchy to the individual
# instances of the building blocks

# CQ -- Currently we force this kind of thing through ParmEd, and we
# artificially label the hierarchy in the way that ParmEd allows.
# This sometimes requires flattening/losing information.
# But after parameter assignment, we largely forget the rich description
# of the hierarchy, and
# #fix# the amount of information that is stored.

# CQ -- We're thinking about the idea of using "tags" instead of 
# rich hierarchy
```


```python
hier_level_names = ['complex','l2_id', 'l3_id', 'l4_id']
for hier_level_1 in mb_complex.hierarchy:
    for atom in hier_level_1:
        atom.metadata['complex'] = ...
    for hier_level_2 in hier_level_1.hierarchy:
        
```
