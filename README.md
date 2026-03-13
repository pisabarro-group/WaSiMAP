# WaSiMap

**WaSiMap** is a Python tool for identifying and mapping persistent water-sites from molecular dynamics (MD) simulations of biomolecular complexes.

The software detects regions where solvent molecules repeatedly occupy similar spatial positions across simulation frames and clusters them into **water-sites**, enabling comparative analysis across multiple simulations.

---
 
## Overview

Water molecules frequently mediate hydrogen bonding and stabilize interaction patterns that extend beyond direct residue–residue contacts in biomolecular recognition. Identifying solvent-mediated interactions is particularly important in protein–protein and protein–DNA interfaces, where conserved water-sites may contribute to binding specificity and stability.

WaSiMap analyzes explicitly solvated MD simulations to:

- Detect water molecules with increased interfacial persistence 
- Determine the heavy atoms contacted by these waters 
- Cluster water positions into spatial **water-site centroids** 
- Enable visual comparison of water-sites across **many** closely related MD simulations

The tool is designed to facilitate post-processing of MD-based HTP mutational series, allowing researchers to investigate how solvent-mediated interactions evolve across closely related molecular systems.

---

## Features

- Detection of persistent water molecules at biomolecular interfaces  
- Identification of contacted heavy atoms for each water molecule  
- Clustering of water coordinates into **water-site centroids**  
- Aggregation of results from multiple MD simulations  
- Generation of an interactive **HTML visualization interface**  
- Compatible with multiple MD trajectory formats (based on [MDTraj](https://github.com/mdtraj/mdtraj))

---

## Installation
The installation of WaSiMap is straightforward.

Install directly from PIP:

```bash
pip install wasimap
```

Or install from source:
```bash
git clone https://github.com/pisabarro-group/WaSiMap
cd wasimap
pip install .
```

## Example usage (with test data)
The following commands download three (3) MD simulations into the current folder (~320 Mb) and perform a WaSiMap analysis on all simulations.

```bash
wasimap --testdata
wasimap --gui
```

Upon completion, the WaSiMap viewer should open automatically in the default browser. Else, manually open http://localhost:8080/wasimap.html

To exit the program, use CTRL-C.


## Basic Usage
WaSiMap searches for MD trajectory files and their corresponding topology files in the folder where it is executed. Each trajectory must share the same base filename as its topology file. For example, an AMBER simulation named SIM1 should exist as SIM1.mdcrd and SIM1.prmtop. In contrast, H5 trajectories contain topology information as embedded metadata and therefore exist as a single file. WaSiMap automatically attempts to identify and pair trajectory files with their matching topologies.

The program supports MD trajectories readable by [MDTraj](https://github.com/mdtraj/mdtraj), including: .".xtc", ".trr", ".dcd", ".nc", ".h5", ".hdf5", ".lh5", ".binpos" - and also topology formats ".pdb", ".gro", ".prmtop", ".parm7", ".psf", and ".top".

Standard usage syntax:
```bash
wasimap --flag --flag ... [input]
```

For a typical automated execution, use:
```bash
wasimap --gui
```

To process a single trajectory, use:
```bash
wasimap --gui <<trajectory_filename>> (name only, no extension)
```

The command "wasimap" accepts several flags, described as follows:
```bash
 -c --cuttoff      H-bonding distance cutoff (default: 0.35 nm)
 -r --persistence  Minimum percentage of frames to consider water as resident (default: 5)
 --gui             Start the WaSiMap explorer after execution
 --onlygui         Do not process anything. Show GUI for existing results.
 --testdata        Download test MD simulations to the current folder
 -h --help         Displays help
 -v --version      Displays current version
```
## Navigating Results

WaSiMap automatically creates a folder called ./wasimap_outputs and stores:
- The 3D structure of the first frame of each trajectory, in PDB format
- A collection of water-site data, centroids and water ids in file wasimap.json
- A copy of the HTML viewer wasimap.html

A small web server binds to all interfaces (port 8080) upon completion, and the viewer is launched on the default browser. This web application embeds a 3D viewer (based on [NGL Viewer](https://github.com/nglviewer/ngl)) that superposes all structures. Water-sites are shown as spheres. The size of each sphere is proportional to the residence time of all waters in the water-site.

- Each structure is assigned a color for ease of identification
- Interfacial residues involved in water-sites are displayed as "ball+sticks"
- Molecular backbones are displayed as ribbons
- An intuitive menu (left panel) allows toggling structures and water-sites on/off
- The stage can be conveniently navigated with the mouse (see NGL documentation)
- Clicking on a water-site sphere shows the water-site table for the given structure.
- The table displays one water-site per row, and also its resident waters.
- Each water shows its internal AtomID
- Hovering the mouse over waters highlights its occurence around other residues.

## Citation
If you use this software in academic work, please cite us.
(DOI PENDING)

## Contributing
Contributions, bug reports, and feature requests are welcome.
Please open an issue or submit a pull request.

## Dependencies

WaSiMap makes use of these great open-source projects:

- [MDTraj](https://github.com/mdtraj/mdtraj) – Trajectory management library
- [NGL Viewer](https://github.com/nglviewer/ngl) – Molecular visualization