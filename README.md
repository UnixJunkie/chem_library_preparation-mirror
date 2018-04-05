# Suite of scripts for chemical library preparation -- JRM March 2018

Licence agreement: Feel free to reproduce it, modify it, add things to it
But if you do so, keep the original reference and most importantly, continue to share it with your colleagues

This folder contains a meta script (library_preparation.sh), launching commands for preparing a chemical library for docking.
It relies heavily on python scripts, and other software:

SOFTWARE REQUIREMENTS
generation of protomers: ChemAxon, with a working cxcalc binary
generation of conformers: RDKit > 2017 
proper mol2 file conversion: corina
generation of parameters: cgenff, rdkit (for absinth)
Python scripts are in ROOT/pythonscripts/
Reference free energy of solvation structure are in ROOT/otherfiles/

It also contains a part of my group meeting of the 4th of April, with benchmark/considerations in ROOT/Slides_librarypreparation_git.pdf

Initial commit 5th April 2018
This version should replace the version shared on Slack.
