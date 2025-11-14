# EVP
Plant based drug discovery pipeline through experimentally validated ligand binding.
This project links COCONUT, a natural products database where one can search from a multitude of species
and receive a list of plant derived compounds, to known ligand interactions captured on ChEMBL.
The purpose is to further elucidate potential medicinal properties in native plants
and identify alternative sources of drugs that are already commonplace within modern medicine.
Also, I am extremely new to programming, so any advice would be greatly appreciated.

How the Experimentally Validated Pipeline works:
1. Configure the plant you want to analyze in `start-cmds.R` (the default is `Digitalis purpurea`).
2. Extract COCONUT compounds and run the validation pipeline from an R session following the examples in
   `start-cmds.R`.