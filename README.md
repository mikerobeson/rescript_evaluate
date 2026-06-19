:construction: 
**Note: this repo is still under develoopment.**
 :construction: 

# RESCRIPt-evaluate

RESCRIPt-eval is a supplementary python 3 package for the [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) plugin that supports cross validation for reference sequence and taxonomy databases. Thus, `rescript-eval` depends upon the [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) plugin within the [rachis-qiime2-*](https://library.qiime2.org/quickstart/qiime2) branch.

We are essentially returning some capabilities that have recently been removed from the main RESCRIPt plugin. Specifically, these are `evaluate-cross-validate` and `evaluate-fit-classifier`. We may add other supplimentary tools in the future. 


## Install RESCRIPt-evaluate

RESCRIPt is now installed as part of [QIIME 2](https://qiime2.org/) and you can install as outlined in the [QIIME 2 Documentation](https://docs.qiime2.org/), for `rachis-qiime2-2026.4` and later.

Once QIIME 2 is installed, you can do he following:

```
conda activate rachis-qiime2-2026.4
pip install git+https://github.com/mikerobeson/rescript_evaluate
qiime dev refresh-cache
```


### Read help documentation

To view a help menu for using `rescript-eval` via the QIIME 2 CLI:
```
qiime rescript-eval --help
```

## Tutorials that contain examples of some components of RESCRIPt-eval

*Note: the tutorials may not be updated to point to this repo for a while. When you see references to use RESCRIPt's actions `evaluate-cross-validate` and `evaluate-fit-classifier`, they will not work unless this repo is installed. In which case you'll run `qiime rescript-evaluate evaluate-cross-validate...`.*

These tutorials demonstrate some of the basic functionality of RESCRIPT and RESCRIPt-eval, via the q2CLI (QIIME 2 command-line interface):
- [General Overview and working with SILVA data](https://forum.qiime2.org/t/rescript-sequence-reference-database-management-tutorial/15494)
- [Getting sequences and taxonomy with get-ncbi-data](https://forum.qiime2.org/t/using-rescript-to-compile-an-sequence-databases-and-taxonomy-classifiers-from-ncbi-genbank/15947)
- [Building a COI database with BOLD sequences](https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129)
- [Using RESCRIPt's 'extract-seq-segments' to extract reference sequences without PCR primer pairs](https://forum.qiime2.org/t/using-rescripts-extract-seq-segments-to-extract-reference-sequences-without-pcr-primer-pairs/23618)
- [How to train a GTDB SSU classifier using RESCRIPt](https://forum.qiime2.org/t/how-to-train-a-gtdb-ssu-classifier-using-rescript/25725)
- [How to train a UNITE classifier using RESCRIPt](https://forum.qiime2.org/t/how-to-train-a-unite-classifier-using-rescript/28285)


## Getting Help

Problem? Suggestion? Technical errors and user support requests can be filed on the [QIIME 2 Forum](https://forum.qiime2.org/).


## Citation

If you use RESCRIPt in your research, please cite the following:

Michael S Robeson II, Devon R O'Rourke, Benjamin D Kaehler, Michal Ziemski, Matthew R Dillon, Jeffrey T Foster, Nicholas A Bokulich. (2021) *RESCRIPt: Reproducible sequence taxonomy reference database management*. PLoS Computational Biology 17 (11): e1009581. doi: [10.1371/journal.pcbi.1009581](http://dx.doi.org/10.1371/journal.pcbi.1009581).


## License

RESCRIPt is released under a BSD-3-Clause license. See LICENSE for more details.
