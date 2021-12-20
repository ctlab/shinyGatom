### Help

This is Shiny GATOM:
a web-service for omics-based identification of regulated metabolic modules in atom transition networks. It is an updade of
[Shiny GAM](https://artyomovlab.wustl.edu/shiny/gam/) web-service (http://www.ncbi.nlm.nih.gov/pmc/articles/pmc4987878/) web-service.

Overall the analysis consists of two major steps:

* Step 1: Creating a metabolic network specific to the provided omics data.
* Step 2: Finding a connected subnetwork that contains the most significant changes.

*Reset all* button allows to start the analysis from the beginning, clearing the input data sets.


#### Example data

Three example datasets are embdeded in the web-service.

*Example DE for genes* and  *Example DE for metabolites* checkboxes load
transcriptional and metabolomic data respectively for M0 vs M1 macrophage activation comparison [Jha et al, Immunity 2015](http://dx.doi.org/10.1016/j.immuni.2015.02.005).
The raw files can be downloaded from [here](https://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.gene.de.tsv) and [here](https://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.met.de.tsv).

"Example DE for lipids" checkbox loads the example lipidomics data for WT vs the ΔSPT
comparison from [ST001127 dataset](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST001127). The raw file can be downloaded [here](https://artyomovlab.wustl.edu/publications/supp_materials/GATOM/ST27.lip.de.tsv).

We also provide example data with differential gene expression for human [MCF10A cell line treated with 2DG](https://artyomovlab.wustl.edu/publications/supp_materials/GAM/MCF10A.Ctrl.vs.2DG.gene.de.tsv)
(from [GSE59228](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59228)).


#### Network construction

First, slect an organism with which you will be working from the *Select an organism*
dropdown menu. Only reactions possible in the selected organism are used.
Currently supported organisms are Homo sapiens and Mus musculus.

The next step is to upload differential expression (DE) data
for genes and/or metabolites. GATOM can be run using
either gene DE data or metabolite DE data or both data types. Each DE dataset
must be in a separate text file (comma-, tab- and space- separated files are
supported, archived files are supported too). The first line of each file must
contain a header with the column names. The following columns are expected:

* "ID": RefSeq mRNA transcript ID, Entrez ID or symbol for genes
  and HMDB or KEGG ID for metbolites.
* "pval": Differential expression p-value (non-adjusted).
* "log2FC": Base 2 logarithm of the fold-change.
* "baseMean" (for genes): average expression level.
* "signal" (optional): ID of the measured entity, such as probe ID for gene expression microarrays and ion ID for metabolite data, if absent, signal column will be generated
based on unique "pval" and "log2FC" combinations.
* "signalRank" (for genes, optional): ranking of the genes by average expression, if absent, the values are generated based on "baseMean" column.

Any other columns will be copied to
a network as node or edge attributes.

After files are uploaded, a file summary is displayed. In the bottom of the summary rows
with IDs that were not mapped to the IDs used in network (e.g. Entrez for mouse) are shown. Verify that the files were parsed correctly.

![Summary of the example data](www/img/data_summary.png)

Next, a network can be constructed. There are thee major network types:
1) based on KEGG database, 2) based on Rhea database, 3) subset of Rhea reactions involving lipids. Additionally, two network topologies are supported: atom- and metabolite-level.
In atom-level topology, network nodes correspond to individual carbon atoms of compounds,
and an edge connects to atoms if there is a reaction with such atom transition.
In metabolite-level topology nodes correspond to metabolites, and two metabolites
are connected if there is a reaction where at least on carbon atom is transitioned between
the two metabolites.

Click *Step 1: Make network* to finish this step.

#### Finding a module

After making a network, you can find a connected module that contains the most
significantly changed genes and metabolites. Internally, this is done by first
scoring nodes and edges based on their p-values in such way that positive scores
correspond to significant p-values and negative scores correspond to
insignificant changes. The problem of finding a connected subgraph with
the maximum weight is solved with Virgo solver (https://github.com/ctlab/virgo-solver).

The size of the module can be controlled by changing "Number of positive genes" and "Number of positive metabolites" parameters. Increasing them will relax scoring thresholds and
will result in a bigger module.

Click *Step 2: Find module* button to find a module in the network. The module will
be shown on the right panel.

![Example of a module](www/img/module.png)

#### Graph legend

We use the following scheme:

* Red nodes and edges are up-regulated (*log2FC* > 0).
* Green nodes and edges are down-regulated (*log2FC* < 0).
* Blue nodes and edges don't have *log2FC* values.
* Bigger size of nodes and width of edges means lower p-values.


#### Saving the module

You can download the module in several formats by clicking corresponding links:
* SVG file contains the vector image of the module, as it appears on the web-site,
* XGMML file contains module that can be imported into Cytoscape (see *Importing to Cytoscape section below*).

#### Interpreting the results

The result of the analysis is a set of connected most regulated reactions. These reactions
are associated with strongly differentially changing enzymes and metabolites or
are closely connected to such reactions which implies regulated flux.
Inside the module one would typically consider small groups of particularly
significant reactions or individual reactions close
to the center of the module. Such groups or individual reactions are implied
to be of biological importance for the considered process.
These can lead to potential experimental designs including labeling,
media perturbation or gene knockout/knockdown experiments. For examples
see [Network integration of parallel metabolic and transcriptional data reveals metabolic modules that regulate macrophage polarization](http://www.ncbi.nlm.nih.gov/pubmed/25786174)
and [Mitochondrial Phosphoenolpyruvate Carboxykinase Regulates Metabolic Adaptation and Enables Glucose-Independent Tumor Growth](http://www.ncbi.nlm.nih.gov/pubmed/26474064)
papers.

#### Importing module into Cytoscape

XGMML file with module can be imported into Cytoscape for consecutive exploring
and processing. For the best experience we recommend to use GAM's
VizMap style for Cytoscape that can be downloaded
[here](https://artyomovlab.wustl.edu/publications/supp_materials/GAM/GAM_VizMap.xml).

#### Running the analysis in R

GATOM analysis can be done manually in R. For this check out [gatom](https://github.com/ctlab/gatom/) R package.



