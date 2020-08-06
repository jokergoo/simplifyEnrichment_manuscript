# simplifyEnrichment_manuscript

Scripts for the analysis in simplifyEnrichment manuscript:

- `collect_sim_mat.R`: It performes gene set enrichment analysis on [EBI Expression Atlas datasets](https://www.ebi.ac.uk/gxa/download) with various ontologies. It also calculates similarity matrices with different measurements on the significant terms.
- `run_example_random_GO.R`: It runs various clustering methods on the similarity matrices that were generated from random GO terms.
- `run_examples_EBI.R`: It runs various clustering methods on the similarity matrices from [EBI Expression Atlas datasets](https://www.ebi.ac.uk/gxa/download).
- `compare_semantic_and_overlap.R`: It compares clusterings from different similarity measurements, i.e., semantic similarity and gene overlap similarities.
- `test_partition_methods.R`: It tests the effect of partitionning methods in the binary cut clustering with rando GO datasets.
- `website.R`: It deploys all the results to GitHub Page.


Many analysis in the scripts were sent to the computing cluster with [the bsub package](https://github.com/jokergoo/bsub).
