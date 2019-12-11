# miRReg
We propose a framework called miRReg to study the miRNA regulation mechanism. For the direct miRNA regulation, we identify miRNA-mRNA regulatory network and further infer functional miRNA regulatory modules in TOF. Meanwhile, we infer mRNA related miRNA sponge interaction network and modules for the indirect miRNA regulation in TOF. The functional analysis shows that the identified direct or indirect miRNA regulatory networks and modules are closely associated with TOF disease. Moreover, several miRNA-mRNA regulatory relationships and miRNA sponge interactions are experimentally confirmed. We believe that miRReg can be also applied to study miRNA regulation in other human diseases. miRReg is released under the GPL-3.0 License.
This framework mainly include three steps:
 
Step 1:Identifying miRNA-mRNA regulatory network and module
Based on matched differential miRNA and mRNA expression profiles, we firstly use the LASSO regression method to construct miRNA-mRNA regression matrix. Then, based on the constructed miRNA-mRNA regression matrix, and miRNA and mRNA expression data, we infer miRNA-mRNA biclique modules using BiModuleWe use the default parameter settings of the BiModule method to infer miRNA-mRNA biclique modules in this work.For each miRNA-mRNA biclique module, all miRNAs are regarded as regulators of all mRNAs. By merging all of possible miRNA-mRNA interactions in each module, we remove duplicate miRNA-mRNA interactions to generate an integrated miRNA-mRNA regulatory network
 
Step 2:Identifying miRNA sponge interaction network and module, and functional analysis
To identify mRNA related miRNA sponge interaction network, we use the positive correlation (pc) method implemented in the miRspongeR R package. In the miRNA sponge interaction network, each mRNA-mRNA pair should have significant sharing of miRNAs (e.g. p-value < 0.05) and significant positive correlation (e.g. p-value < 0.05). To further understand the modularity of miRNA sponges in the network, we use the Markov Cluster Algorithm (MCL)implemented in the miRspongeR  package to infer miRNA sponge modules. The number of miRNA sponges in each module is at least 3.

Step 3: Function analysis
To understand whether the identified miRNA-related networks and modules are functional or not, it is necessary to perform functional analysis of them. For the GO and KEGG enrichment analysis, we use the clusterProfiler  R package in this work. Moreover, we use the experimentally validated miRNA-mRNA interactions from miRTarBase  and miRNA sponge interactions from miRSponge for validation.
More details can be found in the literature(miRReg: a framework for studying miRNA regulation in tetralogy of fallot)
