# miRReg
A framework called miRReg to study miRNA regulatory network and module identification
This framework mainly include three steps: 
step 1: miRNA-mRNA regulatory network and module 
step 2: miRNA sponge network and module 
Step 3: function analysis 
In step 1, we get DE miRNA and mRNA by limma algorithm, then we use the LASSO algorithm get miRNA-mRNA regression matrix.Next, the DE mRNA and DE miRNA were combined with the LASSO regression matrix, input into the biclique algorithm to obtain the biclique module, and a regulatory network was constructed based on the module.
In step 2, the DE mRNA and DE miRNA were combined with miRTarBase to obtain miRNA and mRNA regulatory relationships.
Based on the sponge package miRspongeR to construct a regulatory network.
In step 3,function analysis by using the package ClusterProfiler.
