library(igraph)
info <- read.table("string_interactions.tsv",header = T,sep = "\t")# PPI nfile
links <- info[,c(1,2,ncol(info))]
colnames(links) <- c("from","to","weight")

nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
seed <- read.table("seeds.txt",header = T,sep = "\t",row.names = 1) #input your seed genes file
seed <- as.matrix(seed)
sameGene <- intersect(rownames(seed),nodes$gene)
Val <- NULL
for (id in nodes$gene) {
  if (id %in% sameGene) {
    Val <- c(Val,1)
  }else{
    Val <- c(Val,0)
  }
}
seedtab <- data.frame(gene=nodes$gene,value1=Val,value2=Val)
rownames(seedtab) <- seedtab$gene
seedtab <- seedtab[,-1]
write.table(seedtab,file = "seedtab.txt",quote = F,sep = "\t",col.names = F)


# create network graph
net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)

igraph::V(net)$deg <- igraph::degree(net)
igraph::V(net)$size <- igraph::degree(net)/5
igraph::E(net)$width <- igraph::E(net)$weight/10
#BiocManager::install("dnet")
library(dnet)

PTmatrix <- dRWR(g=net, normalise="laplacian",setSeeds = seedtab,restart=0.85,
                 parallel=TRUE,multicores = 12)

result <- as.matrix(PTmatrix)
rownames(result) <- nodes$gene
write.table(result,file="affinity_score.txt",quote = F,sep = "\t",col.names = F)
