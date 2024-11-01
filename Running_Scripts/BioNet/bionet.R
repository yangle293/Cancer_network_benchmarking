library(graph,lib="~/rlibs")
library(RBGL,lib="~/rlibs")
library(BioNet,lib="~/rlibs")
library(R.matlab,lib="~/rlibs")
args = commandArgs(trailingOnly=TRUE)
ptm<-proc.time()
network_list = c("biogrid","irefindex18","reactome21","string")
data_list = c("corum_0","corum_1","corum_2","corum_3","corum_4","reactome_0","reactome_1","reactome_2","reactome_3","reactome_4")
permu_list = c(1,2,3,4,5,6,7,8,9,10)
beta_list = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11)
comb = expand.grid(network_list,data_list,permu_list,beta_list)

network = as.character(comb[strtoi(args[1]),1])
data = as.character(comb[strtoi(args[1]),2])
permu = comb[strtoi(args[1]),3]
beta = comb[strtoi(args[1]),4] 
print(network)
print(data)
print(permu)
print(beta)

#data = data_list[strtoi(args[1])]

#mutation_fileName = paste("./data/",data,"_score.mat",sep="")
network_edge_list_fileName = paste("./network/",network,"_edge_list",sep="")
network_node_list_fileName = paste("./network/",network,"_index_gene",sep="")
  
# load simulated mutation data
#mutation <- read.table(mutation_fileNam
mutation_filename = paste("./data/",data,"_p_beta_",beta,"_",permu,".txt",sep="")
gene_raw = read.table(mutation_filename)
gene <- as.character(gene_raw$V1)
pvalues = gene_raw$V2
pval = unlist(pvalues)
names(pval) = gene
#print(pval)
fb <- fitBumModel(pval,plot=FALSE)  
# load network
edge_list = read.table(network_edge_list_fileName,header=FALSE,sep="\t",quote="")
node_list = read.table(network_node_list_fileName,header=FALSE,sep="\t",quote="")
node_list = as.matrix(node_list)
nodes = unlist(node_list[,2])

edge_list = as.matrix(edge_list)
tmp = graphNEL(nodes = nodes,edgeL=list())
g = addEdge(nodes[edge_list[,1]],nodes[edge_list[,2]],tmp,edge_list[,3])
proc.time()-ptm
print("Finish graph building") 
for (j in 0.1)
{
  ptm<-proc.time()
  print("Scoring nodes...")
  score <- scoreNodes(g, fb, fdr=j)
  result = list()
  group = list()
  k=0
  while (1)
  {
    print(paste("Searching ",k+1,"-th subnetworks...",sep=""))
    k = k + 1
    edge_file = paste("./tmp/",network,data,"_p_beta",beta,"_",permu,"_",k,"_edge.txt",sep="")
    node_file = paste("./tmp/",network,data,"_p_beta",beta,"_",permu,"_",k,"_node.txt",sep="")
    output_file = paste("./result/",network,data,"_p_beta_",beta,"_",permu,"_",k,"_result.txt",sep="")
    writeHeinzEdges(network = g, file = edge_file,use.score = FALSE)
    writeHeinzNodes(network = g, file = node_file,node.scores = score)
    command <- paste("./heinz -e ", edge_file, " -n ", node_file, " -o ", output_file, sep = "")
    system(command)
    module <- readHeinzGraph(node.file = output_file, network = g)
        #plotModule(module,scores=score)
    result_node = nodes(module)
    if (length(result_node) == 0)
    {
      break
    }
    score[result_node] = -9999.0
  }
      
      
}


