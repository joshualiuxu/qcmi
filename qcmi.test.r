#   qcmi.test.r
#   version 2022.10.28
#	a method to quantify the effects of biotic filtering on the variations of microbial alpha- and beta-diversity

rm(list=ls())

#   计算消耗时间
#   to calculate time cost
t0=Sys.time()  

#   加载相关R包
#   load R packages
library("phyloseq")
library("ggClusterNet")
library("tidyverse")
library("WGCNA")
library("igraph")
library("Matrix")
library("vegan")
library("qcmi")

#   设置工作路径
#   the folder saving the input files
wd="D:\\YANJIUSHENG\\method\\qcmi\\bac"
save.wd="D:\\YANJIUSHENG\\method\\qcmi\\bac"

# 1 # 整理涉及到的数据，计算微生物共存网络
# 1 # Collate the data involved and calculate the microbial correlation network
#   导入数据
#   load  data
setwd(wd)
otu=read.table("otu_table.txt",header=TRUE, row.names=1)
tax=read.table("tax.txt",header=TRUE, row.names=1)
env=read.table("env.txt",header=TRUE, row.names=1)
geo=read.table("geo.txt",header=TRUE, row.names=1)

if(!dir.exists(save.wd)){dir.create(save.wd)}
setwd(save.wd)

#   数据格式转化为phyloseq
#   The data format is transformed into the phyloseq  
ps= trans_ps(otu_table = otu, taxa_table =  tax)

#   按照频率和丰度过滤数据集
#   Filter the data  by frequency and abundance
filteredps= filter_ps(ps= ps, occurrence.threshold=15,  abundance.threshold=100)
filteredps

#   通过不同的方法构建原始的生态网络
#   Through different methods to construct the original ecological network. (e.g., Pearson, Spearman, SpiecEasi)
result_net= cal_network(ps = filteredps,lambda.min.ratio=1e-2,nlambda=20,ncores=1,method ="spieceasi") 

#   评估网络可靠性
#   To evaluate the network stability (around 0.05 is acceptable)
getStability(result_net$result.se)
sum(getRefit(result_net$result.se))/2

#	若已使用其他网络推断方法构建网络，可灵活在此处导入igraph格式文件。
#	If you have used other network inference methods to build your network, you have the flexibility to import the igraph format file here.

#   将igraph格式结果文件赋值
#   Assign the igraph format result file
ig=result_net$igraph

#   查看网络点和边的数量
#   Check the number of network nodes and edges
length(V(ig))#number of nodes
length(E(ig))#number of edges

#	查看关联权重的分布情况
#	View the distribution of association weights
hist(result_net$elist[,3], main='', xlab='edge weights')

edgelist=result_net$elist

#   保存edgelist格式文件
#   Save the edgelist file  
write.csv(edgelist,file="edgelist.csv")

#   导入edgelist格式文件三联表
#   re-import the edgelist as triad table
edgelist=read.csv("edgelist.csv",header = FALSE,sep=' ')

#   得到过滤后的OTU表
#   Get filtered OTU table
otu_table = as.data.frame(t(vegan_otu(filteredps)))

#   扩散限制
#   Classify ecological associations of dispersal limitation
result_dl= assigned_process(link_table_row=edgelist, OTUabd=otu_table, p=0.05 , data= geo,cutoff=0, method="dl")

#	Optional
#   检测环境变量冗余
#   Detect redundancy of environment variables
library(Hmisc)
plot(varclus(as.matrix(env) ))

#   去除 rho>0.5 的变量
#   Remove variables with rho>0.5

#   环境过滤
#   Classify ecological associations of environmental filtering
result_pH= assigned_process(link_table_row=edgelist, OTUabd=otu_table, p=0.05 , data= env['pH'],cutoff=0, method="ef")
result_SM= assigned_process(link_table_row=edgelist, OTUabd=otu_table, p=0.05 , data= env['Salinity'],cutoff=0, method="ef")
result_Elevation= assigned_process(link_table_row=edgelist, OTUabd=otu_table, p=0.05 , data= env['Elevation'],cutoff=0, method="ef")
result_TN= assigned_process(link_table_row=edgelist, OTUabd=otu_table, p=0.05 , data= env['TN'],cutoff=0, method="ef")

#	保存结果
#	Save the results
write.csv(result_pH,"result_pH.csv")
write.csv(result_Salinity,"result_Salinity.csv")
write.csv(result_Elevation,"result_Elevation.csv")
write.csv(result_TN,"result_TN.csv")
write.csv(result_dl,"result_dl.csv")
#	注意: 海拔因素有的研究将其划分为综合环境，也有研究划分为空间限制，此处仅做展示。
#	Note: Some studies classify elevation as a comprehensive environment, while others classify it as dispersal limitation, which is only shown here.

library(dplyr)
total_link=row.names(edgelist)
ef_link=union(as.character(row.names(result_pH)),as.character(row.names(result_Salinity)))
ef_link=union(ef_link,as.character(row.names(result_Elevation)))
ef_link=union(ef_link,as.character(row.names(result_TN)))

bi_link=setdiff(total_link,ef_link)

efedge=edgelist[ef_link,]
write.csv(efedge,"efedge.csv")
biedge=edgelist[bi_link,]
write.csv(biedge,"effilteredlist.csv")

#   手动分类edge归属并导入得到新的网络
#   Manually classify edge ownership and import the new network
biedge <- read.csv("biedge.csv")
biedge=biedge[which( biedge[,4]=='bi'),]

#   转化为igraph格式并赋值weight
#   Convert to igraph format and assign weight  
ig.bi=graph_from_edgelist(as.matrix(biedge[,1:2]) , directed = FALSE)
ig.bi=set_edge_attr(ig.bi, 'weight', index = E(ig.bi),  as.numeric(biedge[,3]) )

#   定量生物互作强度
#   Quantify the local intensity of microbial biotic interaction
result_bi= qcmi(igraph= ig.bi, OTU= otu_table, pers.cutoff=0)

data=cbind(result_bi[[3]],result_bi[[4]],env)


#   检测生物互作与非生物环境因子对微生物多样性分布格局的影响，包括alpha多样性和beta多样性
#   The effects of biotic interactions and abiotic factors on microbial diversity patterns, including alpha diversity and beta diversity, were examined  
#alpha
library(MASS)
library(adespatial)
library(lmerTest)   

#data=scale(data,center=T,scale=T)
data=cbind(result_bi[[3]],result_bi[[4]],env)
colnames(data)=c("NC","PC",colnames(env))
data=as.matrix(data)  
res_forward.sel=cal_alphacon(data,"forward.sel",filteredps)
res_forward.sel

res_olslm=cal_alphacon(data,"olslm",filteredps)
res_olslm


#beta
library(vegan)

res_mantel= cal_betacon(data,otu_table,"spearman")

#   计算所用时间
#   Calculate total time
(t=format(Sys.time()-t0))

