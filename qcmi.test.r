#   qcmi.test.r
#   version 2021.12.24


rm(list=ls())

#   计算消耗时间
#   to calculate time cost
t0=Sys.time()  

#   加载相关R包
#   load R packages and data
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
wd="D:\\YANJIUSHENG\\method\\qcmi\\test"
save.wd="D:\\YANJIUSHENG\\method\\qcmi\\test\\result"

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
#   Filter the data  by occurrence and abundance
filteredps= filter_ps(ps= ps, occurrence.threshold=60 ,  abundance.threshold= 99)

#   通过WGCNA包计算相关性网络
#   The correlation network is calculated using the WGCNA package
result_net= cal_network(ps = filteredps, N = 100,r.threshold=0.6, p.threshold=0.05,method ="spearman") #    通常N设置为0，为计算速度在此设置前100个otu进行计算

#   将igraph格式结果文件赋值
#   Assign the igraph format result file
ig.wgcna=result_net$igraph

#   查看网络点和边的数量
#   Check the number of network nodes and edges
length(V(ig.wgcna))#number of nodes
length(E(ig.wgcna))#number of edges

#   保存edgelist格式文件
#   Save the edgelist file  
write.graph(ig.wgcna,file="edgelist.csv",format="ncol")

#   导入edgelist格式文件三联表
#   re-import the edgelist as triad table
edgelist=read.csv("edgelist.csv",header = FALSE,sep=' ')

#   得到过滤后的OTU表
#   Get filtered OTU table
ps_sub = filter_OTU_ps(ps = filteredps,Top = 100)
otu_table = as.data.frame(t(vegan_otu(ps_sub)))

#   扩散限制
#   classifies ecological associations of dispersal limitation
result_dl= assigned_process(link_table_row=edgelist, OTUabd=otu_table, p=0.05 , data= geo,cutoff=0, method="dl")

#   检测环境变量冗余
#   Detect redundancy of environment variables
library(Hmisc)
plot(varclus(as.matrix(env) ))

#   去除 rho>0.5 的变量
#   Remove variables with rho>0.5
env=env[,-3]
env=env[,-5]
plot(varclus(as.matrix(env) ))

#   环境过滤
#   classifies ecological associations of environmental filtering
result_ef= assigned_process(link_table_row=edgelist, OTUabd=otu_table, p=0.05 , data= env,cutoff=0, method="ef")

#   手动分类edge归属并导入得到新的网络
#   Manually classify edge ownership and import the new network
filteredlist <- read.csv("filteredlist.csv")
filteredlist=filteredlist[which( filteredlist[,4]=='bi'),]

#   转化为igraph格式并赋值weight
#   Convert to igraph format and assign weight  
ig=graph_from_edgelist(as.matrix(filteredlist[,1:2]) , directed = FALSE)
ig=set_edge_attr(ig, 'weight', index = E(ig),  as.numeric(filteredlist[,3]) )

#   定量生物互作强度
#   Quantify the local intensity of microbial biotic interaction
result_bi= qcmi(igraph= ig, OTU= otu_table, pers.cutoff=0)

data=cbind(env,result_bi[[3]],result_bi[[4]])


#   检测生物互作与非生物环境因子对微生物多样性分布格局的影响，包括alpha多样性和beta多样性
#   The effects of biotic interactions and abiotic factors on microbial diversity patterns, including alpha diversity and beta diversity, were examined  
#alpha
library(MASS)
library(adespatial)
library(lmerTest)   
#data=scale(data,center=T,scale=T)
data=as.matrix(data)  
res_forward.sel=cal_alphacon(data,"forward.sel",ps)
res_olslm=cal_alphacon(data,"olslm",ps)

#beta
library(vegan)
data=cbind(env,result_bi[[3]],result_bi[[4]])
otu_abund=as.data.frame(t(vegan_otu(ps)))
res_mantel= cal_betacon(data,otu_abund,"spearman")

#   计算所用时间
#   Calculate total time
(t=format(Sys.time()-t0))

