#! /home/sunyan/miniconda3/envs/r/bin/R

args <- commandArgs(T)


library(ape)
library(ggtree)
library(distanceR)
tree<-calc_tree(args[1])
write.tree(tree, file=args[2])
