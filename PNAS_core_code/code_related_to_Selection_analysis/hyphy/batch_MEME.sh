#!/usr/bin/env bash
#hyphy-2.3.14

#args
codon_based_MSA=$1
tree=$2



inputRedirect = {};
inputRedirect["01"]="Universal"; // genetic code
inputRedirect["02"]="${codon_based_MSA}"; // codon data
inputRedirect["03"]="${tree}"; // tree
inputRedirect["04"]="All"; // Test for selection on a branch
inputRedirect["05"]=""; // complete selection

ExecuteAFile('./res/TemplateBatchFiles/SelectionAnalyses/MEME.bf', inputRedirect);

