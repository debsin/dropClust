#!/bin/bash
# $0 is the script name, $1 id the first ARG, $2 is second...
#sh louvain_script <louvainpath> <inputfile> <outputfile>
wd=`pwd`

louvain_path="$1"
edgelistfile="$wd/$2"
outputfile=$wd"/$3"
echo $edgelistfile


#NewName="myFileIs${NAME}and that is all"
cd $louvain_path
echo "Converting to graph object ..."
./convert -i $edgelistfile -o  graph.bin

echo "Building Community tree..."
./louvain graph.bin -l -1 -q 0 > graph.tree
levels=`./hierarchy graph.tree | head -n 1|cut -d' ' -f4`
levels=`expr $levels - 1`

echo "Writing community labels..."
./hierarchy -l $levels graph.tree > $outputfile


rm graph.tree
rm graph.bin
cd $wd
