#!/bin/bash
template='1s/.*/set output "^.png"/'
num=$1
template="${template/\^/$num}"
echo $template
sed -i "$template" ./show.gnu
gnuplot ./show.gnu
#sed -i '1s/.*/set output 'x.png'/' ./show.gnu
