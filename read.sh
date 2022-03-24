f=$1
n=${2:-"10"}
xxd $f | head -n $n
