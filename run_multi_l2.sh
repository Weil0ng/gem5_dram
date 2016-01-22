echo "Running $1 $2 $3 $4..."
echo "stats: $5"

#testing 8 core
time build/ALPHA/gem5.debug --stats-file=$5 configs/example/spec.py --num-cpus=8 --cpu-type=detailed --mem-type=DDR3_1600_x64 --mem-size=32GB --mem-channels=4 --caches --l1i_size=128kB --l1d_size=128kB --l1i_assoc=4 --l1d_assoc=4 --benchmark=$1\;$2\;$3\;$4\;$1\;$2\;$3\;$4 --standard-switch 1000000000 --maxinsts=100000000

