v0=0.1
for c in 0.1 1.0 10.0
do
./dissipativespring "$c" "$v0" > output_c"$c"
done
./dampedosc.gp
