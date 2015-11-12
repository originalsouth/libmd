tests:
	make test-regressiontester
	make test-forcedipole
	make test-autodiff
	make test-autodiff-yukawa-benchmark
	make test-geodesics
	make test-stresstensor
	make test-dissipativenetwork
	make test-dissipativespring
	make test-skewedbox

test-regressiontester:
	make -C ./regressiontester

test-forcedipole:
	make -C ./forcedipole

test-autodiff:
	make -C ./autodiff

test-autodiff-yukawa-benchmark:
	make -C ./autodiff-yukawa-benchmark

test-geodesics:
	make -C ./geodesics

test-stresstensor:
	make -C ./stresstensor
	
test-dissipativenetwork:
	make -C ./dissipativenetwork
	
test-dissipativespring:
	make -C ./dissipativespring
	
test-skewedbox:
	make -C ./skewedbox

clean:
	make cleanregressiontester
	make cleanforcedipole 
	make cleanautodiff
	make cleanautodiff-yukawa-benchmark
	make cleangeodesics
	make cleanstresstensor
	make cleandissipativespring
	make cleandissipativenetwork
	make cleanskewedbox

cleanregressiontester:
	make -C ./regressiontester clean

cleanforcedipole:
	make -C ./forcedipole clean

cleanautodiff:
	make -C ./autodiff clean 

cleanautodiff-yukawa-benchmark:
	make -C ./autodiff-yukawa-benchmark clean 

cleangeodesics:
	make -C ./geodesics clean

cleanstresstensor:
	make -C ./stresstensor clean

cleandissipativenetwork:
	make -C ./dissipativenetwork clean

cleandissipativespring:
	make -C ./dissipativespring clean

cleanskewedbox:
	make -C ./skewedbox clean
