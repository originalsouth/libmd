tests:
	make test-regressiontester
	make test-boxshear
	make test-forcedipole
	make test-autodiff
	make test-geodesics
	make test-stresstensor
	make test-dissipativenetwork
	make test-dissipativespring
	make test-skewedbox

test-regressiontester:
	make -C ./regressiontester

test-boxshear:
	make -C ./boxshear

test-forcedipole:
	make -C ./forcedipole

test-autodiff:
	make -C ./autodiff

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
	make cleanboxshear
	make cleanforcedipole 
	make cleanautodiff
	make cleangeodesics
	make cleanstresstensor
	make cleandissipativespring
	make cleandissipativenetwork
	make cleanskewedbox

cleanregressiontester:
	make -C ./regressiontester clean

cleanboxshear:
	make -C ./boxshear clean

cleanforcedipole:
	make -C ./forcedipole clean

cleanautodiff:
	make -C ./autodiff clean 

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
