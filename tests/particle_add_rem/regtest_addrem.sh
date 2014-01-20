#!/bin/bash

make CC="/data/misc/simulshare/compilers/gcc-4.8.2/bin/g++ -static" test_add_rem
./test_add_rem >/dev/null
cat sim*.pts sim*.bds | md5sum > __foo
cat __foo hash_check
diff __foo hash_check
rm __foo

