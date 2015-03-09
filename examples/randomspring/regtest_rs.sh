#!/bin/bash

make CC="g++ -static" randomspring
./randomspring >/dev/null
cat sim*.pts | md5sum > __foo
cat __foo hash_check
diff __foo hash_check
rm __foo

