#!/bin/bash

make STD=c++0x shearxy
shearxy > __out_std
cat __out_std | md5sum > __out_hash
diff __out_hash hash_shearxy
rm __out_std __out_hash

