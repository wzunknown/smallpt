# smallpt
rewrite https://www.kevinbeason.com/smallpt/

## build
```bash
mkdir build && cd build
cmake .. && make
# or build yaml-cpp as shared library
cmake -DBUILD_SHARED_LIBS=ON .. && make
```

## usage
```bash
./mysmallpt path_to_config.yaml
```

50000spp:

<img src="img/image50000.png" width="512">
