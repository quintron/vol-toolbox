cd build
cmake --build . --config Release --target install
ctest -C Release
cd..
