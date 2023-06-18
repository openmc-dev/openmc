cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 16
sudo make install # update our old build with our new build