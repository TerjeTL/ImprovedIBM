#!/bin/bash

printf "-----------------------------------------\nRun Task: Release\n-----------------------------------------\n\n"

# cmake generate
printf "\n-----------------------------------------\nGenerating cmake ...\n-----------------------------------------\n\n"
cmake -GNinja \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -DCMAKE_C_COMPILER=clang \
  -DCMAKE_CXX_COMPILER=clang++ \
  -S . -B ./build/release

# copy compile_commands.json to /build/
printf "Moving compile_commands.json ...\n"
cp ./build/release/compile_commands.json ./build/compile_commands.json

# build
printf "\n-----------------------------------------\nBuilding ...\n-----------------------------------------\n\n"
cd ./build/release/
ninja
cd ./../../

# run
printf "\n-----------------------------------------\nRunning Application ...\n-----------------------------------------\n\n"
./build/release/ibm_app

