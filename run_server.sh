#!/bin/bash

export JULIA_NUM_THREADS=10
export JULIA_DEPOT_PATH="/rscratch/als217/.julia"

source bin/activate && jupyter lab --no-browser --port=8008
