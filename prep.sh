#!/bin/bash
# prep.sh
prep-gas thermally-perfect-air.inp tpair.lua
e4shared --prep --job=heatshield
