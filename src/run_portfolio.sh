#!/bin/bash

inst_ ="port1"
solver_ = "cplex"

julia portfolio.jl >> ../report/out_${inst_}.txt
