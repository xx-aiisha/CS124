CS 1240 Programming Assignment 1
Author: Aisha Ahmed
Spring 2026

Files Included:
- randmst.cpp
- Makefile
- report.pdf (or report.tex)

Compilation:
Run:
    make randmst

This produces an executable named:
    randmst

Usage:
    ./randmst 0 numpoints numtrials dimension

Arguments:
    0            : required flag (unused in this implementation)
    numpoints    : number of vertices (n)
    numtrials    : number of independent trials
    dimension    : graph type

Dimension values:
    0 : Complete graph with edge weights uniform in [0,1]
    1 : Hypercube graph with edge weights uniform in [0,1]
    2 : Complete graph on random points in [0,1]^2
    3 : Complete graph on random points in [0,1]^3
    4 : Complete graph on random points in [0,1]^4

Output format:
    average numpoints numtrials dimension

Notes:
- For dimensions 0 and 1, edge weights are generated deterministically
  from a hash of (seed, u, v) to avoid storing all edges.
- Dimension 1 uses heap-based Prim’s algorithm (sparse graph).
- Other dimensions use array-based Prim’s algorithm (dense graph).
