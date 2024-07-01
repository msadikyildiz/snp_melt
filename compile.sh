#!/usr/bin/bash
ENV="/home2/s414024/micromamba/envs/toprak-ngs-310"
x86_64-conda-linux-gnu-gcc -o snp_melt          \
    -I$ENV/include         \
    -L$ENV/lib             \
    -pthread -O3 -Wall -Wextra -std=c99 -lhts   \
    snp_melt.c
