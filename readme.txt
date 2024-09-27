BSDGNS --- Block Structural Discontinuous Galerkin for Navier-Stokes equations

分区结构网格求解器，采用LGL配点单元。
先从一维开始实现，高维复用一维。
高维使用分块结构四边形网格，在每个块上使用原先的算法。