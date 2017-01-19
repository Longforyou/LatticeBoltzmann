#! /usr/bin/env julia

abstract FullPeriodicStreaming <: Streaming

# =========== Steaming

function compute_streaming(lbm::Lattice_Boltzmann_2D{Cells.D2Q9})

  # Distribution direction
  lbm.grid.f_prop[:, :, 1] = circshift(lbm.grid.f_temp[:, : ,1],[ 0  1])
  lbm.grid.f_prop[:, :, 2] = circshift(lbm.grid.f_temp[:, :, 2],[ 1  0])
  lbm.grid.f_prop[:, :, 3] = circshift(lbm.grid.f_temp[:, :, 3],[ 0 -1])
  lbm.grid.f_prop[:, :, 4] = circshift(lbm.grid.f_temp[:, :, 4],[-1  0])
  lbm.grid.f_prop[:, :, 5] = circshift(lbm.grid.f_temp[:, :, 5],[ 1  1])
  lbm.grid.f_prop[:, :, 6] = circshift(lbm.grid.f_temp[:, :, 6],[ 1 -1])
  lbm.grid.f_prop[:, :, 7] = circshift(lbm.grid.f_temp[:, :, 7],[-1 -1])
  lbm.grid.f_prop[:, :, 8] = circshift(lbm.grid.f_temp[:, :, 8],[-1  1])
  lbm.grid.f_prop[:, :, 9] = lbm.grid.f_temp[:, :, 9]
  
end