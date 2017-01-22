#! /usr/bin/env julia

immutable FullPeriodicStreaming_2D <: Streaming
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}

    function FullPeriodicStreaming_2D(grid::Grid_2D)
        new(Array{Int64, 1}(1:grid.width),
            Array{Int64, 1}(1:grid.length))
    end

end

abstract InnerStreaming <: Streaming

# =========== Streaming
function compute_streaming(grid::Grid_2D{D2Q9},
                           stream::FullPeriodicStreaming_2D)

 / # Distribution direction
    grid.f_prop[stream.rows, stream.cols, 1] =
        circshift( grid.f_temp[:, : ,1],[ 0  1])
    grid.f_prop[stream.rows, stream.cols, 2] =
        circshift( grid.f_temp[stream.rows, stream.cols, 2],[ 1  0])
    grid.f_prop[stream.rows, stream.cols, 3] =
        circshift( grid.f_temp[stream.rows, stream.cols, 3],[ 0 -1])
    grid.f_prop[stream.rows, stream.cols, 4] =
        circshift( grid.f_temp[stream.rows, stream.cols, 4],[-1  0])
    grid.f_prop[stream.rows, stream.cols, 5] =
        circshift( grid.f_temp[stream.rows, stream.cols, 5],[ 1  1])
    grid.f_prop[stream.rows, stream.cols, 6] =
        circshift( grid.f_temp[stream.rows, stream.cols, 6],[ 1 -1])
    grid.f_prop[stream.rows, stream.cols, 7] =
        circshift( grid.f_temp[stream.rows, stream.cols, 7],[-1 -1])
    grid.f_prop[stream.rows, stream.cols, 8] =
        circshift( grid.f_temp[stream.rows, stream.cols, 8],[-1  1])
    grid.f_prop[stream.rows, stream.cols, 9] =
        grid.f_temp[stream.rows, stream.cols, 9]
  
end

