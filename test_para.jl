JULIA_NUM_THREADS=4
import Pkg
Pkg.add("Base")
using Base.Threads
a = zeros(10)

@threads for i in 1:10
  # your code here
  a(i) = 1 + i
end
