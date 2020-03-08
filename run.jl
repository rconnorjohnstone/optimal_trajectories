include("hw1/prob4.jl")
include("hw1/prob5.jl")

using .hw1_4
using .hw1_5

const μ = 3.986004e5 #km^3/s^2
r1 = 8000. #km
r2 = 10r1

# Plot the non-zoomed bi-elliptic plane change
# hw1_4.plot_maneuver(r1, r2, π/2, μ, π/2)

# Plot the zoomed bi-elliptic plane change
# hw1_4.plot_maneuver(r1, r2, π/2, μ, 0.1)
