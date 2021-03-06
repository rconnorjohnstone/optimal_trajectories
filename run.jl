include("hw1/prob4.jl")

using .hw1_4

const μ = 3.986004e5 #km^3/s^2
r1 = 8000. #km
r2 = 10r1

# println(hw1_4.oe_to_xyz([8000, 0, 0, 0, 0, 0], μ))
hw1_4.plot_maneuver(r1, r2, π/2, μ, π/2)
# hw1_4.plot_maneuver(r1, r2, π/2, μ, 0.1)
