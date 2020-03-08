module hw1_4

using LinearAlgebra, PlotlyJS, ORCA

include("../utilities/utilities.jl")
using .util

export total_manuever, plot_maneuver

function total_manuever(r1::Float64, 
                        r2::Float64,
                        i1::Float64, 
                        i2::Float64, 
                        i3::Float64, 
                        i4::Float64, 
                        μ::Float64)
  if i2 > i3
    return Inf
  else
    a1 = r1
    a2 = (r1+r2)/2
    e1 = 0
    e2 = (r2-r1)/(r2+r1)
    v1_pre = util.oe_to_xyz([a1, 0., i1, 0., 0., 0.], μ)[4:6]
    v1_post = util.oe_to_xyz([a2, e2, i2, 0., 0., 0.], μ)[4:6]
    v2_pre = util.oe_to_xyz([a2, e2, i2, 0., 0., π], μ)[4:6]
    v2_post = util.oe_to_xyz([a2, e2, i3, 0., 0., π], μ)[4:6]
    v3_pre = util.oe_to_xyz([a2, e2, i3, 0., 0., 0.], μ)[4:6]
    v3_post = util.oe_to_xyz([a1, e1, i4, 0., 0., 0.], μ)[4:6]
    Δv1 = norm(v1_post-v1_pre)
    Δv2 = norm(v2_post-v2_pre)
    Δv3 = norm(v3_post-v3_pre)
    return Δv1 + Δv2 + Δv3
  end
end

function plot_maneuver(r1::Float64, r2::Float64, Δ_i::Float64,  μ::Float64, lim::Float64)
  N = 200
  node1_inc = node3_inc = collect(range(0, length=N, stop=lim))
  z = zeros(Float64, N, N)
  base = zeros(Float64, N, N, 2)
  for i in 1:N
    for j in 1:N
      z[i,j] = total_manuever(r1, r2, 0., node1_inc[i], Δ_i-node3_inc[j], Δ_i, μ) 
      base[i,j,:] = [node1_inc[i], node3_inc[j]]
    end
  end

  data = contour(;z=z, x=node1_inc, y=node3_inc, 
                 name="ΔV", contours=attr(size=0.2),
                 colorbar=attr(title=attr(text="ΔV (km/s)")))
  if lim < Δ_i
    layout = Layout(
      ;title="Total ΔV vs. Δi₁ and Δi₃ for a 90° Bi-Elliptic Transfer (zoomed)",
      width=920,
      height=920,
      xaxis=attr(title=attr(text="Δi₁ (radians)")),
      yaxis=attr(title=attr(text="Δi₃ (radians)")))
    p = Plot(data, layout)
    savefig(p::Union{Plot,PlotlyJS.SyncPlot}, "hw1_4_zoomed.png")
  else
    layout = Layout(
      ;title="Total ΔV vs. Δi₁ and Δi₃ for a 90° Bi-Elliptic Transfer",
      width=920,
      height=920,
      xaxis=attr(title=attr(text="Δi₁ (radians)")),
      yaxis=attr(title=attr(text="Δi₃ (radians)")))
    p = Plot(data, layout)
    savefig(p::Union{Plot,PlotlyJS.SyncPlot}, "hw1_4.png")
  end
  display(plot(p))

  min, position = findmin(z)
  angles = base[position, :]
  println("Found min: ", min)
  println("At: ", 180/π*angles[1], "°, ", 180/π*angles[2], "°")
  while true
    println("Ctrl+C to close")
    sleep(10)
  end
end

end
