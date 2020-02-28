module hw1_4

using LinearAlgebra, PlotlyJS, ORCA

export total_manuever, plot_maneuver

function apoapsis_change(r1::Float64, r2::Float64, μ::Float64)
  if r2 > r1
    circular_speed = √(μ/r1)
    eccentric_speed = √((2*μ)/(r1+r2)*(r2/r1))
    ΔV = abs(eccentric_speed - circular_speed)
  else
    circular_speed = √(μ/r2)
    eccentric_speed = √((2*μ)/(r1+r2)*(r1/r2))
    ΔV = abs(eccentric_speed - circular_speed)
  end
  return ΔV
end

function plane_change(r1::Float64, r2::Float64, μ::Float64, Δi::Float64)
  ra = r2 > r1 ? r2 : r1
  rp = r2 > r1 ? r1 : r2
  a = (ra+rp)/2
  e = (ra-rp)/(ra+rp)
  n = √(μ/a^3)
  return (2*sin(Δi/2)*√(1-e^2)*n*a)/(1+e)
end

function total_manuever(r1::Float64, r2::Float64,
                        Δ_i1::Float64, Δ_i3::Float64, Δ_i::Float64,  μ::Float64)
  if Δ_i1 + Δ_i3 > Δ_i
    return NaN
  else
    plane1 = plane_change(r1, r2, μ, Δ_i1)
    raise1 = apoapsis_change(r1, r2, μ)
    Δ_i2 = Δ_i - Δ_i1 - Δ_i3
    plane2 = plane_change(r1, r2, μ, Δ_i2)
    lower3 = apoapsis_change(r2, r1, μ)
    plane3 = plane_change(r1, r2, μ, Δ_i3)
    return plane1 + raise1 + plane2 + lower3 + plane3
  end
end

function plot_maneuver(r1::Float64, r2::Float64, Δ_i::Float64,  μ::Float64)
  N = 100
  node1_inc = node3_inc = collect(range(0, length=N, stop=Δ_i))
  z = zeros(Float64, N, N)
  for i in 1:N
    for j in 1:N
      z[i,j] = total_manuever(r1, r2, node1_inc[i], node3_inc[j], Δ_i, μ) 
    end
  end

  data = contour(;z=z, x=node1_inc, y=node3_inc)
  p = Plot(data)
  display(plot(p))
  savefig(p::Union{Plot,PlotlyJS.SyncPlot}, "hw1_4.png")
  while true
    println("Ctrl+C to close")
    sleep(10)
  end
end


end
