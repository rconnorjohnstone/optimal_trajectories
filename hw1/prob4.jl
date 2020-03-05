module hw1_4

using LinearAlgebra, PlotlyJS, ORCA

export total_manuever, plot_maneuver

function oe_to_rθh(oe::Vector,μ::Real) :: Vector

  a,e,i,Ω,ω,ν = oe
  return [a*(1-e^2)/(1+e*cos(ν)),
  0,
  0,
  (μ/sqrt(μ*a*(1-e^2)))*e*sin(ν),
  (μ/sqrt(μ*a*(1-e^2)))*(1+e*cos(ν)),
  0]

end

function rθh_to_xyz(rθh_vec::Vector,oe::Vector)

  a,e,i,Ω,ω,ν = oe
  θ = ω+ν
  cΩ,sΩ,ci,si,cθ,sθ = cos(Ω),sin(Ω),cos(i),sin(i),cos(θ),sin(θ)
  DCM = [cΩ*cθ-sΩ*ci*sθ -cΩ*sθ-sΩ*ci*cθ sΩ*si;
  sΩ*cθ+cΩ*ci*sθ -sΩ*sθ+cΩ*ci*cθ -cΩ*si;
  si*sθ si*cθ ci]
  DCM = kron(Matrix(I,2,2),DCM)
  return DCM*rθh_vec

end

function oe_to_xyz(oe::Vector,μ::Real)

  return rθh_to_xyz(oe_to_rθh(oe,μ),oe)

end

function xyz_to_oe(cart_vec::Vector,μ::Real)

  r_xyz, v_xyz = cart_vec[1:3],cart_vec[4:6]
  r = norm(r_xyz)
  h_xyz = cross(r_xyz,v_xyz) #km^2/s
  h = norm(h_xyz) #km^2/s
  ξ = dot(v_xyz,v_xyz)/2 - μ/r #km^2 s^-2
  a = -μ/(2ξ) #km
  e = sqrt(1 + (2h^2*ξ)/μ^2)
  e_xyz = cross(v_xyz,h_xyz)/μ - r_xyz/r
  i = acos(h_xyz[3]/h) #rad
  n_xyz = cross([0,0,1],h_xyz)
  Ω = acos(dot(n_xyz,[1,0,0])/norm(n_xyz))
  ω = acos((dot(n_xyz,e_xyz)/(norm(n_xyz)*e)))
  ν = acos((dot(r_xyz,e_xyz))/(r*norm(e_xyz)))
  Ω = dot(n_xyz,[0,1,0]) > 0. ? Ω : -Ω
  ω = dot(e_xyz,[0,0,1]) > 0. ? ω : -ω
  ν = dot(r_xyz,v_xyz) > 0. ? ν : -ν
  return [a,e,i,Ω,ω,ν]

end

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
    v1_pre = oe_to_xyz([a1, 0., i1, 0., 0., 0.], μ)[4:6]
    v1_post = oe_to_xyz([a2, e2, i2, 0., 0., 0.], μ)[4:6]
    v2_pre = oe_to_xyz([a2, e2, i2, 0., 0., π], μ)[4:6]
    v2_post = oe_to_xyz([a2, e2, i3, 0., 0., π], μ)[4:6]
    v3_pre = oe_to_xyz([a2, e2, i3, 0., 0., 0.], μ)[4:6]
    v3_post = oe_to_xyz([a1, e1, i4, 0., 0., 0.], μ)[4:6]
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
