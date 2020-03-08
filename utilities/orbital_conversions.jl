using LinearAlgebra

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
