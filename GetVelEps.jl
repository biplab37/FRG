module GetVelEps

export get_velocity, get_dielectric

@doc raw"""
    get_velocity(velocity::Array{Float64,2}, k1::Float64, k2::Float64, m::Int64, n::Int64, i::Int64)

Returns the velocities at momenta k1 and k2 at the previous cutoff (indexed m-i+2).
"""
function get_velocity(velocity::Array{Float64,2}, k1::Float64, k2::Float64, m::Int64, n::Int64, i::Int64)

	index1 = Int64(round(k1*n))
    index2 = Int64(round(k2*n))

    vel1::Float64 = velocity[index1,m-i+2]

    if index2<n
        vel2::Float64 = velocity[index2,m-i+2]
    else
        vel2 = 1.0 # if the index goes out of the boundary take velocity to be v_F.
    end

    return vel1,vel2
end

@doc raw"""
    get_dielctric(dielctric::Array{Float64,2}, k1::Float64, k2::Float64, m::Int64, n::Int64, i::Int64)

Returns the dielectric functions at momenta k1 and k2 at the previous cutoff (indexed m-i+2).
"""
function get_dielectric(dielectric::Array{Float64,2}, k1::Float64, k2::Float64, m::Int64, n::Int64, i::Int64)

    index1 = Int64(round(k1*n))
    index2 = Int64(round(k2*n))

    eps1::Float64 = dielectric[index1,m-i+2]

    if index2<n
        eps2::Float64 = dielectric[index2,m-i+2]
    else
        eps2 = 1.0 # if the index goes out of the boundary take dielectric to be the free space one.
    end

    return eps1,eps2
end

end