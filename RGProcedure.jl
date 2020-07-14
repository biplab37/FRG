@doc raw"""
# RGProcedure
This module provides the main function which solves the exact FRG equations.
"""
module RGProcedure

using QuadGK

export rg_procedure

@doc raw"""
    rg_procedure(velocity, dielectric, m, n)
This function solves the exact FRG equations for velocity and dielectric function.
"""
function rg_procedure(velocity, dielectric, velocity_integrand, dielectric_integrand, m::Int64, n::Int64)

	dcutoff::Float64 = 1.0/m

    for i in 2:m

        cutoff = Float64(m - i +1)/m

        for j in 1:n

            momentum = Float64(j)/n

            velocity_integrand_phi(phi) = velocity_integrand(dielectric,momentum,cutoff,phi,m,n,i)
            dielectric_integrand_phi(phi) = dielectric_integrand(velocity,momentum,cutoff,phi,m,n,i)

            ## Solving ODE's using Euler method
            velocity[j,m-i+1] = velocity[j,m-i+2] + dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
            dielectric[j,m-i+1] = dielectric[j,m-i+2] + dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
        end
    end
end

end