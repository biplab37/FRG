#########################################################################################
###    This code solves the exact FRG equations numerically for graphene near Dirac   ###
###    points. Here we have introduced a cutoff for the Bososnic momenta as well.     ###
#########################################################################################

include("GetVelEps.jl")
include("RGProcedure.jl")
include("Plotting.jl")

using .GetVelEps, .RGProcedure, .Plotting
using HDF5

## Initialisation

@show const m = 198 # number of cutoffs
@show const n = 209 # number of momenta

velocity = zeros(n,m)
dielectric = zeros(n,m)

@doc raw"""
    velocity_integrand(momentum, cutoff, phi, m, n, i)

This function returns the integrand of the FRG equation for the velocity renormlisation.
## Args
    dielectric (Array) : Array containing the dielectric values
    momentum (Float64) : momentum value
    cutoff   (FLoat64) : running cutoff
    phi      (Float64) : angular coordinate
    m          (Int64) : number of cutoffs
    n          (Int64) : number of momenta
    i          (Int64) : index for the running cutoff
"""
function velocity_integrand(dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)
    ## Theta function implementation with conditional
    if cos(phi)<=1 - 2*cutoff/momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi)*momentum

        eps1,eps2 = get_dielectric(dielectric, k1, k2, m, n, i)

        return 2.2*((momentum^2 - k1^2 + k2^2)/(momentum^2*eps1) + (momentum^2 + k1^2 - k2^2)/(momentum^2*eps2))/(2.0*pi*sqrt((k1+k2)^2 - momentum^2))
    end
end

@doc raw"""
    dielectric_integrand(momentum, cutoff, phi, m, n, i)

This function returns the integrand of the FRG equation for the dielectric function renormlisation.
## Args
    velocity   (Array) : The velocity array
    momentum (Float64) : momentum value
    cutoff   (FLoat64) : running cutoff
    phi      (Float64) : angular coordinate
    m          (Int64) : number of cutoffs
    n          (Int64) : number of momenta
    i          (Int64) : index for the running cutoff
"""
function dielectric_integrand(velocity::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)
    ## Theta function implementation
    if cos(phi)<=1 - 2*cutoff/momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi)*momentum

        vel1, vel2 = get_velocity(velocity, k1, k2, m, n, i)

        return 4.4*momentum*sin(phi)^2/(pi*(k1*vel1 + k2*vel2)*sqrt((k1+k2)^2 - momentum^2))
    end
end

## Boundary values initialisation
velocity[:,m] .= 1.0
dielectric[:,m] .= 1.0

## solving exact FRG using an user defined function in the module RGProcedure

rg_procedure(velocity,dielectric,velocity_integrand, dielectric_integrand ,m,n)


## Plots using user defined functions defined in the module Plotting

plot_velocity(velocity[:,1],"Renormalised_velocity.pdf")

plot_dielectric(dielectric[:,1],"Renormalised_dielectric.pdf")

## Save the data for future usage
# using JLD
# save("coupled.jld","velocity",velocity,"dielectric",dielectric)

# using HDF5
# h5write("coupled.h5","velocity",velocity)
# h5write("coupled.h5","dielectric",dielectric)
