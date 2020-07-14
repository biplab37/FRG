#########################################################################################
###    This code solves the exact FRG equations numerically for graphene near Dirac   ###
###    points. Here we have introduced a cutoff for the Bososnic momenta as well.     ###
#########################################################################################

## Initialisation

const m::Int64 = 198 # number of cutoffs
const n::Int64 = 209 # number of momenta

velocity = zeros(n,m)
dielectric = zeros(n,m)

@doc raw"""
    velocity_integrand(momentum, cutoff, phi, m, n, i)

This function returns the integrand of the FRG equation for the velocity renormlisation.
## Args
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

        index1::Int64 = Int(round(k1*n))
        index2::Int64 = Int(round(k2*n))

        eps1::Float64 = dielectric[index1,m-i+2]

        if index2<n
            eps2::Float64 = dielectric[index2,m-i+2]
        else
            eps2::Float64 = 1.0 # if the index goes out of the boundary take dielectric to be free space one.
        end

        return 2.2*((momentum^2 - k1^2 + k2^2)/(momentum^2*eps1) + (momentum^2 + k1^2 - k2^2)/(momentum^2*eps2))/(2.0*pi*sqrt((k1+k2)^2 - momentum^2))
    end
end

@doc raw"""
    dielectric_integrand(momentum, cutoff, phi, m, n, i)

This function returns the integrand of the FRG equation for the dielectric function renormlisation.
## Args
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

        index1::Int64 = Int64(round(k1*n))
        index2::Int64 = Int64(round(k2*n))

        vel1::Float64 = velocity[index1,m-i+2]

        if index2<n
            vel2::Float64 = velocity[index2,m-i+2]
        else
            vel2::Float64 = 1.0
        end

        return 4.4*momentum*sin(phi)^2/(pi*(k1*vel1 + k2*vel2)*sqrt((k1+k2)^2 - momentum^2))
    end
end

## Boundary values initialisation
velocity[:,m] .= 1.0
dielectric[:,m] .= 1.0

## solving exact FRG using an user defined function
using .RGProcedure

rg_procedure(velocity,dielectric,velocity_integrand, dielectric_integrand ,m,n)

## Plots using a user defined module
using .Plotting

plot_velocity(velocity[:,1],"Renormalised_velocity.pdf")

plot_dielectric(dielectric[:,1],"Renormalised_dielectric.pdf")