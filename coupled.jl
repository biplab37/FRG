#########################################################################################
###    This code solves the exact FRG equations numerically for graphene near Dirac   ###
###    points. Here we have introduced a cutoff for the Bososnic momenta as well.     ###
#########################################################################################

## package for numerical integration
using QuadGK

## Initialisation
const m = 198 # number of cutoffs
const n = 209 # number of momenta

const dcutoff = 1.0/m

velocity = zeros(n,m)
dielectric = zeros(n,m)

"""
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
function velocity_integrand(momentum,cutoff,phi,m,n,i)
    ## Theta function implementation with conditional
    if cos(phi)<=1 - 2*cutoff/momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi)*momentum

        index1 = Int(round(k1*n))
        index2 = Int(round(k2*n))

        eps1 = dielectric[index1,m-i+2]

        if index2<n
            eps2 = dielectric[index2,m-i+2]
        else
            eps2 = 1.0 # if the index goes out of the boundary take dielectric to be free space one.
        end

        return 2.2*((momentum^2 - k1^2 + k2^2)/(momentum^2*eps1) + (momentum^2 + k1^2 - k2^2)/(momentum^2*eps2))/(2.0*pi*sqrt((k1+k2)^2 - momentum^2))
    end
end

"""
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
function dielectric_integrand(momentum,cutoff,phi,m,n,i)
    ## Theta function implementation
    if cos(phi)<=1 - 2*cutoff/momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi)*momentum
        index1 = Int(round(k1*n))
        index2 = Int(round(k2*n))

        vel1 = velocity[index1,m-i+2]

        if index2<n
            vel2 = velocity[index2,m-i+2]
        else
            vel2 = 1.0
        end

        return 4.4*momentum*sin(phi)^2/(pi*(k1*vel1 + k2*vel2)*sqrt((k1+k2)^2 - momentum^2))
    end
end

## boundary values initialisation
velocity[:,m] .= 1.0
dielectric[:,m] .= 1.0

## The main for loops that implements the RG procedure
for i in 2:m

    cutoff = Float64(m - i +1)/m

    for j in 1:n

        momentum = Float64(j)/n

        velocity_integrand_phi(phi) = velocity_integrand(momentum,cutoff,phi,m,n,i)
        dielectric_integrand_phi(phi) = dielectric_integrand(momentum,cutoff,phi,m,n,i)

        ## Solving ODE's using Euler method
        velocity[j,m-i+1] = velocity[j,m-i+2] + dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
        dielectric[j,m-i+1] = dielectric[j,m-i+2] + dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
    end
end

## Plots using a user defined module

using .Plotting

plot_velocity(velocity[:,1],"Renormalised_velocity.pdf")

plot_dielectric(dielectric[:,1],"Renormalised_dielectric.pdf")