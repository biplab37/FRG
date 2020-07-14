@doc raw"""
# Plotting
===
    
This module contain the plotting functions which automatically plots the data and saves in a nice format.

## Examples
```julia-repl
julia> plot_velocity(data, "filename")

julia> plot_dielectric(data, "filename")
```
"""
module Plotting


using PyPlot # imports the matplotlib plotting routine

export plot_velocity,plot_dielectric

@doc raw"""
    plot_velocity(velocity[,name])

Plots the velocity as a function of momentum. Automatically labels the graph. If name is given then saves the figure with the given name.

## Example
```julia-repl
julia> plot_velocity(velocity,"renormalised_velocity.pdf")
```

"""
function plot_velocity(velocity,name="")
    
    figure()
    plot(range(0.0,stop=1.0,length=length(velocity)),velocity,label="velocity")
    title("Renormalised Velocity")
    xlabel(L"$k/\Lambda_0$")
    ylabel(L"$\dfrac{v_{\Lambda \to 0}(k)}{v_F}$")
    legend()
    if name != ""
        savefig(name)
    end
end

@doc raw"""
    plot_dielectric(dielectric[,name])

Plots the dielectric as a function of momentum and automatically labels the graph. If name is given then saves the figure with the given name.
    
## Example
```julia-repl
julia> plot_dielectric(dielectric,"renormalised_dielectric.pdf")
```
"""
function plot_dielectric(dielectric,name="")
    figure()
    plot(range(0.,stop=1.,length=length(dielectric)),dielectric,label="dielectric")
    title("Renormalised Dielectric Function")
    xlabel(L"$q/\Lambda_0$")
    ylabel(L"$\epsilon_{\Lambda \to 0}(q)$")
    legend()
    if name != ""
        savefig(name)
    end
end

end