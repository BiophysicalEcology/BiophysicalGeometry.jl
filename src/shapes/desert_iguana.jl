"""
    DesertIguana <: AbstractShape

A lizard-shaped organism based on the desert iguana (Porter and Tracy 1984).

**Surface areas** (total, ventral, silhouette) are computed from empirical allometric
functions (Porter and Tracy 1984), *not* from ellipsoid or other geometric formulae.

**Body dimensions** for convection and conduction are derived by assuming a cylinder
with length L = 2D = 4R, so that volume = π R² L = 4π R³, giving R = (V/4π)^(1/3).
"""
struct DesertIguana{M,D} <: AbstractShape
    mass::M
    density::D
end

function geometry(shape::DesertIguana, ::Naked)
    volume = body_volume(shape)
    # Cylinder approximation for body dimensions: assume L=2D=4R, so vol=4π R³
    r_skin = cbrt(volume / (4 * π))
    # Surface areas from allometric functions (Porter and Tracy 1984)
    total, ventral = surface_area(shape)
    return Geometry(volume, (; r_skin), SurfaceAreas(; total, ventral))
end

# Surface area

function surface_area(shape::DesertIguana)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    total = Unitful.uconvert(u"m^2", (10.4713 * Unitful.ustrip(mass_g) ^ 0.688)u"cm^2")
    ventral = Unitful.uconvert(u"m^2", (0.425 * Unitful.ustrip(mass_g) ^ 0.85)u"cm^2")
    return (; total, ventral)
end

# Silhouette area

# Required by the generic silhouette_area(body::AbstractBody) dispatch in geometry.jl,
# which calls silhouette_area(shape, insulation, body) and expects a (; normal, parallel) NamedTuple.
function silhouette_area(shape::DesertIguana, ::Naked, ::AbstractBody)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    normal   = Unitful.uconvert(u"m^2", (3.798 * Unitful.ustrip(mass_g) ^ 0.683)u"cm^2")
    parallel = Unitful.uconvert(u"m^2", (0.694 * Unitful.ustrip(mass_g) ^ 0.743)u"cm^2")
    return (; normal, parallel)
end

function silhouette_area(shape::DesertIguana, ::NormalToSun)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    return Unitful.uconvert(u"m^2", (3.798 * Unitful.ustrip(mass_g) ^ 0.683)u"cm^2")
end
function silhouette_area(shape::DesertIguana, ::ParallelToSun)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    return Unitful.uconvert(u"m^2", (0.694 * Unitful.ustrip(mass_g) ^ 0.743)u"cm^2")
end
function silhouette_area(shape::DesertIguana, ::Intermediate)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    normal = Unitful.uconvert(u"m^2", (3.798 * Unitful.ustrip(mass_g) ^ 0.683)u"cm^2")
    parallel = Unitful.uconvert(u"m^2", (0.694 * Unitful.ustrip(mass_g) ^ 0.743)u"cm^2")
    return (normal + parallel) * 0.5
end

# Radius accessor — DesertIguana only supports Naked

_skin_radius(::DesertIguana, length) = length.r_skin
