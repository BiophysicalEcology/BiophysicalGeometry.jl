"""
    Ellipsoid <: AbstractShape

An ellipsoidal organism shape.
"""
mutable struct Ellipsoid{M,D,B,C} <: AbstractShape
    mass::M
    density::D
    b::B
    c::C
end

function geometry(shape::Ellipsoid, ::Naked)
    volume = shape.mass / shape.density
    b_semi_minor_skin = ((3 / 4) * volume / (π * shape.b)) ^ (1 / 3)
    c_semi_minor_skin = b_semi_minor_skin
    a_semi_major_skin = b_semi_minor_skin * shape.b
    e = ((ustrip(u"m", a_semi_major_skin) ^ 2 - ustrip(u"m", c_semi_minor_skin) ^ 2) ^ (1 / 2)) / ustrip(u"m", a_semi_major_skin)
    total = surface_area(shape, ustrip(u"m", a_semi_major_skin), ustrip(u"m", b_semi_minor_skin), ustrip(u"m", b_semi_minor_skin), e)
    characteristic_dimension = volume^(1 / 3) # b_semi_minor_skin * 2
    return Geometry(volume, characteristic_dimension, (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin), (; total))
end

function geometry(shape::Ellipsoid, fur::Fur)
    volume = shape.mass / shape.density
    b_semi_minor_skin = ((3 / 4) * volume / (π * shape.b)) ^ (1 / 3)
    c_semi_minor_skin = b_semi_minor_skin
    a_semi_major_skin = b_semi_minor_skin * shape.b
    e = ((ustrip(u"m", a_semi_major_skin) ^ 2 - ustrip(u"m", c_semi_minor_skin) ^ 2) ^ (1 / 2)) / ustrip(u"m", a_semi_major_skin)
    a_semi_major_fur = a_semi_major_skin + fur.thickness
    b_semi_minor_fur = b_semi_minor_skin + fur.thickness
    c_semi_minor_fur = c_semi_minor_skin + fur.thickness
    e = ((a_semi_major_skin ^ 2 - c_semi_minor_skin ^ 2) ^ (1 / 2)) / a_semi_major_skin
    e_fur = ((a_semi_major_fur ^ 2 - c_semi_minor_fur ^ 2) ^ (1 / 2)) / a_semi_major_fur
    total = surface_area(shape, ustrip(u"m", a_semi_major_fur),  ustrip(u"m", b_semi_minor_fur),  ustrip(u"m", c_semi_minor_fur),  e_fur)
    skin = surface_area(shape, ustrip(u"m", a_semi_major_skin), ustrip(u"m", b_semi_minor_skin), ustrip(u"m", c_semi_minor_skin), e)
    area_hair = hair_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair 
    characteristic_dimension = volume^(1 / 3) + fur.thickness # b_semi_minor_fur * 2
    return Geometry(volume, characteristic_dimension, (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin, a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur), (; total, skin, convection))
end

function geometry(shape::Ellipsoid, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    b_flesh = (((3 / 4) * flesh_volume) / (π * shape.b)) ^ (1 / 3)
    c_flesh = b_flesh # assuming c = b
    a_flesh = shape.b * b_flesh
    fat = prolate_fat_layer(
        ustrip(u"m^3", flesh_volume), 
        ustrip(u"m^3", fat_volume), 
        shape.b, 
        ustrip(u"m", b_flesh))
    if fat <= 0.0u"m"
        b_semi_minor_skin = (((3 / 4) * volume) / (π * shape.b)) ^ (1 / 3)
        c_semi_minor_skin = b_semi_minor_skin # assuming c = b
        a_semi_major_skin = shape.b * b_semi_minor_skin
    else
        a_semi_major_skin = a_flesh + fat
        b_semi_minor_skin = b_flesh + fat
        c_semi_minor_skin = c_flesh + fat
    end
    e = ((a_semi_major_skin ^ 2 - c_semi_minor_skin ^ 2) ^ (1 / 2)) / a_semi_major_skin
    total = surface_area(shape, ustrip(u"m", a_semi_major_skin), ustrip(u"m", b_semi_minor_skin), ustrip(u"m", c_semi_minor_skin), e)
    characteristic_dimension = volume^(1 / 3) # b_semi_minor_skin * 2
    return Geometry(volume, characteristic_dimension, (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin, fat), (; total))
end

function geometry(shape::Ellipsoid, fur::Fur, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    b_flesh = (((3 / 4) * flesh_volume) / (π * shape.b)) ^ (1 / 3)
    c_flesh = b_flesh # assuming c = b
    a_flesh = shape.b * b_flesh
    fat = prolate_fat_layer(
        ustrip(u"m^3", flesh_volume), 
        ustrip(u"m^3", fat_volume), 
        shape.b, 
        ustrip(u"m", b_flesh))
    if fat <= 0.0u"m"
        b_semi_minor_skin = (((3 / 4) * volume) / (π * shape.b)) ^ (1 / 3)
        c_semi_minor_skin = b_semi_minor_skin # assuming c = b
        a_semi_major_skin = shape.b * b_semi_minor_skin
    else
        a_semi_major_skin = a_flesh + fat
        b_semi_minor_skin = b_flesh + fat
        c_semi_minor_skin = c_flesh + fat
    end
    e = ((a_semi_major_skin ^ 2 - c_semi_minor_skin ^ 2) ^ (1 / 2)) / a_semi_major_skin
    a_semi_major_fur = a_semi_major_skin + fur.thickness
    b_semi_minor_fur = b_semi_minor_skin + fur.thickness
    c_semi_minor_fur = c_semi_minor_skin + fur.thickness
    e = ((a_semi_major_skin ^ 2 - c_semi_minor_skin ^ 2) ^ (1 / 2)) / a_semi_major_skin
    e_fur = ((a_semi_major_fur ^ 2 - c_semi_minor_fur ^ 2) ^ (1 / 2)) / a_semi_major_fur
    total = surface_area(shape, ustrip(u"m", a_semi_major_fur),  ustrip(u"m", b_semi_minor_fur),  ustrip(u"m", c_semi_minor_fur),  e_fur)
    skin = surface_area(shape, ustrip(u"m", a_semi_major_skin), ustrip(u"m", b_semi_minor_skin), ustrip(u"m", c_semi_minor_skin), e)
    area_hair = hair_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair 
    characteristic_dimension = volume^(1 / 3) + fur.thickness # b_semi_minor_fur * 2
    return Geometry(volume, characteristic_dimension, (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin, a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur, fat), (; total, skin, convection))
end

# fat thickness calculation

function prolate_fat_layer(
    flesh_volume,
    fat_volume,
    shape_b,
    semi_minor_flesh
)
    # Flesh is approximated as a prolate spheroid:
    # Volume = 4/3 π * A * B * C
    # C = B, A = shape_b * B
    # → B = ((3 * volume) / (4 * shape_b * π))^(1/3)
    # Fat thickness X is root of cubic: A X^3 + B X^2 + C X + D = 0
    A = 1.0
    B = shape_b * semi_minor_flesh + 2 * semi_minor_flesh
    C = 2 * shape_b * semi_minor_flesh^2 + semi_minor_flesh^2
    D = shape_b * semi_minor_flesh^3 - (( (fat_volume + flesh_volume) * 3.0 ) / (4.0 * π))

    # Components of cubic formula

    T1a = (-B)^3 / (27 * A^3)
    T1b = (B * C) / (6 * A^2)
    T1c = D / (2 * A)
    T1  = T1a + T1b - T1c

    T2a = T1^2
    T2b = ( (C / (3*A)) - (B^2) / (9 * A^2) )^3

    # Prevent sqrt of negative number
    T2 = (T2a + T2b >= 0) ? sqrt(T2a + T2b) : 0.0

    T3 = B / (3*A)

    # Cube roots with sign handling

    function signed_cuberoot(x)
        x < 0 ? -((-x)^(1/3)) : x^(1/3)
    end

    root1 = signed_cuberoot(T1 + T2)
    root2 = signed_cuberoot(T1 - T2)

    fat = (root1 + root2 - T3) * u"m"

    # If negative, not enough fat to cover spheroid
    return max(0.0u"m", fat)
end

# surface area functions

function surface_area(shape::Ellipsoid, body::AbstractBody)
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    c = body.geometry.length.c_semi_minor_skin
    return surface_area(shape, a, b, c)
end

function surface_area(shape::Ellipsoid, a, b, c)
    #e = ((a ^ 2 - c ^ 2) ^ 0.5 ) / a # eccentricity
    #2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)
    p =  1.6075
    return(4 * π * (((a ^ p * b ^ p + a ^ p * c ^ p + b ^ p * c ^ p)) / 3) ^ (1 / p))
end

function surface_area(shape::Ellipsoid, a, b, c, e)
    (2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)) * u"m^2"
end

# silhouette area functions

"""
    silhouette_area

Calculates the silhouette (projected) area of a prolate spheroid.
"""
function silhouette_area(shape::Ellipsoid, a, b, c, θ)
    a2 = cos(90) ^ 2 * (cos(θ) ^ 2 / a ^ 2 + sin(θ) ^ 2 / b ^ 2) + sin(90) ^ 2 / c ^ 2
    twohh = 2 * cos(90) * sin(90) * cos(θ) * (1 / b ^ 2 - 1 / a ^ 2)
    b2 = sin(θ) ^ 2 / a ^ 2 + cos(θ) ^ 2 / b ^ 2
    θ2 = 0.5 * atan(twohh, a2 - b2)
    sps = sin(θ2)
    cps = cos(θ2)
    a3 = cps * (a2 * cps + twohh * sps) + b2 * sps * sps
    b3 = sps * (a2 * sps - twohh * cps) + b2 * cps * cps
    semax1 = 1 / sqrt(a3)
    semax2 = 1 / sqrt(b3)
    π * semax1  * semax2
end

function silhouette_area(shape::Ellipsoid, insulation::Union{Naked,Fat}, body::AbstractBody)
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    c = body.geometry.length.c_semi_minor_skin
    normal = π * a * b
    parallel = π * b * c
    return (; normal, parallel)
end

function silhouette_area(shape::Ellipsoid, insulation::Union{Fur,CompositeInsulation}, body::AbstractBody)
    a = body.geometry.length.a_semi_major_fur
    b = body.geometry.length.b_semi_minor_fur
    c = body.geometry.length.c_semi_minor_fur
    normal = π * a * b
    parallel = π * b * c
    return (; normal, parallel)
end

function silhouette_area(shape::Ellipsoid, ::Naked, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    c = body.geometry.length.c_semi_minor_skin
    return silhouette_area(shape, a, b, c, θ)
end

function silhouette_area(shape::Ellipsoid, insulation::Fur, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major_fur
    b = body.geometry.length.b_semi_minor_fur
    c = body.geometry.length.c_semi_minor_fur
    return silhouette_area(shape, a, b, c, θ)
end

function silhouette_area(shape::Ellipsoid, insulation::Fat, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    c = body.geometry.length.c_semi_minor_skin
    return silhouette_area(shape, a, b, c, θ)
end

function silhouette_area(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major_fur
    b = body.geometry.length.b_semi_minor_fur
    c = body.geometry.length.c_semi_minor_fur
    return silhouette_area(shape, a, b, c, θ)
end

# area and radii functions

# naked
total_area(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.area.total
skin_area(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.area.total
evaporation_area(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.area.total

skin_radius(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.length.b_semi_minor_skin
insulation_radius(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.length.b_semi_minor_skin
flesh_radius(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.length.b_semi_minor_skin

# fur
total_area(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.area.total
skin_area(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.area.skin
evaporation_area(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.area.convection

skin_radius(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.length.b_semi_minor_skin
insulation_radius(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.length.b_semi_minor_fur
flesh_radius(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.length.b_semi_minor_skin

# fat
total_area(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.area.total
skin_area(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.area.total
evaporation_area(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.area.total

skin_radius(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.length.b_semi_minor_skin
insulation_radius(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.length.b_semi_minor_skin
flesh_radius(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.length.b_semi_minor_skin - body.geometry.length.fat

# fur and fat
total_area(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.area.total
skin_area(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.area.skin
evaporation_area(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.area.convection

skin_radius(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.length.b_semi_minor_skin
insulation_radius(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.length.b_semi_minor_fur
flesh_radius(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.length.b_semi_minor_skin - body.geometry.length.fat