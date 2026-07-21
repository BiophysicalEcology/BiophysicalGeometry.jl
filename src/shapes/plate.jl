"""
    Plate <: AbstractShape

A flat plate-shaped organism shape.
"""
mutable struct Plate{M,D,B,C} <: AbstractShape
    mass::M
    density::D
    b::B
    c::C
end

function geometry(shape::Plate, ::Naked)
    volume = shape.mass / shape.density
    length_skin = (volume * shape.b * shape.c)^(1 / 3)
    width_skin = length_skin / shape.b
    height_skin = length_skin / shape.c 
    total = surface_area(shape, length_skin, width_skin, height_skin)
    characteristic_dimension = volume^(1 / 3) # width_skin * 2
    return Geometry(volume, characteristic_dimension, (; length_skin, width_skin, height_skin), SurfaceAreas(; total))
end
function geometry(shape::Plate, fur::Fur)
    volume = shape.mass / shape.density
    length_skin = (volume * shape.b * shape.c)^(1 / 3)
    width_skin = length_skin / shape.b
    height_skin = length_skin / shape.c 
    length_fur = length_skin + fur.thickness * 2
    width_fur = width_skin + fur.thickness * 2
    height_fur = height_skin + fur.thickness * 2
    total = surface_area(shape, length_fur, width_fur, height_fur)
    skin = surface_area(shape, length_skin, width_skin, height_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    fat = 0.0u"m"
    characteristic_dimension = volume^(1 / 3) + fur.thickness # width_fur * 2
    return Geometry(volume, characteristic_dimension, (; length_skin, width_skin, height_skin, length_fur, width_fur, height_fur, fat), SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Plate, fat::Fat)
    volume = shape.mass / shape.density
    length_skin = (volume * shape.b * shape.c)^(1 / 3)
    width_skin = length_skin / shape.b
    height_skin = length_skin / shape.c 
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    width_flesh = (flesh_volume * shape.b * shape.c)^(1 / 3) / shape.b
    fat = (width_skin - width_flesh) / 2
    total = surface_area(shape, length_skin, width_skin, height_skin)
    characteristic_dimension = volume^(1 / 3) # width_skin * 2 
    return Geometry(volume, characteristic_dimension, (; length_skin, width_skin, height_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Plate, fur::Fur, fat::Fat)
    volume = shape.mass / shape.density
    length_skin = (volume * shape.b * shape.c)^(1 / 3)
    width_skin = length_skin / shape.b
    height_skin = length_skin / shape.c
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    width_flesh = (flesh_volume * shape.b * shape.c)^(1 / 3) / shape.b
    fat = (width_skin - width_flesh) / 2
    length_fur = length_skin + fur.thickness * 2
    width_fur = width_skin + fur.thickness * 2
    height_fur = height_skin + fur.thickness * 2
    total = surface_area(shape, length_fur, width_fur, height_fur)
    skin = surface_area(shape, length_skin, width_skin, height_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness #width_fur * 2 
    return Geometry(volume, characteristic_dimension, (; length_skin, width_skin, height_skin, length_fur, width_fur, height_fur, fat), SurfaceAreas(; total, skin, convection))
end

# Surface area

function surface_area(shape::Plate, body)
    length = body.geometry.length.length_skin
    width = body.geometry.length.width_skin
    height = body.geometry.length.height_skin
    surface_area(shape, length, width, height)
end
surface_area(shape::Plate, length, width, height) = length * width * 2 + length * height * 2 + width * height * 2

# Silhouette area

function silhouette_area(shape::Plate, insulation::Union{Naked,Fat}, body::AbstractBody)
    length = body.geometry.length.length_skin
    width = body.geometry.length.width_skin
    height = body.geometry.length.height_skin
    normal = max(length * width, length * height, height * width)
    parallel = min(length * width, length * height, height * width)
    return (; normal, parallel)
end
function silhouette_area(shape::Plate, insulation::Union{Fur,CompositeInsulation}, body::AbstractBody)
    length = body.geometry.length.length_fur
    width = body.geometry.length.width_fur
    height = body.geometry.length.height_fur
    normal = max(length * width, length * height, height * width)
    parallel = min(length * width, length * height, height * width)
    return (; normal, parallel)
end

# Radius

skin_radius(shape::Plate, insulation::AbstractInsulation, body) = body.geometry.length.width_skin / 2

# naked
insulation_radius(shape::Plate, insulation::Naked, body) = body.geometry.length.width_skin / 2
flesh_radius(shape::Plate, insulation::Naked, body) = body.geometry.length.width_skin / 2

# fur
insulation_radius(shape::Plate, insulation::Fur, body) = body.geometry.length.width_fur / 2
flesh_radius(shape::Plate, insulation::Fur, body) = body.geometry.length.width_skin / 2

# fat
insulation_radius(shape::Plate, insulation::Fat, body) = body.geometry.length.width_skin / 2
flesh_radius(shape::Plate, insulation::Fat, body) = body.geometry.length.width_skin / 2 - body.geometry.length.fat

# fur plus fat
insulation_radius(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width_fur / 2
flesh_radius(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width_skin / 2 - body.geometry.length.fat

# Composition

attachment_surfaces(::Plate) = (Top, Bottom, SideA, SideB, SideC, SideD)

# Outer (insulation-aware) dimensions (L, W, H). Insulation-dispatched.
outer_dims(sh::Plate, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::Plate, ::Union{Naked,Fat}, body::AbstractBody) =
    (L = body.geometry.length.length_skin,
     W = body.geometry.length.width_skin,
     H = body.geometry.length.height_skin)
outer_dims(::Plate, ::Fur, body::AbstractBody) =
    (L = body.geometry.length.length_fur,
     W = body.geometry.length.width_fur,
     H = body.geometry.length.height_fur)

# Skin-level dimensions — used for flesh-anchored attachment positions.
function _plate_skin(body::AbstractBody)
    gl = body.geometry.length
    (gl.length_skin, gl.width_skin, gl.height_skin)
end

function surface_area(sh::Plate, body::AbstractBody, ::Top)
    d = outer_dims(sh, body); d.L * d.W
end
surface_area(sh::Plate, body::AbstractBody, ::Bottom) = surface_area(sh, body, Top())
function surface_area(sh::Plate, body::AbstractBody, ::SideA)
    d = outer_dims(sh, body); d.W * d.H
end
surface_area(sh::Plate, body::AbstractBody, ::SideB) = surface_area(sh, body, SideA())
function surface_area(sh::Plate, body::AbstractBody, ::SideC)
    d = outer_dims(sh, body); d.L * d.H
end
surface_area(sh::Plate, body::AbstractBody, ::SideD) = surface_area(sh, body, SideC())

function validate_range(::Plate, body::AbstractBody, loc::Union{Top,Bottom})
    L, W, _ = _plate_skin(body)
    abs(loc.x) ≤ L/2 || error("$(typeof(loc)) x out of range ±$(L/2): $(loc.x)")
    abs(loc.y) ≤ W/2 || error("$(typeof(loc)) y out of range ±$(W/2): $(loc.y)")
end
function validate_range(::Plate, body::AbstractBody, loc::Union{SideA,SideB})
    _, W, H = _plate_skin(body)
    abs(loc.y) ≤ W/2 || error("$(typeof(loc)) y out of range ±$(W/2): $(loc.y)")
    abs(loc.z) ≤ H/2 || error("$(typeof(loc)) z out of range ±$(H/2): $(loc.z)")
end
function validate_range(::Plate, body::AbstractBody, loc::Union{SideC,SideD})
    L, _, H = _plate_skin(body)
    abs(loc.x) ≤ L/2 || error("$(typeof(loc)) x out of range ±$(L/2): $(loc.x)")
    abs(loc.z) ≤ H/2 || error("$(typeof(loc)) z out of range ±$(H/2): $(loc.z)")
end

function surface_point(::Plate, body::AbstractBody, loc::Top)
    _, _, H = _plate_skin(body); (loc.x, loc.y, H/2)
end
function surface_point(::Plate, body::AbstractBody, loc::Bottom)
    _, _, H = _plate_skin(body); (loc.x, loc.y, -H/2)
end
function surface_point(::Plate, body::AbstractBody, loc::SideA)
    L, _, _ = _plate_skin(body); (L/2, loc.y, loc.z)
end
function surface_point(::Plate, body::AbstractBody, loc::SideB)
    L, _, _ = _plate_skin(body); (-L/2, loc.y, loc.z)
end
function surface_point(::Plate, body::AbstractBody, loc::SideC)
    _, W, _ = _plate_skin(body); (loc.x, W/2, loc.z)
end
function surface_point(::Plate, body::AbstractBody, loc::SideD)
    _, W, _ = _plate_skin(body); (loc.x, -W/2, loc.z)
end

surface_normal(::Plate, ::AbstractBody, ::Top)    = ( 0.0, 0.0, 1.0)
surface_normal(::Plate, ::AbstractBody, ::Bottom) = ( 0.0, 0.0, -1.0)
surface_normal(::Plate, ::AbstractBody, ::SideA)  = ( 1.0, 0.0, 0.0)
surface_normal(::Plate, ::AbstractBody, ::SideB)  = (-1.0, 0.0, 0.0)
surface_normal(::Plate, ::AbstractBody, ::SideC)  = ( 0.0, 1.0, 0.0)
surface_normal(::Plate, ::AbstractBody, ::SideD)  = ( 0.0, -1.0, 0.0)

function surface_centroid(::Plate, body::AbstractBody, ::Top)
    L, _, H = _plate_skin(body); (zero(L), zero(L), H/2)
end
function surface_centroid(::Plate, body::AbstractBody, ::Bottom)
    L, _, H = _plate_skin(body); (zero(L), zero(L), -H/2)
end
function surface_centroid(::Plate, body::AbstractBody, ::SideA)
    L, _, _ = _plate_skin(body); (L/2, zero(L), zero(L))
end
function surface_centroid(::Plate, body::AbstractBody, ::SideB)
    L, _, _ = _plate_skin(body); (-L/2, zero(L), zero(L))
end
function surface_centroid(::Plate, body::AbstractBody, ::SideC)
    _, W, _ = _plate_skin(body); (zero(W), W/2, zero(W))
end
function surface_centroid(::Plate, body::AbstractBody, ::SideD)
    _, W, _ = _plate_skin(body); (zero(W), -W/2, zero(W))
end
surface_centroid_normal(::Plate, ::AbstractBody, ::Top)    = ( 0.0, 0.0, 1.0)
surface_centroid_normal(::Plate, ::AbstractBody, ::Bottom) = ( 0.0, 0.0, -1.0)
surface_centroid_normal(::Plate, ::AbstractBody, ::SideA)  = ( 1.0, 0.0, 0.0)
surface_centroid_normal(::Plate, ::AbstractBody, ::SideB)  = (-1.0, 0.0, 0.0)
surface_centroid_normal(::Plate, ::AbstractBody, ::SideC)  = ( 0.0, 1.0, 0.0)
surface_centroid_normal(::Plate, ::AbstractBody, ::SideD)  = ( 0.0, -1.0, 0.0)
