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

attachment_surfaces(::Plate) = (:top, :bottom, :side_a, :side_b, :side_c, :side_d)

function _plate_outer(body::AbstractBody)
    gl = body.geometry.length
    if haskey(gl, :length_fur)
        (gl.length_fur, gl.width_fur, gl.height_fur)
    else
        (gl.length_skin, gl.width_skin, gl.height_skin)
    end
end

# Skin-level dimensions — used for flesh-anchored attachment positions.
function _plate_skin(body::AbstractBody)
    gl = body.geometry.length
    (gl.length_skin, gl.width_skin, gl.height_skin)
end

function surface_area(::Plate, body::AbstractBody, ::Val{:top})
    L, W, _ = _plate_outer(body); L * W
end
surface_area(sh::Plate, body::AbstractBody, ::Val{:bottom}) = surface_area(sh, body, Val(:top))
function surface_area(::Plate, body::AbstractBody, ::Val{:side_a})
    _, W, H = _plate_outer(body); W * H
end
surface_area(sh::Plate, body::AbstractBody, ::Val{:side_b}) = surface_area(sh, body, Val(:side_a))
function surface_area(::Plate, body::AbstractBody, ::Val{:side_c})
    L, _, H = _plate_outer(body); L * H
end
surface_area(sh::Plate, body::AbstractBody, ::Val{:side_d}) = surface_area(sh, body, Val(:side_c))

function validate_position(::Plate, body::AbstractBody, ::Val{S}, pos) where {S}
    L, W, H = _plate_skin(body)
    if S === :top || S === :bottom
        issetequal(keys(pos), (:x, :y)) || error(":$S needs (x, y); got $(keys(pos))")
        abs(pos.x) ≤ L/2 || error(":$S x out of range ±$(L/2): $(pos.x)")
        abs(pos.y) ≤ W/2 || error(":$S y out of range ±$(W/2): $(pos.y)")
    elseif S === :side_a || S === :side_b
        issetequal(keys(pos), (:y, :z)) || error(":$S needs (y, z); got $(keys(pos))")
        abs(pos.y) ≤ W/2 || error(":$S y out of range ±$(W/2): $(pos.y)")
        abs(pos.z) ≤ H/2 || error(":$S z out of range ±$(H/2): $(pos.z)")
    elseif S === :side_c || S === :side_d
        issetequal(keys(pos), (:x, :z)) || error(":$S needs (x, z); got $(keys(pos))")
        abs(pos.x) ≤ L/2 || error(":$S x out of range ±$(L/2): $(pos.x)")
        abs(pos.z) ≤ H/2 || error(":$S z out of range ±$(H/2): $(pos.z)")
    else
        error("Plate has no surface :$S; valid: $(attachment_surfaces(shape(body)))")
    end
end

function surface_point(::Plate, body::AbstractBody, ::Val{:top}, pos)
    _, _, H = _plate_skin(body); (pos.x, pos.y, H/2)
end
function surface_point(::Plate, body::AbstractBody, ::Val{:bottom}, pos)
    _, _, H = _plate_skin(body); (pos.x, pos.y, -H/2)
end
function surface_point(::Plate, body::AbstractBody, ::Val{:side_a}, pos)
    L, _, _ = _plate_skin(body); (L/2, pos.y, pos.z)
end
function surface_point(::Plate, body::AbstractBody, ::Val{:side_b}, pos)
    L, _, _ = _plate_skin(body); (-L/2, pos.y, pos.z)
end
function surface_point(::Plate, body::AbstractBody, ::Val{:side_c}, pos)
    _, W, _ = _plate_skin(body); (pos.x, W/2, pos.z)
end
function surface_point(::Plate, body::AbstractBody, ::Val{:side_d}, pos)
    _, W, _ = _plate_skin(body); (pos.x, -W/2, pos.z)
end

surface_normal(::Plate, ::AbstractBody, ::Val{:top},    _) = ( 0.0,  0.0,  1.0)
surface_normal(::Plate, ::AbstractBody, ::Val{:bottom}, _) = ( 0.0,  0.0, -1.0)
surface_normal(::Plate, ::AbstractBody, ::Val{:side_a}, _) = ( 1.0,  0.0,  0.0)
surface_normal(::Plate, ::AbstractBody, ::Val{:side_b}, _) = (-1.0,  0.0,  0.0)
surface_normal(::Plate, ::AbstractBody, ::Val{:side_c}, _) = ( 0.0,  1.0,  0.0)
surface_normal(::Plate, ::AbstractBody, ::Val{:side_d}, _) = ( 0.0, -1.0,  0.0)

function surface_centroid(::Plate, body::AbstractBody, ::Val{:top})
    L, _, H = _plate_skin(body); (zero(L), zero(L), H/2)
end
function surface_centroid(::Plate, body::AbstractBody, ::Val{:bottom})
    L, _, H = _plate_skin(body); (zero(L), zero(L), -H/2)
end
function surface_centroid(::Plate, body::AbstractBody, ::Val{:side_a})
    L, _, _ = _plate_skin(body); (L/2, zero(L), zero(L))
end
function surface_centroid(::Plate, body::AbstractBody, ::Val{:side_b})
    L, _, _ = _plate_skin(body); (-L/2, zero(L), zero(L))
end
function surface_centroid(::Plate, body::AbstractBody, ::Val{:side_c})
    _, W, _ = _plate_skin(body); (zero(W), W/2, zero(W))
end
function surface_centroid(::Plate, body::AbstractBody, ::Val{:side_d})
    _, W, _ = _plate_skin(body); (zero(W), -W/2, zero(W))
end
surface_centroid_normal(::Plate, ::AbstractBody, ::Val{:top})    = ( 0.0,  0.0,  1.0)
surface_centroid_normal(::Plate, ::AbstractBody, ::Val{:bottom}) = ( 0.0,  0.0, -1.0)
surface_centroid_normal(::Plate, ::AbstractBody, ::Val{:side_a}) = ( 1.0,  0.0,  0.0)
surface_centroid_normal(::Plate, ::AbstractBody, ::Val{:side_b}) = (-1.0,  0.0,  0.0)
surface_centroid_normal(::Plate, ::AbstractBody, ::Val{:side_c}) = ( 0.0,  1.0,  0.0)
surface_centroid_normal(::Plate, ::AbstractBody, ::Val{:side_d}) = ( 0.0, -1.0,  0.0)
