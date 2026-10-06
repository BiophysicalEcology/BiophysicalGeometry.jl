"""
    Plate <: AbstractShape

A flat plate-shaped organism shape.
"""
struct Plate{M,D,B,C} <: AbstractSlab
    mass::M
    density::D
    axis_ratio_b::B
    axis_ratio_c::C
end

function _skin_level(shape::Plate, volume)
    length_skin = cbrt(volume * shape.axis_ratio_b * shape.axis_ratio_c)
    width_skin = length_skin / shape.axis_ratio_b
    height_skin = length_skin / shape.axis_ratio_c
    (; dims = (; length_skin, width_skin, height_skin),
       area = surface_area(shape, length_skin, width_skin, height_skin))
end
function _fibrous_level(shape::Plate, skin, thickness)
    length_fibrous = skin.length_skin + thickness * 2
    width_fibrous = skin.width_skin + thickness * 2
    height_fibrous = skin.height_skin + thickness * 2
    (; dims = (; length_fibrous, width_fibrous, height_fibrous),
       area = surface_area(shape, length_fibrous, width_fibrous, height_fibrous))
end
function _fat_thickness(shape::Plate, skin, flesh_volume, fat_volume)
    width_flesh = cbrt(flesh_volume * shape.axis_ratio_b * shape.axis_ratio_c) / shape.axis_ratio_b
    (skin.width_skin - width_flesh) / 2
end

# Surface area

surface_area(shape::Plate, length, width, height) =
    length * width * 2 + length * height * 2 + width * height * 2

# Silhouette area

# A box's shadow along a unit direction d is |d_x|·A_yz + |d_y|·A_xz + |d_z|·A_xy.
# As for the other shapes θ runs from the long axis: the sun moves in the plane
# of the smallest face's normal (θ = 0, `parallel`) and the largest face's
# normal (θ = π/2, `normal`), so only those two faces cast shadow.
function _plate_faces(ins, body)
    (length, width, height) = _plate_outer_dims(ins, body)
    sides = (length * width, length * height, height * width)
    (min(sides...), max(sides...))
end

function silhouette(shape::Plate, ins::AbstractInsulationLayer, body::AbstractBody, θ)
    (smallest, largest) = _plate_faces(ins, body)
    abs(cos(θ)) * smallest + abs(sin(θ)) * largest
end
function silhouette(shape::Plate, ins::AbstractInsulationLayer, body::AbstractBody)
    (smallest, largest) = _plate_faces(ins, body)
    return (; normal=largest, parallel=smallest)
end

_plate_outer_dims(::Union{Naked,FatLayer}, body) = _plate_skin_outer_lengths(body.geometry.length)
_plate_outer_dims(::Union{FibrousLayer,CompositeInsulation}, body) = _plate_fibrous_outer_lengths(body.geometry.length)

_plate_skin_outer_lengths(length) = (length.length_skin, length.width_skin, length.height_skin)
_plate_fibrous_outer_lengths(length) = (length.length_fibrous, length.width_fibrous, length.height_fibrous)

# Characteristic dimension

shortest_outer_dim(::Plate, ::Union{Naked,FatLayer}, geom) =
    min(_plate_skin_outer_lengths(geom.length)...)
shortest_outer_dim(::Plate, ::Union{FibrousLayer,CompositeInsulation}, geom) =
    min(_plate_fibrous_outer_lengths(geom.length)...)

# Radius accessors — Plate uses width/2 as the equivalent radius

_skin_radius(::Plate, length) = length.width_skin / 2
_fibrous_radius(::Plate, length) = length.width_fibrous / 2

# Composition

attachment_surfaces(::Plate) = (Top, Bottom, SideA, SideB, SideC, SideD)

# Outer (insulation-aware) dimensions (L, W, H). Insulation-dispatched.
outer_dims(sh::Plate, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::Plate, ::Union{Naked,FatLayer}, body::AbstractBody) =
    (L = body.geometry.length.length_skin,
     W = body.geometry.length.width_skin,
     H = body.geometry.length.height_skin)
outer_dims(::Plate, ::FibrousLayer, body::AbstractBody) =
    (L = body.geometry.length.length_fibrous,
     W = body.geometry.length.width_fibrous,
     H = body.geometry.length.height_fibrous)

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
