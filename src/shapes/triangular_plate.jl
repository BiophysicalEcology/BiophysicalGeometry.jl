"""
    TriangularPlate(; mass, density, volume, length, width, height, axis_ratio_b, axis_ratio_c) <: AbstractShape

A plate cut in half along its diagonal: a right-triangular prism, for wings and
ears. The right angle sits at the origin, with the two legs `length` along `x`
and `width` along `y`, the long `Diagonal` edge joining their ends, and the
thickness `height` along `z`, centred on `z = 0`. `axis_ratio_b` is length /
width and `axis_ratio_c` length / height, as for [`Plate`](@ref). Give any
sufficient set of keywords and the rest is solved for. Dimensions are at skin
level; a fibrous layer offsets every face outward by its thickness.
"""
struct TriangularPlate{M,D,B,C} <: AbstractSlab
    mass::M
    density::D
    axis_ratio_b::B
    axis_ratio_c::C
    TriangularPlate(::Resolved, mass::M, density::D, axis_ratio_b::B, axis_ratio_c::C) where {M,D,B,C} =
        new{M,D,B,C}(mass, density, axis_ratio_b, axis_ratio_c)
end

# volume = length·width·height / 2; ratios as for the box.
const TRIANGLE_SPEC = ShapeSpec((:length, :width, :height), (1, 1, 1), log(1 / 2),
                                  (:axis_ratio_b => (1, 2, 1.0), :axis_ratio_c => (1, 3, 1.0)))

function TriangularPlate(; kw...)
    s = _resolve_shape("TriangularPlate", TRIANGLE_SPEC, NamedTuple(kw))
    TriangularPlate(RESOLVED, s.mass, s.density, s.axis_ratio_b, s.axis_ratio_c)
end

_diagonal(length, width) = sqrt(length^2 + width^2)
# Inradius of the right triangle (legs L, W): area / semi-perimeter.
inradius(length, width) = length * width / (length + width + _diagonal(length, width))

surface_area(::TriangularPlate, length, width, height) =
    length * width + (length + width + _diagonal(length, width)) * height

function _skin_level(shape::TriangularPlate, volume)
    length_skin = cbrt(2 * volume * shape.axis_ratio_b * shape.axis_ratio_c)
    width_skin = length_skin / shape.axis_ratio_b
    height_skin = length_skin / shape.axis_ratio_c
    (; dims = (; length_skin, width_skin, height_skin),
       area = surface_area(shape, length_skin, width_skin, height_skin))
end
# Offsetting every edge of a triangle outward by t gives the similar triangle
# scaled about its incentre by (r + t)/r, r the inradius; the thickness grows by
# t on each face. The right-angle corner moves to (-t, -t).
function _fibrous_level(shape::TriangularPlate, skin, thickness)
    scale = 1 + thickness / inradius(skin.length_skin, skin.width_skin)
    length_fibrous = skin.length_skin * scale
    width_fibrous = skin.width_skin * scale
    height_fibrous = skin.height_skin + 2 * thickness
    (; dims = (; length_fibrous, width_fibrous, height_fibrous),
       area = surface_area(shape, length_fibrous, width_fibrous, height_fibrous))
end
# As for `Plate`, the flesh is the same prism scaled to the flesh volume, and the
# fat thickness is the difference in equivalent (in)radius.
function _fat_thickness(shape::TriangularPlate, skin, flesh_volume, fat_volume)
    flesh = _skin_level(shape, flesh_volume).dims
    inradius(skin.length_skin, skin.width_skin) - inradius(flesh.length_skin, flesh.width_skin)
end

# Silhouette. A convex polytope's shadow along a unit direction d is half the sum
# of |n·d|·area over its faces. For the prism that is
#     |d_z|·L·W/2 + H·(W·|d_x| + L·|d_y| + |W·d_x + L·d_y|)/2,
# the last term the diagonal face (normal (W, L, 0)/√(L² + W²)). As for the other
# shapes the sun moves in the x–z plane, d = (cos θ, 0, sin θ): θ = π/2 looks down
# on the triangle (`normal`), θ = 0 along the length (`parallel`).
function _prism_silhouette(length, width, height, d)
    abs(d[3]) * length * width / 2 +
        height * (width * abs(d[1]) + length * abs(d[2]) + abs(width * d[1] + length * d[2])) / 2
end
function silhouette(sh::TriangularPlate, ::AbstractInsulationLayer, body::AbstractBody, θ)
    o = outer_dims(sh, body)
    _prism_silhouette(o.length, o.width, o.height, (cos(θ), 0.0, sin(θ)))
end
function silhouette(sh::TriangularPlate, ins::AbstractInsulationLayer, body::AbstractBody)
    (; normal = silhouette(sh, ins, body, π / 2), parallel = silhouette(sh, ins, body, 0.0))
end

# Radius accessors — the equivalent radius is the triangle's inradius.
_skin_radius(::TriangularPlate, l) = inradius(l.length_skin, l.width_skin)
_fibrous_radius(::TriangularPlate, l) = inradius(l.length_fibrous, l.width_fibrous)

# Composition: the two triangular faces, the two leg faces and the diagonal.

attachment_surfaces(::TriangularPlate) = (Top, Bottom, SideB, SideD, Diagonal)

# Outer (insulation-aware) dimensions, as for `Plate`.
outer_dims(sh::TriangularPlate, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::TriangularPlate, ::Union{Naked,FatLayer}, body::AbstractBody) =
    (length = body.geometry.length.length_skin,
     width = body.geometry.length.width_skin,
     height = body.geometry.length.height_skin)
outer_dims(::TriangularPlate, ::FibrousLayer, body::AbstractBody) =
    (length = body.geometry.length.length_fibrous,
     width = body.geometry.length.width_fibrous,
     height = body.geometry.length.height_fibrous)

function _triangle_skin(body::AbstractBody)
    gl = body.geometry.length
    (gl.length_skin, gl.width_skin, gl.height_skin)
end

# Surface areas are outer (insulation-aware), as elsewhere.
function surface_area(sh::TriangularPlate, body::AbstractBody, ::Union{Top,Bottom})
    o = outer_dims(sh, body); o.length * o.width / 2
end
function surface_area(sh::TriangularPlate, body::AbstractBody, ::SideB)
    o = outer_dims(sh, body); o.width * o.height
end
function surface_area(sh::TriangularPlate, body::AbstractBody, ::SideD)
    o = outer_dims(sh, body); o.length * o.height
end
function surface_area(sh::TriangularPlate, body::AbstractBody, ::Diagonal)
    o = outer_dims(sh, body); _diagonal(o.length, o.width) * o.height
end

# Attachment positions are at skin level.
function validate_range(::TriangularPlate, body::AbstractBody, loc::Union{Top,Bottom})
    L, W, _ = _triangle_skin(body)
    loc.x ≥ zero(loc.x) && loc.y ≥ zero(loc.y) && loc.x / L + loc.y / W ≤ 1 + 1e-9 ||
        error("$(nameof(typeof(loc))) ($(loc.x), $(loc.y)) is outside the triangle")
end
function validate_range(::TriangularPlate, body::AbstractBody, loc::SideB)
    _, W, H = _triangle_skin(body)
    zero(W) ≤ loc.y ≤ W || error("SideB y out of range [0, $W]: $(loc.y)")
    abs(loc.z) ≤ H / 2 || error("SideB z out of range ±$(H/2): $(loc.z)")
end
function validate_range(::TriangularPlate, body::AbstractBody, loc::SideD)
    L, _, H = _triangle_skin(body)
    zero(L) ≤ loc.x ≤ L || error("SideD x out of range [0, $L]: $(loc.x)")
    abs(loc.z) ≤ H / 2 || error("SideD z out of range ±$(H/2): $(loc.z)")
end
function validate_range(::TriangularPlate, body::AbstractBody, loc::Diagonal)
    L, W, H = _triangle_skin(body)
    D = _diagonal(L, W)
    zero(D) ≤ loc.position ≤ D || error("Diagonal position out of range [0, $D]: $(loc.position)")
    abs(loc.z) ≤ H / 2 || error("Diagonal z out of range ±$(H/2): $(loc.z)")
end

function surface_point(::TriangularPlate, body::AbstractBody, loc::Top)
    _, _, H = _triangle_skin(body); (loc.x, loc.y, H / 2)
end
function surface_point(::TriangularPlate, body::AbstractBody, loc::Bottom)
    _, _, H = _triangle_skin(body); (loc.x, loc.y, -H / 2)
end
surface_point(::TriangularPlate, body::AbstractBody, loc::SideB) = (zero(loc.y), loc.y, loc.z)
surface_point(::TriangularPlate, body::AbstractBody, loc::SideD) = (loc.x, zero(loc.x), loc.z)
function surface_point(::TriangularPlate, body::AbstractBody, loc::Diagonal)
    L, W, _ = _triangle_skin(body)
    f = loc.position / _diagonal(L, W)
    (L * (1 - f), W * f, loc.z)
end

function _diagonal_normal(body)
    L, W, _ = _triangle_skin(body)
    D = _diagonal(L, W)
    (W / D, L / D, 0.0)
end
surface_normal(::TriangularPlate, ::AbstractBody, ::Top) = (0.0, 0.0, 1.0)
surface_normal(::TriangularPlate, ::AbstractBody, ::Bottom) = (0.0, 0.0, -1.0)
surface_normal(::TriangularPlate, ::AbstractBody, ::SideB) = (-1.0, 0.0, 0.0)
surface_normal(::TriangularPlate, ::AbstractBody, ::SideD) = (0.0, -1.0, 0.0)
surface_normal(::TriangularPlate, body::AbstractBody, ::Diagonal) = _diagonal_normal(body)

function surface_centroid(::TriangularPlate, body::AbstractBody, ::Top)
    L, W, H = _triangle_skin(body); (L / 3, W / 3, H / 2)
end
function surface_centroid(::TriangularPlate, body::AbstractBody, ::Bottom)
    L, W, H = _triangle_skin(body); (L / 3, W / 3, -H / 2)
end
function surface_centroid(::TriangularPlate, body::AbstractBody, ::SideB)
    _, W, _ = _triangle_skin(body); (zero(W), W / 2, zero(W))
end
function surface_centroid(::TriangularPlate, body::AbstractBody, ::SideD)
    L, _, _ = _triangle_skin(body); (L / 2, zero(L), zero(L))
end
function surface_centroid(::TriangularPlate, body::AbstractBody, ::Diagonal)
    L, W, _ = _triangle_skin(body); (L / 2, W / 2, zero(L))
end
surface_centroid_normal(sh::TriangularPlate, body::AbstractBody, loc::Union{Top,Bottom,SideB,SideD}) =
    surface_normal(sh, body, loc)
surface_centroid_normal(::TriangularPlate, body::AbstractBody, ::Diagonal) = _diagonal_normal(body)
