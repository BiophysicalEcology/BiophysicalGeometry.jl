# ── Half shapes ──────────────────────────────────────────────────────────
#
# A `Half{S}` (see the type in geometry.jl) is the full shape `S` of double
# mass, cut through the centre and mirror-joined at a flat face. It inherits
# the parent's dimension math and thermal family; only surface area, the flat
# join face, and the mesh are half-specific. Two `Half`s of equal mass joined
# at their flat faces reconstruct the full shape of double mass.

"""
    HalfCylinder(mass, density, axis_ratio_b) -> Half{<:Cylinder}

Dorsal/ventral half-cylinder: `Half(Cylinder(2mass, density, axis_ratio_b))`.
Axis along `+z`, curved surface in `y ≥ 0`, flat face at `y = 0`.
"""
HalfCylinder(mass, density, axis_ratio_b) = Half(Cylinder(2mass, density, axis_ratio_b))

"""
    HalfEllipsoid(mass, density, axis_ratio_b, axis_ratio_c) -> Half{<:Ellipsoid}

Dorsal/ventral half-ellipsoid: `Half(Ellipsoid(2mass, density, b, c))`. Long
axis along `+x`, dome in `z ≥ 0`, flat elliptical face at `z = 0`.
"""
HalfEllipsoid(mass, density, axis_ratio_b, axis_ratio_c) =
    Half(Ellipsoid(2mass, density, axis_ratio_b, axis_ratio_c))

"""
    HalfSphere(mass, density) -> Half{<:Sphere}

Hemisphere: `Half(Sphere(2mass, density))`. Wrapping `Sphere` (not an equal-axis
ellipsoid) means the geometry uses the fast spherical closed forms rather than
the eccentricity/asin and fat-cubic of the ellipsoidal path.
"""
HalfSphere(mass, density) = Half(Sphere(2mass, density))

# ── Geometry: delegate dimensions to the parent, override area + volume ───

geometry(h::Half, ins::Naked)        = _halve(h, geometry(h.parent, ins))
geometry(h::Half, fat::FatLayer)     = _halve(h, geometry(h.parent, fat))
geometry(h::Half, fur::FibrousLayer) = _halve(h, geometry(h.parent, fur), fur)
geometry(h::Half, fur::FibrousLayer, fat::FatLayer) =
    _halve(h, geometry(h.parent, fur, fat), fur)

# Assemble a half Geometry from the parent's: same `length`, half the volume.
function _halfgeom(full, total)
    vol = full.volume / 2
    Geometry(vol, vol^(1 / 3), full.length, SurfaceAreas(; total))
end
function _halfgeom(full, total, skin, fur::FibrousLayer)
    vol = full.volume / 2
    convection = skin - insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    Geometry(vol, vol^(1 / 3) + fur.thickness, full.length, SurfaceAreas(; total, skin, convection))
end

# A half's surface is the parent's, halved, plus the cut face — the mirror plane
# through the centre where the two halves join. The parent area is already
# computed (`full.area.*`); only this cut face is family-specific.
_cut_face_area(::Half{<:AbstractCylindrical}, l) = 2 * l.radius_skin * l.length_skin
_cut_face_area(::Half{<:AbstractEllipsoidal}, l) = π * l.b_semi_minor_skin * l.c_semi_minor_skin
_cut_face_area(::Half{<:AbstractSpherical},   l) = π * l.radius_skin^2

function _halve(h::Half, full)
    cut_face = _cut_face_area(h, full.length)
    _halfgeom(full, full.area.total / 2 + cut_face)
end
function _halve(h::Half, full, fur::FibrousLayer)
    cut_face = _cut_face_area(h, full.length)
    _halfgeom(full, full.area.total / 2 + cut_face, full.area.skin / 2 + cut_face, fur)
end

# ── Route-around-the-wrapper forwards ────────────────────────────────────
# The half's `length` NamedTuple is the parent's, so these read identically.

skin_radius(h::Half, ins, body)       = skin_radius(h.parent, ins, body)
insulation_radius(h::Half, ins, body) = insulation_radius(h.parent, ins, body)
flesh_radius(h::Half, ins, body)      = flesh_radius(h.parent, ins, body)
outer_dims(h::Half, body::AbstractBody) = outer_dims(h.parent, body)

# The half's silhouette is exactly half the parent's, in every family.
silhouette(h::Half, ins::AbstractInsulationLayer, body::AbstractBody, θ) =
    silhouette(h.parent, ins, body, θ) / 2
function silhouette(h::Half, ins::AbstractInsulationLayer, body::AbstractBody)
    s = silhouette(h.parent, ins, body)
    (; normal = s.normal / 2, parallel = s.parallel / 2)
end

# ── Composition: surfaces are half-specific (extra flat join face) ────────
#
# Attachment positions are anchored at flesh (skin) level so joins meet
# flesh-to-flesh; the flat face is reported at skin level so mixed-insulation
# halves join cleanly under FullCover.

# Cylindrical half: EndA, EndB (semicircular), Lateral (half tube), Flat (rectangle).
attachment_surfaces(::Half{<:AbstractCylindrical}) = (EndA, EndB, Lateral, Flat)

surface_area(::Half{<:AbstractCylindrical}, body::AbstractBody, ::EndA) =
    π * insulation_radius(body)^2 / 2
surface_area(sh::Half{<:AbstractCylindrical}, body::AbstractBody, ::EndB) =
    surface_area(sh, body, EndA())
function surface_area(sh::Half{<:AbstractCylindrical}, body::AbstractBody, ::Lateral)
    d = outer_dims(sh, body)
    π * d.r * d.L
end
surface_area(::Half{<:AbstractCylindrical}, body::AbstractBody, ::Flat) =
    2 * body.geometry.length.radius_skin * body.geometry.length.length_skin

function validate_range(::Half{<:AbstractCylindrical}, body::AbstractBody, loc::EndA)
    R = body.geometry.length.radius_skin
    loc.r ≥ zero(loc.r) && loc.r ≤ R || error("EndA r out of range [0, $R]: $(loc.r)")
    0 ≤ loc.φ ≤ π || error("EndA φ out of range [0, π]: $(loc.φ)")
end
validate_range(sh::Half{<:AbstractCylindrical}, body::AbstractBody, loc::EndB) =
    validate_range(sh, body, EndA(loc.r, loc.φ))
function validate_range(::Half{<:AbstractCylindrical}, body::AbstractBody, loc::Lateral)
    L = body.geometry.length.length_skin
    loc.z ≥ zero(loc.z) && loc.z ≤ L || error("Lateral z out of range [0, $L]: $(loc.z)")
    0 ≤ loc.φ ≤ π || error("Lateral φ out of range [0, π]: $(loc.φ)")
end
function validate_range(::Half{<:AbstractCylindrical}, body::AbstractBody, loc::Flat)  # (u=z, v=x)
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    loc.u ≥ zero(loc.u) && loc.u ≤ L || error("Flat z out of range [0, $L]: $(loc.u)")
    abs(loc.v) ≤ r || error("Flat x out of range ±$r: $(loc.v)")
end

surface_point(::Half{<:AbstractCylindrical}, body::AbstractBody, loc::EndA) =
    (loc.r * cos(loc.φ), loc.r * sin(loc.φ), zero(loc.r))
surface_point(::Half{<:AbstractCylindrical}, body::AbstractBody, loc::EndB) =
    (loc.r * cos(loc.φ), loc.r * sin(loc.φ), body.geometry.length.length_skin)
function surface_point(::Half{<:AbstractCylindrical}, body::AbstractBody, loc::Lateral)
    R = body.geometry.length.radius_skin
    (R * cos(loc.φ), R * sin(loc.φ), loc.z)
end
surface_point(::Half{<:AbstractCylindrical}, body::AbstractBody, loc::Flat) =
    (loc.v, zero(loc.v), loc.u)

surface_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::EndA) = (0.0, 0.0, -1.0)
surface_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::EndB) = (0.0, 0.0, 1.0)
surface_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, loc::Lateral) =
    (cos(loc.φ), sin(loc.φ), 0.0)
surface_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::Flat) = (0.0, -1.0, 0.0)

function surface_centroid(::Half{<:AbstractCylindrical}, body::AbstractBody, ::EndA)
    R = body.geometry.length.radius_skin; (zero(R), R / 2, zero(R))
end
function surface_centroid(::Half{<:AbstractCylindrical}, body::AbstractBody, ::EndB)
    R = body.geometry.length.radius_skin; L = body.geometry.length.length_skin
    (zero(R), R / 2, L)
end
function surface_centroid(::Half{<:AbstractCylindrical}, body::AbstractBody, ::Lateral)
    R = body.geometry.length.radius_skin; L = body.geometry.length.length_skin
    (zero(R), R, L / 2)
end
function surface_centroid(::Half{<:AbstractCylindrical}, body::AbstractBody, ::Flat)
    L = body.geometry.length.length_skin; (zero(L), zero(L), L / 2)
end
surface_centroid_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::EndA) = (0.0, 0.0, -1.0)
surface_centroid_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::EndB) = (0.0, 0.0, 1.0)
surface_centroid_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::Lateral) = (0.0, 1.0, 0.0)
surface_centroid_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::Flat) = (0.0, -1.0, 0.0)

# Domed half (ellipsoidal or spherical): Dome + Flat (elliptical disc at z = 0).
# A sphere is the a=b=c case, so both families share one parametrization; only
# the skin semi-axes are read differently.
const HalfDomed = Half{<:Union{AbstractEllipsoidal,AbstractSpherical}}

_domed_semiaxes(::Half{<:AbstractEllipsoidal}, body) =
    (body.geometry.length.a_semi_major_skin, body.geometry.length.b_semi_minor_skin, body.geometry.length.c_semi_minor_skin)
function _domed_semiaxes(::Half{<:AbstractSpherical}, body)
    r = body.geometry.length.radius_skin; (r, r, r)
end

attachment_surfaces(::HalfDomed) = (Dome, Flat)

surface_area(::HalfDomed, body::AbstractBody, ::Dome) = body.geometry.area.total
function surface_area(sh::HalfDomed, body::AbstractBody, ::Flat)
    _, b, c = _domed_semiaxes(sh, body); π * b * c
end

function validate_range(::HalfDomed, ::AbstractBody, loc::Dome)
    0 ≤ loc.α ≤ π || error("Dome α out of range [0, π]: $(loc.α)")
    0 ≤ loc.β ≤ π || error("Dome β out of range [0, π]: $(loc.β)")
end
function validate_range(sh::HalfDomed, body::AbstractBody, loc::Flat)  # (u=x, v=y)
    a, b, _ = _domed_semiaxes(sh, body)
    (loc.u / a)^2 + (loc.v / b)^2 ≤ 1 + 1e-9 || error("Flat (x,y) outside boundary ellipse")
end

function surface_point(sh::HalfDomed, body::AbstractBody, loc::Dome)
    a, b, c = _domed_semiaxes(sh, body)
    (a * cos(loc.α), b * sin(loc.α) * cos(loc.β), c * sin(loc.α) * sin(loc.β))
end
surface_point(::HalfDomed, body::AbstractBody, loc::Flat) =
    (loc.u, loc.v, zero(loc.u))

function surface_normal(sh::HalfDomed, body::AbstractBody, loc::Dome)
    a, b, c = _domed_semiaxes(sh, body)
    nx = cos(loc.α) * (b * c)
    ny = sin(loc.α) * cos(loc.β) * (a * c)
    nz = sin(loc.α) * sin(loc.β) * (a * b)
    n = sqrt(nx^2 + ny^2 + nz^2)
    (nx / n, ny / n, nz / n)
end
surface_normal(::HalfDomed, ::AbstractBody, ::Flat) = (0.0, 0.0, -1.0)

function surface_centroid(sh::HalfDomed, body::AbstractBody, ::Dome)
    _, _, c = _domed_semiaxes(sh, body); (zero(c), zero(c), c)
end
function surface_centroid(sh::HalfDomed, body::AbstractBody, ::Flat)
    a, = _domed_semiaxes(sh, body); (zero(a), zero(a), zero(a))
end
surface_centroid_normal(::HalfDomed, ::AbstractBody, ::Dome) = (0.0, 0.0, 1.0)
surface_centroid_normal(::HalfDomed, ::AbstractBody, ::Flat) = (0.0, 0.0, -1.0)
