# ══════════════════════════════════════════════════════════════════════════════
# MULTI-PART ORGANISM TYPES
# ══════════════════════════════════════════════════════════════════════════════

"""
    BodyPart{B<:AbstractBody}

    BodyPart(body, n)

A body segment with a multiplicity count.

# Fields
- `body`: the `Body` describing shape and insulation of this segment
- `n`: number of this segment (1 for head/torso, 2 for paired limbs/wings)
"""
struct BodyPart{B<:AbstractBody}
    body::B
    n::Int
end

"""
    Join

    Join(parent, child)
    Join(parent, child, attachment)
    Join(parent, child, attachment, neck_fraction)
    Join(parent, child, attachment, neck_fraction, dihedral_angle)

Directed connection from a child body segment to a parent segment.

# Fields
- `parent`: Symbol key of the parent part in the `Organism`
- `child`: Symbol key of the child part
- `attachment`: how the child is oriented relative to the parent for plotting.
  `:top` — child above parent along spine;
  `:bottom` — child below, with small lateral offset for paired parts;
  `:lateral` — child to the side at ~70% of parent length (arms, wings);
  `:forward` — child in front along spine (head of a prone organism);
  `:back` — child behind along spine
- `neck_fraction`: for `Ellipsoid` children, the truncation plane expressed as a
  fraction of the outer equatorial radius `b_semi_minor`. 0.0 = auto (uses
  `min(r_parent, b_child)`); 0.5 = cut at half the equatorial radius ("neck").
- `dihedral_angle`: angle (radians) between the parent surface plane and the
  child surface plane at the junction. Default π/2 (perpendicular). For
  butterfly wings: angle from the dorsal body surface to the wing (PHI in
  NicheMapR's BTRFLGEOM.f).
"""
struct Join
    parent::Symbol
    child::Symbol
    attachment::Symbol
    neck_fraction::Float64
    dihedral_angle::Float64
    spine_fraction::Float64  # for :lateral attachment, fraction of parent half-length along spine (default 0.7; negative = toward tail)
end
Join(parent, child, attachment, neck_fraction, dihedral_angle) = Join(parent, child, attachment, neck_fraction, dihedral_angle, 0.7)
Join(parent, child, attachment, neck_fraction) = Join(parent, child, attachment, neck_fraction, π/2, 0.7)
Join(parent, child, attachment)               = Join(parent, child, attachment, 0.0, π/2, 0.7)
Join(parent, child)                           = Join(parent, child, :top, 0.0, π/2, 0.7)

"""
    OrganismPosture

Abstract supertype for organism posture/orientation used in visualisation.
"""
abstract type OrganismPosture end

"""
    Upright <: OrganismPosture

Spine aligned with the z-axis (vertical). Suitable for standing humans.
"""
struct Upright <: OrganismPosture end

"""
    Prone <: OrganismPosture

Spine aligned with the x-axis (horizontal). Suitable for quadrupeds, lizards,
butterflies etc.
"""
struct Prone <: OrganismPosture end

"""
    Organism{P<:NamedTuple, J<:Tuple}

    Organism(; parts, joins)

A multi-part organism composed of named [`BodyPart`](@ref)s connected by
[`Join`](@ref)s.

# Fields
- `parts`: `NamedTuple` of `BodyPart` values, e.g.
  `(head = BodyPart(..., 1), arms = BodyPart(..., 2))`
- `joins`: `Tuple` of `Join` connections

# Example
```julia
human = Organism(
    parts = (
        head  = BodyPart(Body(Ellipsoid(mass * 0.076, density, 1.6, 1.0), Naked()), 1),
        torso = BodyPart(Body(Cylinder(mass * 0.501, density, 1.9),        Naked()), 1),
        arms  = BodyPart(Body(Cylinder(mass * 0.049, density, 12.0),       Naked()), 2),
        legs  = BodyPart(Body(Cylinder(mass * 0.162, density, 7.0),        Naked()), 2),
    ),
    joins = (
        Join(:torso, :head,  :top,     0.5),
        Join(:torso, :arms,  :lateral, 0.0),
        Join(:torso, :legs,  :bottom,  0.0),
    )
)
```
"""
struct Organism{P<:NamedTuple, J<:Tuple}
    parts::P
    joins::J
end

function Organism(; parts, joins)
    _validate_organism(parts, joins)
    Organism(parts, joins)
end

function _validate_organism(parts, joins)
    part_names = Set(keys(parts))
    for j in joins
        j.parent ∈ part_names || error(
            "Join references unknown parent :$(j.parent). " *
            "Available parts: $(sort(collect(part_names)))")
        j.child ∈ part_names || error(
            "Join references unknown child :$(j.child). " *
            "Available parts: $(sort(collect(part_names)))")
        j.parent != j.child || error(
            "Join cannot connect a part to itself (:$(j.parent))")
    end
end

"""
    OrganismGeometry

Result of [`geometry`](@ref) called on an [`Organism`](@ref).

# Fields
- `total_area`: total outer surface area after all junction subtractions
- `total_volume`: sum of part volumes × multiplicity
- `per_part`: `NamedTuple` of individual [`Geometry`](@ref) structs (one per
  unique part name, not expanded by multiplicity `n`)
- `join_areas`: `NamedTuple` mapping each child part symbol to its junction
  area (the area of one circular contact — before multiplying by `n`)
- `area_fractions`: `NamedTuple` giving each part's net area as a fraction of
  `total_area`
"""
struct OrganismGeometry{A,V,P<:NamedTuple,JA<:NamedTuple,AF<:NamedTuple}
    total_area::A
    total_volume::V
    per_part::P
    join_areas::JA
    area_fractions::AF
end

# ══════════════════════════════════════════════════════════════════════════════
# JUNCTION RADIUS HELPERS
# ══════════════════════════════════════════════════════════════════════════════

"""
    junction_radius(body, neck_fraction=0.0)

Return the outer radius at the connection point for a body part.
Dispatches on shape type via `_junction_radius`.
"""
junction_radius(body::AbstractBody, neck_fraction::Float64=0.0) =
    _junction_radius(shape(body), insulation(body), body, neck_fraction)

# Cylinder: neck_fraction is irrelevant — the cap radius is always r
_junction_radius(::Cylinder, ::AbstractInsulation, body, _) = insulation_radius(body)

# Sphere: similarly, always outer radius
_junction_radius(::Sphere, ::AbstractInsulation, body, _) = insulation_radius(body)

# Ellipsoid: neck_fraction > 0 means truncate at that fraction of b_semi_minor
_junction_radius(::Ellipsoid, ::AbstractInsulation, body, neck_frac) =
    neck_frac > 0.0 ? neck_frac * insulation_radius(body) : insulation_radius(body)

function _junction_radius(sh::AbstractShape, ::AbstractInsulation, body, _)
    error(
        "junction_radius is not defined for shape type $(typeof(sh)). " *
        "Define `_junction_radius(::$(typeof(sh)), ins, body, neck_frac)` or " *
        "implement `junction_area` directly for this shape pair.")
end

# ══════════════════════════════════════════════════════════════════════════════
# JUNCTION AREA  (multiple dispatch on parent/child shape types)
# ══════════════════════════════════════════════════════════════════════════════

"""
    junction_area(parent, child, join) -> area

Area of the circular (or face) contact zone where `child` attaches to `parent`.
Both parts lose this area in surface-area accounting.

Dispatch table:
- `Cylinder` child → `π × r_child_outer²`  (end-cap defines junction)
- `Cylinder`–`Cylinder` → `π × min(r_parent, r_child)²`
- `Plate` child → edge face area depending on `join.attachment`
- Generic (Ellipsoid, …) child → `π × r²` where `r` uses `join.neck_fraction`
"""

# Cylinder child attaching to ANY parent — cap radius
function junction_area(parent::AbstractBody,
                       child::Body{<:Cylinder},
                       join::Join)
    π * insulation_radius(child)^2
end

# Cylinder–Cylinder — use minimum radius (larger cylinder "hosts" the smaller)
function junction_area(parent::Body{<:Cylinder},
                       child::Body{<:Cylinder},
                       join::Join)
    π * min(insulation_radius(parent), insulation_radius(child))^2
end

# Plate child (wing, fin, ear, etc.) — face area depends on attachment direction
function junction_area(parent::AbstractBody,
                       child::Body{<:Plate},
                       join::Join)
    gl = child.geometry.length
    if join.attachment ∈ (:top, :bottom)
        # Flat face (length × width) is the junction
        gl.length_skin * gl.width_skin
    else
        # Lateral/inboard edge (width × height) for wings, fins, and ears
        gl.width_skin * gl.height_skin
    end
end

# Generic fallback — Ellipsoid child or any unspecialised shape
function junction_area(parent::AbstractBody, child::AbstractBody, join::Join)
    r_child  = junction_radius(child, join.neck_fraction)
    r_parent = insulation_radius(parent)
    r = join.neck_fraction > 0.0 ? r_child : min(r_parent, r_child)
    π * r^2
end

# ══════════════════════════════════════════════════════════════════════════════
# GEOMETRY COMPUTATION
# ══════════════════════════════════════════════════════════════════════════════

"""
    geometry(org::Organism) -> OrganismGeometry

Compute surface areas, volumes, and junction corrections for a multi-part
organism.

For each `Join(parent, child)` with child multiplicity `n`:
- Both `parent` and `child` lose `junction_area × n` from their gross surface area.
"""
function geometry(org::Organism)
    parts = org.parts

    # Step 1: compute junction area for each join
    junction_triples = map(org.joins) do j
        parent_bp = parts[j.parent]
        child_bp  = parts[j.child]
        ja   = junction_area(parent_bp.body, child_bp.body, j)
        ja_n = ja * child_bp.n   # total area to remove from each connected part
        (parent=j.parent, child=j.child, area=ja, area_n=ja_n)
    end

    # Step 2: accumulate per-part deductions
    deduction_dict = Dict{Symbol, Any}()
    for t in junction_triples
        jn = t.area_n
        deduction_dict[t.parent] = get(deduction_dict, t.parent, zero(jn)) + jn
        deduction_dict[t.child]  = get(deduction_dict, t.child,  zero(jn)) + jn
    end

    # Step 3: net areas as NamedTuple
    part_names  = keys(parts)
    net_areas   = NamedTuple{part_names}(
        Tuple(
            parts[n].body.geometry.area.total * parts[n].n -
            get(deduction_dict, n, zero(parts[n].body.geometry.area.total))
            for n in part_names
        )
    )

    total_area_val = sum(values(net_areas))
    total_vol      = sum(bp -> bp.body.geometry.volume * bp.n, values(parts))
    per_part_geom  = map(bp -> bp.body.geometry, parts)
    # join_areas keyed by child symbol; last value wins for repeated children
    join_areas_nt  = NamedTuple(j.child => t.area for (j, t) in zip(org.joins, junction_triples))
    area_fracs     = map(a -> a / total_area_val, net_areas)

    return OrganismGeometry(total_area_val, total_vol, per_part_geom, join_areas_nt, area_fracs)
end

# ══════════════════════════════════════════════════════════════════════════════
# ACCESSOR FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

"""Total outer surface area of the organism after junction corrections."""
total_area(org::Organism) = geometry(org).total_area

"""Total volume of all parts (multiplicity-weighted)."""
total_volume(org::Organism) = geometry(org).total_volume

"""Total skin surface area after junction corrections."""
function skin_area(org::Organism)
    parts = org.parts
    junction_triples = map(org.joins) do j
        parent_bp = parts[j.parent]
        child_bp  = parts[j.child]
        ja   = junction_area(parent_bp.body, child_bp.body, j)
        ja_n = ja * child_bp.n
        (parent=j.parent, child=j.child, area_n=ja_n)
    end
    deduction_dict = Dict{Symbol, Any}()
    for t in junction_triples
        deduction_dict[t.parent] = get(deduction_dict, t.parent, zero(t.area_n)) + t.area_n
        deduction_dict[t.child]  = get(deduction_dict, t.child,  zero(t.area_n)) + t.area_n
    end
    sum(
        parts[n].body.geometry.area.skin * parts[n].n -
        get(deduction_dict, n, zero(parts[n].body.geometry.area.skin))
        for n in keys(parts)
    )
end

# ══════════════════════════════════════════════════════════════════════════════
# PART POSITIONS  (used by Makie extension for 3-D visualisation)
# ══════════════════════════════════════════════════════════════════════════════

"""
    PartPosition

3-D position and axis direction of a body part in world coordinates.

# Fields
- `center`: `(x, y, z)` coordinates of the part centre in metres
- `axis`: unit vector along the part's long axis
"""
struct PartPosition
    center::NTuple{3, Float64}
    axis::NTuple{3, Float64}
end

# ── helpers: half-size along long axis (metres, unit-stripped) ──────────────
_half_axis_m(body::Body{Cylinder})  = ustrip(u"m", body.geometry.length.length_skin) / 2.0
_half_axis_m(body::Body{Ellipsoid}) = ustrip(u"m", body.geometry.length.a_semi_major_skin)
_half_axis_m(body::Body{Sphere})    = ustrip(u"m", insulation_radius(body))
_half_axis_m(body::Body{Plate})     = ustrip(u"m", body.geometry.length.length_skin) / 2.0
_half_axis_m(body::AbstractBody)    = ustrip(u"m", insulation_radius(body))  # fallback

# ── helpers: transverse radius (metres, unit-stripped) ──────────────────────
_transverse_r_m(body::Body{Cylinder})  = ustrip(u"m", insulation_radius(body))
_transverse_r_m(body::Body{Ellipsoid}) = ustrip(u"m", insulation_radius(body))
_transverse_r_m(body::Body{Sphere})    = ustrip(u"m", insulation_radius(body))
_transverse_r_m(body::Body{Plate})     = ustrip(u"m", body.geometry.length.width_skin) / 2.0
_transverse_r_m(body::AbstractBody)    = ustrip(u"m", insulation_radius(body))  # fallback

# ── posture direction vectors ────────────────────────────────────────────────
_spine_dir(::Upright) = (0.0, 0.0,  1.0)   # z up
_spine_dir(::Prone)   = (1.0, 0.0,  0.0)   # x forward/spine
_side_dir(::Upright)  = (1.0, 0.0,  0.0)   # x lateral
_side_dir(::Prone)    = (0.0, 1.0,  0.0)   # y lateral
_fwd_dir(::Upright)   = (0.0, 1.0,  0.0)   # y forward
_fwd_dir(::Prone)     = (1.0, 0.0,  0.0)   # x forward (same as spine)
_down_dir(::Upright)  = (0.0, 0.0, -1.0)
_down_dir(::Prone)    = (0.0, 0.0, -1.0)   # gravity always -z

"""
    part_positions(org::Organism, posture::OrganismPosture) -> NamedTuple

Compute the 3-D centre position and axis direction for each body part for
rendering by the Makie extension. Positions are in metres.

The root part (the one that is not a child in any `Join`) is placed at the
origin. All other parts are placed relative to their parent using the
`attachment` field of the `Join` and the organism's `posture`.
"""
function part_positions(org::Organism, posture::OrganismPosture)
    parts = org.parts
    joins = org.joins

    # Identify root (not a child in any join)
    child_set = Set(j.child for j in joins)
    root_name = first(k for k in keys(parts) if k ∉ child_set)

    spine = _spine_dir(posture)
    side  = _side_dir(posture)
    down  = _down_dir(posture)

    pos_dict = Dict{Symbol, PartPosition}()
    pos_dict[root_name] = PartPosition((0.0, 0.0, 0.0), spine)

    for j in joins
        parent_pos  = pos_dict[j.parent]
        parent_body = parts[j.parent].body
        child_body  = parts[j.child].body

        pc           = parent_pos.center
        parent_half  = _half_axis_m(parent_body)
        child_half_full = _half_axis_m(child_body)
        # For Ellipsoid children with neck truncation: the cut plane is at
        # a * sqrt(1 - nc²) below the ellipsoid centre, so the effective
        # half-height (centre → cut plane) replaces the full a_semi_major.
        child_half = if child_body.shape isa Ellipsoid && j.neck_fraction > 0.0
            nc = j.neck_fraction
            child_half_full * sqrt(1.0 - nc^2)
        else
            child_half_full
        end
        parent_r     = _transverse_r_m(parent_body)
        child_r      = _transverse_r_m(child_body)
        lat_spacing  = parent_r + child_r

        center, ax = if j.attachment == :top
            c = (pc[1] + spine[1]*(parent_half + child_half),
                 pc[2] + spine[2]*(parent_half + child_half),
                 pc[3] + spine[3]*(parent_half + child_half))
            (c, spine)

        elseif j.attachment == :bottom
            # Below spine and slightly to one side (for n=2 pairs)
            c = (pc[1] + side[1]*child_r + down[1]*(parent_half + child_half),
                 pc[2] + side[2]*child_r + down[2]*(parent_half + child_half),
                 pc[3] + side[3]*child_r + down[3]*(parent_half + child_half))
            (c, spine)

        elseif j.attachment == :lateral
            # To side, at spine_fraction of parent half-length (positive = toward head/front)
            attach_along_spine = parent_half * j.spine_fraction
            c = (pc[1] + side[1]*lat_spacing  + spine[1]*attach_along_spine
                       - spine[1]*child_half,
                 pc[2] + side[2]*lat_spacing  + spine[2]*attach_along_spine
                       - spine[2]*child_half,
                 pc[3] + side[3]*lat_spacing  + spine[3]*attach_along_spine
                       - spine[3]*child_half)
            # For Prone, limbs (Cylinder) hang downward; Plate wings/fins stay level
            c_down = if posture isa Prone && !(child_body.shape isa Plate)
                (c[1], c[2],
                 c[3] + down[3]*(parent_r + child_r))
            else
                c
            end
            (c_down, spine)

        elseif j.attachment == :forward
            fwd = _fwd_dir(posture)
            c = (pc[1] + fwd[1]*(parent_half + child_half),
                 pc[2] + fwd[2]*(parent_half + child_half),
                 pc[3] + fwd[3]*(parent_half + child_half))
            (c, fwd)

        elseif j.attachment == :back
            fwd = _fwd_dir(posture)
            c = (pc[1] - fwd[1]*(parent_half + child_half),
                 pc[2] - fwd[2]*(parent_half + child_half),
                 pc[3] - fwd[3]*(parent_half + child_half))
            (c, (.-fwd...,))

        elseif j.attachment == :ear
            # Attach near top of parent, to the side; child extends further along spine.
            # Suitable for ears (Plate), which sit on top of the head and point upward.
            lat_off = parent_r + child_r
            c = (pc[1] + side[1]*lat_off   + spine[1]*(parent_half + child_half),
                 pc[2] + side[2]*lat_off   + spine[2]*(parent_half + child_half),
                 pc[3] + side[3]*lat_off   + spine[3]*(parent_half + child_half))
            (c, spine)

        else  # default: top
            c = (pc[1] + spine[1]*(parent_half + child_half),
                 pc[2] + spine[2]*(parent_half + child_half),
                 pc[3] + spine[3]*(parent_half + child_half))
            (c, spine)
        end

        pos_dict[j.child] = PartPosition(center, ax)
    end

    return NamedTuple{keys(parts)}(Tuple(pos_dict[n] for n in keys(parts)))
end

# ══════════════════════════════════════════════════════════════════════════════
# BUTTERFLY CONFIGURATION FACTORS STUB
# ══════════════════════════════════════════════════════════════════════════════

"""
    butterfly_config_factors(org::Organism, wing_angle) -> NamedTuple

Compute radiation view/configuration factors between body and wing surfaces
for a butterfly-like organism.

## TODO
Port algorithm from `NicheMapR/src/BTRFLGEOM.f`. Key variables:
- `SIDE2` = body dorsal width; `SIDE1 = SIDE3` = wing width; `BODLEN` = body length
- `PHI` (`wing_angle`): angle from dorsal surface to wing
  - `PHI = π/2`: wings horizontal; `PHI > π/2`: V-open; `PHI < π/2`: tent/closed
- Returns: `(F12, F21, F13, F24, F25, F26, …)` between wing (1,3), body back (2),
  and imaginary closure surfaces (4,5,6)
"""
function butterfly_config_factors(org::Organism, wing_angle::Real)
    error(
        "butterfly_config_factors is not yet implemented. " *
        "See NicheMapR/src/BTRFLGEOM.f for the reference Fortran algorithm.")
end
