const LABEL_WIDTH = 22  # field-label column width

# ── Marker types ──────────────────────────────────────────────────────────────

Base.show(io::IO, ::MIME"text/plain", ::Naked)              = print(io, "Naked()")
Base.show(io::IO, ::MIME"text/plain", ::NormalToSun)        = print(io, "NormalToSun()")
Base.show(io::IO, ::MIME"text/plain", ::ParallelToSun)      = print(io, "ParallelToSun()")
Base.show(io::IO, ::MIME"text/plain", ::Intermediate)       = print(io, "Intermediate()")
Base.show(io::IO, ::MIME"text/plain", ::ZenithAngleVarying) = print(io, "ZenithAngleVarying()")

# ── Shape types ───────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", s::Sphere)
    println(io, "Sphere")
    println(io, rpad("  mass:", LABEL_WIDTH), s.mass)
    print(io,   rpad("  density:", LABEL_WIDTH), s.density)
end

function Base.show(io::IO, ::MIME"text/plain", s::Cylinder)
    println(io, "Cylinder")
    println(io, rpad("  mass:", LABEL_WIDTH), s.mass)
    println(io, rpad("  density:", LABEL_WIDTH), s.density)
    print(io,   rpad("  axis_ratio_b:", LABEL_WIDTH), s.axis_ratio_b)
end

function Base.show(io::IO, ::MIME"text/plain", s::Ellipsoid)
    println(io, "Ellipsoid")
    println(io, rpad("  mass:", LABEL_WIDTH), s.mass)
    println(io, rpad("  density:", LABEL_WIDTH), s.density)
    println(io, rpad("  axis_ratio_b:", LABEL_WIDTH), s.axis_ratio_b)
    print(io,   rpad("  axis_ratio_c:", LABEL_WIDTH), s.axis_ratio_c)
    iszero(s.pole_a_truncation) || print(io, "\n", rpad("  pole_a_truncation:", LABEL_WIDTH), s.pole_a_truncation)
end

function Base.show(io::IO, ::MIME"text/plain", s::Plate)
    println(io, "Plate")
    println(io, rpad("  mass:", LABEL_WIDTH), s.mass)
    println(io, rpad("  density:", LABEL_WIDTH), s.density)
    println(io, rpad("  axis_ratio_b:", LABEL_WIDTH), s.axis_ratio_b)
    print(io,   rpad("  axis_ratio_c:", LABEL_WIDTH), s.axis_ratio_c)
end

function Base.show(io::IO, ::MIME"text/plain", s::TriangularPlate)
    println(io, "TriangularPlate")
    println(io, rpad("  mass:", LABEL_WIDTH), s.mass)
    println(io, rpad("  density:", LABEL_WIDTH), s.density)
    println(io, rpad("  axis_ratio_b:", LABEL_WIDTH), s.axis_ratio_b)
    print(io,   rpad("  axis_ratio_c:", LABEL_WIDTH), s.axis_ratio_c)
end

function Base.show(io::IO, ::MIME"text/plain", s::Cone)
    println(io, "Cone")
    println(io, rpad("  mass:", LABEL_WIDTH), s.mass)
    println(io, rpad("  density:", LABEL_WIDTH), s.density)
    println(io, rpad("  axis_ratio_b:", LABEL_WIDTH), s.axis_ratio_b)
    print(io,   rpad("  top_ratio:", LABEL_WIDTH), s.top_ratio)
end

function Base.show(io::IO, mime::MIME"text/plain", h::Half)
    # The parent is the full shape of double mass; show the half's own mass first.
    println(io, "Half (mass $(mass(h))) of")
    show(io, mime, h.parent)
end

# ── Insulation types ──────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", f::FibrousLayer)
    println(io, "FibrousLayer")
    println(io, rpad("  thickness:", LABEL_WIDTH), f.thickness)
    println(io, rpad("  fibre_diameter:", LABEL_WIDTH), f.fibre_diameter)
    print(io,   rpad("  fibre_density:", LABEL_WIDTH), f.fibre_density)
end

function Base.show(io::IO, ::MIME"text/plain", f::FatLayer)
    println(io, "FatLayer")
    println(io, rpad("  fraction:", LABEL_WIDTH), f.fraction)
    print(io,   rpad("  density:", LABEL_WIDTH), f.density)
end

function Base.show(io::IO, ::MIME"text/plain", c::CompositeInsulation)
    n = length(c.layers)
    println(io, "CompositeInsulation ($n layer$(n == 1 ? "" : "s"))")
    for (i, layer) in enumerate(c.layers)
        print(io, "  [$i] ")
        show(io, MIME"text/plain"(), layer)
        i < n && println(io)
    end
end

# ── Internal helper ───────────────────────────────────────────────────────────

function _show_area_lines(io, a::SurfaceAreas, indent, w)
    items = Pair{String,Any}["total" => a.total]
    a.skin !== a.total       && push!(items, "skin"       => a.skin)
    a.convection !== a.total && push!(items, "convection" => a.convection)
    !isnothing(a.ventral)    && push!(items, "ventral"    => a.ventral)
    for (i, (k, v)) in enumerate(items)
        i < length(items) ? println(io, rpad("$indent$k:", w), v) : print(io, rpad("$indent$k:", w), v)
    end
end

# ── SurfaceAreas ──────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", a::SurfaceAreas)
    println(io, "Surface Areas")
    _show_area_lines(io, a, "  ", LABEL_WIDTH)
end

# ── Geometry ──────────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", g::Geometry)
    println(io, "Geometry")
    println(io, rpad("  volume:", LABEL_WIDTH), g.volume)
    println(io, "  Dimensions:")
    for (k, v) in pairs(g.length)
        println(io, rpad("    $k:", LABEL_WIDTH + 2), v)
    end
    println(io, "  Surface Areas:")
    _show_area_lines(io, g.area, "    ", LABEL_WIDTH + 2)
end

# ── Body ──────────────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", b::Body)
    header = "Body{$(nameof(typeof(b.shape))), $(nameof(typeof(b.insulation)))}"
    println(io, header)
    println(io, "─" ^ length(header))

    print(io, "Shape:      ")
    show(io, MIME"text/plain"(), b.shape)
    println(io)
    println(io)

    print(io, "Insulation: ")
    show(io, MIME"text/plain"(), b.insulation)
    println(io)
    println(io)

    println(io, "Geometry:")
    println(io, rpad("  volume:", LABEL_WIDTH), b.geometry.volume)
    println(io, "  Dimensions:")
    for (k, v) in pairs(b.geometry.length)
        println(io, rpad("    $k:", LABEL_WIDTH + 2), v)
    end
    println(io, "  Surface Areas:")
    _show_area_lines(io, b.geometry.area, "    ", LABEL_WIDTH + 2)
end
