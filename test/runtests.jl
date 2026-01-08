using BiophysicalGeometry
using Aqua
using SafeTestsets
using Test

@testset "Quality assurance" begin
    Aqua.test_unbound_args(BiophysicalGeometry)
    Aqua.test_stale_deps(BiophysicalGeometry)
    Aqua.test_undefined_exports(BiophysicalGeometry)
    Aqua.test_project_extras(BiophysicalGeometry)
    Aqua.test_deps_compat(BiophysicalGeometry)
end

@safetestset "geometry" begin include("geometry.jl") end
