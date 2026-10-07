# Serve the "Build an animal" app.
#
#   julia --project=app app/run.jl
#
# PORT (default 8080), HOST (default 127.0.0.1; 0.0.0.0 in the container) and PROXY_URL (the public address, when
# served behind a proxy) are read from the environment.

using BiophysicalGeometry, Bonito, WGLMakie
import BiologicalScaling   # for legs sized by elastic or geometric similarity

server = BiophysicalGeometry.app(;
    host = get(ENV, "HOST", "127.0.0.1"),
    port = parse(Int, get(ENV, "PORT", "8080")),
    open = get(ENV, "OPEN", "true") == "true",
    proxy_url = get(ENV, "PROXY_URL", nothing),
)
wait(server)
