# Build an animal

The interactive "Build an animal" app. The controls run in the browser; the animal is built in Julia by
BiophysicalGeometry.jl, from the code the app shows.

Run it locally:

```sh
julia --project=app -e 'using Pkg; Pkg.instantiate()'
julia --project=app app/run.jl
```

Or in a container, from the repository root:

```sh
docker build -f app/Dockerfile -t animal-builder .
docker run -p 8080:8080 animal-builder
```

| Variable | Default | Meaning |
|:--|:--|:--|
| `PORT` | `8080` | port to listen on |
| `HOST` | `127.0.0.1` (`0.0.0.0` in the container) | interface to listen on |
| `OPEN` | `true` (`false` in the container) | open a browser on start |
| `PROXY_URL` | none | the public address, when served behind a proxy, e.g. `https://example.org/builder/` |

To show the app in the documentation, build the docs with `BUILDER_URL` set to its public address.
