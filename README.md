[![ci](https://github.com/griswaldbrooks/marky/actions/workflows/ci.yml/badge.svg)](https://github.com/griswaldbrooks/marky/actions/workflows/ci.yml)
[![docker](https://github.com/griswaldbrooks/marky/actions/workflows/docker.yaml/badge.svg)](https://github.com/griswaldbrooks/marky/actions/workflows/docker.yaml)
[![Coverage Status](https://coveralls.io/repos/github/griswaldbrooks/marky/badge.svg?branch=)](https://coveralls.io/github/griswaldbrooks/marky?branch=)
# marky
Landmark based localization library

# Development Container
Build a new development image
```shell
mkdir -p ~/.marky/ccache
export UID=$(id -u) export GID=$(id -g); docker compose -f compose.dev.yml build
```
Start an interactive development container
```shell
docker compose -f compose.dev.yml run development
```
Build the repository in the container
```shell
cmake -S src/marky/ -B build
cmake --build build
```
Test the repository in the container
```shell
ctest --test-dir build
```
# Run
```shell
username@marky-dev:~/ws$ ./build/marky_example
```
