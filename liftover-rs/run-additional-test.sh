#!/bin/bash

set -xeu -o pipefail

cargo test --release -F test_all -- "${@}"
