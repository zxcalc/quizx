#!/bin/bash

cargo build --release
cp target/release/libquizx.so .
