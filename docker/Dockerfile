FROM rust:1 AS build
RUN apt-get update && apt-get install git musl-dev musl-tools
RUN rustup update stable && rustup target add x86_64-unknown-linux-musl
WORKDIR /project
RUN git clone https://github.com/informationsea/transanno.git
WORKDIR /project/transanno
RUN git checkout ${SOURCE_COMMIT}
RUN ./prepare-test-files.sh && cargo test --release
RUN cargo build --release --target=x86_64-unknown-linux-musl

FROM alpine:3.10
RUN apk add --no-cache bash
COPY --from=build /project/transanno/target/x86_64-unknown-linux-musl/release/transanno /usr/local/bin/