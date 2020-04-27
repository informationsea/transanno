# This script takes care of testing your crate

set -ex

# TODO This is the "test phase", tweak it as you see fit
main() {
    ./prepare-test-files.sh  

    rustup target add $TARGET
    cargo build --target $TARGET --release

    if [ ! -z $DISABLE_TESTS ]; then
        return
    fi

    pushd liftover-rs
    cargo test --target $TARGET --release
    popd    

    cargo test --target $TARGET --release
    rm -rf target/test-output
}

# we don't run the "test phase" when doing deploys
if [ -z $TRAVIS_TAG ]; then
    main
fi
