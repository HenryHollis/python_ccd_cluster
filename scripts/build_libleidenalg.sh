leiden_ccd_VERSION=0.11.0

ROOT_DIR=`pwd`
echo "Using root dir ${ROOT_DIR}"

# Create source directory
if [ ! -d "${ROOT_DIR}/build-deps/src" ]; then
  echo ""
  echo "Make directory ${ROOT_DIR}/build-deps/src"
  mkdir -p ${ROOT_DIR}/build-deps/src
fi

cd ${ROOT_DIR}/build-deps/src
if [ ! -d "leiden_ccd" ]; then
  echo ""
  echo "Cloning leiden_ccd into ${ROOT_DIR}/build-deps/src/leiden_ccd"
  # Clone repository if it does not exist yet
  # git clone --branch ${leiden_ccd_VERSION} https://github.com/HenryHollis/leiden_ccd.git --single-branch
    git clone https://github.com/HenryHollis/leiden_ccd.git --single-branch

fi

# # Make sure the git repository points to the correct version
# echo ""
# echo "Checking out ${leiden_ccd_VERSION} in ${ROOT_DIR}/build-deps/src/leiden_ccd"
# cd ${ROOT_DIR}/build-deps/src/leiden_ccd
# git fetch origin tag ${leiden_ccd_VERSION} --no-tags
# git checkout ${leiden_ccd_VERSION}

# Make build directory
if [ ! -d "${ROOT_DIR}/build-deps/build/leiden_ccd" ]; then
  echo ""
  echo "Make directory ${ROOT_DIR}/build-deps/build/leiden_ccd"
  mkdir -p ${ROOT_DIR}/build-deps/build/leiden_ccd
fi

# Configure, build and install
cd ${ROOT_DIR}/build-deps/build/leiden_ccd

echo ""
echo "Configure leiden_ccd build"
cmake ${ROOT_DIR}/build-deps/src/leiden_ccd \
    -DCMAKE_INSTALL_PREFIX=${ROOT_DIR}/build-deps/install/ \
    -DBUILD_SHARED_LIBS=ON \
    -Digraph_ROOT=${ROOT_DIR}/build-deps/install/lib/cmake/igraph/ \
    ${EXTRA_CMAKE_ARGS}

echo ""
echo "Build leiden_ccd"
cmake --build .

echo ""
echo "Install leiden_ccd to ${ROOT_DIR}/build-deps/install/"
cmake --build . --target install
