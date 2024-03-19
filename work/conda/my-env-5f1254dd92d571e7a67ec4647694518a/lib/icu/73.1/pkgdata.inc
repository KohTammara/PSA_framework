GENCCODE_ASSEMBLY_TYPE=-a gcc
SO=so
SOBJ=so
A=a
LIBPREFIX=lib
LIB_EXT_ORDER=.73.1
COMPILE=/croot/icu_1692293016303/_build_env/bin/x86_64-conda-linux-gnu-cc -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /home/tammara/Documents/PSA_framework/work/conda/my-env-5f1254dd92d571e7a67ec4647694518a/include -D_REENTRANT  -DU_HAVE_ELF_H=1 -DU_HAVE_STRTOD_L=1 -DU_HAVE_XLOCALE_H=1 -DU_HAVE_STRING_VIEW=1  -DU_ATTRIBUTE_DEPRECATED= -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/tammara/Documents/PSA_framework/work/conda/my-env-5f1254dd92d571e7a67ec4647694518a/include -fdebug-prefix-map=/croot/icu_1692293016303/work=/usr/local/src/conda/icu-73.1 -fdebug-prefix-map=/home/tammara/Documents/PSA_framework/work/conda/my-env-5f1254dd92d571e7a67ec4647694518a=/usr/local/src/conda-prefix -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -c
LIBFLAGS=-I/home/tammara/Documents/PSA_framework/work/conda/my-env-5f1254dd92d571e7a67ec4647694518a/include -DPIC -fPIC
GENLIB=/croot/icu_1692293016303/_build_env/bin/x86_64-conda-linux-gnu-cc -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/tammara/Documents/PSA_framework/work/conda/my-env-5f1254dd92d571e7a67ec4647694518a/include -fdebug-prefix-map=/croot/icu_1692293016303/work=/usr/local/src/conda/icu-73.1 -fdebug-prefix-map=/home/tammara/Documents/PSA_framework/work/conda/my-env-5f1254dd92d571e7a67ec4647694518a=/usr/local/src/conda-prefix -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,/home/tammara/Documents/PSA_framework/work/conda/my-env-5f1254dd92d571e7a67ec4647694518a/lib -Wl,-rpath-link,/home/tammara/Documents/PSA_framework/work/conda/my-env-5f1254dd92d571e7a67ec4647694518a/lib -L/home/tammara/Documents/PSA_framework/work/conda/my-env-5f1254dd92d571e7a67ec4647694518a/lib  -shared -Wl,-Bsymbolic
LDICUDTFLAGS=-nodefaultlibs -nostdlib
LD_SONAME=-Wl,-soname -Wl,
RPATH_FLAGS=
BIR_LDFLAGS=-Wl,-Bsymbolic
AR=/croot/icu_1692293016303/_build_env/bin/x86_64-conda-linux-gnu-ar
ARFLAGS=r
RANLIB=/croot/icu_1692293016303/_build_env/bin/x86_64-conda-linux-gnu-ranlib
INSTALL_CMD=/usr/bin/install -c
