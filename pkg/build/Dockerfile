FROM quay.io/pypa/manylinux_2_34_x86_64:latest
ENV CONAN_HOME=/var/cache/conan2

ADD conanfile.txt                 /src/conanfile.txt
ADD vendor/ncbi-cxx-toolkit-conan /src/ncbi-cxx-toolkit-conan

RUN  yum install perl perl-IPC-Cmd -y \
 &&  /opt/python/cp313-cp313/bin/pip install conan \
 &&  /opt/python/cp313-cp313/bin/conan export /src/ncbi-cxx-toolkit-conan --version 29.6.0 \
 &&  /opt/python/cp313-cp313/bin/conan profile detect \
 &&  /opt/python/cp313-cp313/bin/conan profile show \
 &&  /opt/python/cp313-cp313/bin/conan install /src/conanfile.txt --build=missing -s compiler.cppstd=gnu20 -s build_type=Release \
 &&  yum remove perl -y \
 &&  yum clean all -y

