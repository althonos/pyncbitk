FROM quay.io/pypa/manylinux_2_34_x86_64:latest
ENV CONAN_HOME=/var/cache/conan2

ADD vendor /src/vendor
ADD conanfile.txt /src/conanfile.txt

RUN  yum install perl -y \
 &&  /opt/python/cp313-cp313/bin/pip install conan \
 &&  /opt/python/cp313-cp313/bin/conan export /src/vendor/ncbi-cxx-toolkit-conan --version 28.0.11 \
 &&  /opt/python/cp313-cp313/bin/conan profile detect \
 &&  /opt/python/cp313-cp313/bin/conan profile show \
 &&  /opt/python/cp313-cp313/bin/conan install /src/conanfile.txt --build=missing -s build_type=Release \
 &&  yum remove perl -y \
 &&  yum clean all -y

