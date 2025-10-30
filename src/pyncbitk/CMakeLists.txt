add_subdirectory(objects)
add_subdirectory(algo)

# cython_extension(serial   DEPENDS general xser)
# cython_extension(objmgr   DEPENDS seq xncbi xobjmgr)
# cython_extension(objtools DEPENDS blastdb seqdb writedb xncbi xobjread pystreambuf)
# cython_extension(algo     DEPENDS ncbi_algo blast seq xnetblast xncbi xblast xblastformat pystreambuf)

cython_extension(serial)
cython_extension(objmgr)
cython_extension(objtools DEPENDS pystreambuf)

