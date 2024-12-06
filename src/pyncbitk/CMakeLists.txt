NCBI_add_subdirectory(objects)

cython_extension(serial   DEPENDS general xser LINKS ncbi_general)
# cython_extension(objmgr   DEPENDS seq xncbi xobjmgr)
cython_extension(objmgr   DEPENDS seq xncbi xobjmgr LINKS ncbi_general ncbi_seqext ncbi_core)
# cython_extension(objtools DEPENDS blastdb seqdb writedb xncbi xobjread pystreambuf)
# cython_extension(algo     DEPENDS blast seq xnetblast xncbi xblast)
