from libcpp cimport bool

from .ddumpable cimport CDebugDumpable

cdef extern from "corelib/ncbiobj.hpp" namespace "ncbi" nogil:

    cppclass CObject(CDebugDumpable):
        CObject()
        CObject(const CObject& src)

    cppclass CRef[C]:
        ctypedef C element_type
        ctypedef element_type TObjectType

        CRef() noexcept
        CRef(TObjectType* ptr)
        CConstRef(CRef[C]& ref)

        operator!() noexcept

        bool Empty() noexcept
        bool NotEmpty() noexcept
        bool IsNull() noexcept
        bool NotNull() noexcept
        
        void Swap(CRef[C]& ref)
        void Reset()
        void Reset(TObjectType* newPtr)
        C* ReleaseOrNull()
        C* Release()
        void AtomicResetFrom(const CRef[C]& ref)
        void AtomicReleaseTo(const CRef[C]& ref)

        C* GetNonNullPointer() except +
        C* GetPointerOrNull() noexcept
        C* GetPointer() noexcept
        C& GetObject() except +
        C& operator*() except +

        C& operator* () except +
       

    cppclass CConstRef[C]:
        ctypedef C element_type
        ctypedef element_type TObjectType

        CConstRef() noexcept
        CConstRef(TObjectType* ptr)
        CConstRef(CRef[C]& ref)
        CConstRef(CConstRef[C]& ref)

        bool Empty() noexcept
        bool NotEmpty() noexcept
        bool IsNull() noexcept
        bool NotNull() noexcept

        void Swap(CConstRef[C]& ref)
        void Reset()
        void Reset(TObjectType* newPtr)
        C* ReleaseOrNull()
        C* Release()
        void AtomicResetFrom(const CConstRef[C]& ref)
        void AtomicReleaseTo(const CConstRef[C]& ref)

        const C* GetNonNullPointer() except +
        const C* GetPointerOrNull() noexcept
        const C* GetPointer() const
        C& GetObject() except +

        TObjectType& operator*() const