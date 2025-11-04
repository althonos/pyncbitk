from libcpp cimport bool
from libcpp.list cimport list
from libcpp.string cimport string

from ...serial.serialbase cimport CSerialObject
from ..seqfeat.org_ref cimport COrg_ref


cdef extern from "objects/taxon1/Taxon2_data_.hpp" namespace "ncbi::objects::CTaxon2_data_base" nogil:

    ctypedef COrg_ref TOrg
    ctypedef list[string] TBlast_name
    ctypedef bool TIs_uncultured
    ctypedef bool TIs_species_level


cdef extern from "objects/taxon1/Taxon2_data_.hpp" namespace "ncbi::objects" nogil:

    cppclass CTaxon2_data_Base(CSerialObject):
        CTaxon2_data_Base()

        bool IsSetOrg() const
        bool CanGetOrg() const
        void ResetOrg() except+
        const TOrg& GetOrg() except+
        TOrg& GetOrgMut "SetOrg" () except+
        void SetOrg(TOrg& value) except+

        bool IsSetBlast_name() const
        bool CanGetBlast_name() const
        void ResetBlast_name() except+
        const TBlast_name& GetBlast_name() except+
        # TBlast_name& SetBlast_name() except +

        bool IsSetIs_uncultured() const
        bool CanGetIs_uncultured() const
        void ResetIs_uncultured() except+
        TIs_uncultured GetIs_uncultured() except+
        void SetIs_uncultured(TIs_uncultured value) except+
        # TIs_uncultured& SetIs_uncultured() except +

        bool IsSetIs_species_level() const
        bool CanGetIs_species_level() const
        void ResetIs_species_level() except+
        TIs_species_level GetIs_species_level() except+
        void SetIs_species_level(TIs_species_level value) except+
        # TIs_species_level& SetIs_species_level() except+

        void Reset() except+


cdef extern from "objects/taxon1/Taxon2_data.hpp" namespace "ncbi::objects" nogil:

    cppclass CTaxon2_data(CTaxon2_data_Base):
        CTaxon2_data()

        # void SetProperty( const string& name, const string& value );
        # void SetProperty( const string& name, int value );
        # void SetProperty( const string& name, bool value );

        # bool GetProperty( const string& name, string& value ) const;
        # bool GetProperty( const string& name, int& value ) const;
        # bool GetProperty( const string& name, bool& value ) const;

        # void ResetProperty( const string& name );