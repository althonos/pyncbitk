Object Model (``pyncbitk.objects``)
===================================

.. currentmodule:: pyncbitk.objects

.. automodule:: pyncbitk.objects


General (``pyncbitk.objects.general``)
--------------------------------------

.. currentmodule:: pyncbitk.objects.general

.. autosummary::

    ObjectId
    DBTag

.. toctree::
   :caption: General
   :maxdepth: 1
   :hidden:

   general


Sequence (``pyncbitk.objects.seq``)
-----------------------------------

.. currentmodule:: pyncbitk.objects.seq

.. autosummary::

    BioSeq

.. toctree::
   :caption: Sequence
   :maxdepth: 1
   :hidden:

   seq


Sequence Instance (``pyncbitk.objects.seqinst``)
------------------------------------------------

.. currentmodule:: pyncbitk.objects.seqinst

.. autosummary::

    SeqInst
    VirtualInst
    ContinuousInst
    SegmentedInst
    ConstructedInst
    RefInst
    ConsensusInst
    MapInst
    DeltaInst

.. toctree::
   :caption: Sequence Instance
   :maxdepth: 1
   :hidden:

   seqinst


Sequence Data (``pyncbitk.objects.seqdata``)
--------------------------------------------

.. currentmodule:: pyncbitk.objects.seqdata

.. autosummary::

    SeqData
    SeqAaData
    SeqNaData
    IupacNaData
    Ncbi2NaData
    Ncbi4NaData
    Ncbi8NaData
    NcbiPNaData
    IupacAaData
    Ncbi8AaData
    NcbiEAaData
    NcbiPAaData
    NcbiStdAa
    GapData

.. toctree::
   :caption: Sequence Data
   :maxdepth: 1
   :hidden:

   seqdata


Sequence Identifier (``pyncbitk.objects.seqid``)
------------------------------------------------

.. currentmodule:: pyncbitk.objects.seqid

.. autosummary::

    SeqId
    LocalId
    RefSeqId
    GenBankId
    ProteinDataBankId
    GeneralId
    OtherId
    TextSeqId

.. toctree::
   :caption: Sequence Identifier
   :maxdepth: 1
   :hidden:

   seqid


Sequence Description (``pyncbitk.objects.seqdesc``)
---------------------------------------------------

.. currentmodule:: pyncbitk.objects.seqdesc

.. autosummary::

    SeqDesc
    NameDesc
    TitleDesc
    RegionDesc
    SeqDescSet

.. toctree::
   :caption: Sequence Description
   :maxdepth: 1
   :hidden:

   seqdesc


Sequence Location (``pyncbitk.objects.seqloc``)
-----------------------------------------------

.. currentmodule:: pyncbitk.objects.seqloc

.. autosummary::

    SeqLoc
    NullLoc
    EmptySeqLoc
    WholeSeqLoc
    SeqIntervalLoc
    PackedSeqLoc
    PointLoc
    PackedPointsLoc
    MixLoc
    EquivalentLoc
    BondLoc
    FeatureLoc

.. toctree::
   :caption: Sequence Location
   :maxdepth: 1
   :hidden:

   seqloc


Sequence Alignments (``pyncbitk.objects.seqalign``)
---------------------------------------------------

.. currentmodule:: pyncbitk.objects.seqalign

.. autosummary::

    SeqAlign
    GlobalSeqAlign
    DiagonalSeqAlign
    PartialSeqAlign
    DiscontinuousSeqAlign
    AlignSegments
    DenseSegments
    DenseSegmentsData
    SeqAlignSet
    AlignRow
    SeqAlignScore

.. toctree::
   :caption: Sequence Location
   :maxdepth: 1
   :hidden:

   seqalign