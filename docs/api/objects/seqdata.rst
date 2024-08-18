Sequence Data (``pyncbitk.objects.seqdata``)
============================================

.. currentmodule:: pyncbitk.objects.seqdata

.. automodule:: pyncbitk.objects.seqdata

Base Classes
------------

.. autoclass:: SeqData(Serial)
   :special-members: __init__
   :members:

.. autoclass:: SeqAaData(SeqData)
   :special-members: __init__
   :members:

.. autoclass:: SeqNaData(SeqData)
   :special-members: __init__
   :members:


Nucleotide Data
---------------

.. autoclass:: IupacNaData(SeqNaData)
   :special-members: __init__
   :members:

.. autoclass:: Ncbi2NaData(SeqNaData)
   :special-members: __init__
   :members:

.. autoclass:: Ncbi4NaData(SeqNaData)
   :special-members: __init__
   :members:

.. autoclass:: Ncbi8NaData(SeqNaData)
   :special-members: __init__
   :members:

.. autoclass:: NcbiPNaData(SeqNaData)
   :special-members: __init__
   :members:


Protein Data
------------

.. autoclass:: IupacAaData(SeqAaData)
   :special-members: __init__
   :members:

.. autoclass:: Ncbi8AaData(SeqNaData)
   :special-members: __init__
   :members:

.. autoclass:: NcbiEAaData(IupacAaData)
   :special-members: __init__
   :members:

.. autoclass:: NcbiPAaData(SeqNaData)
   :special-members: __init__
   :members:

.. autoclass:: NcbiStdAa(SeqNaData)
   :special-members: __init__
   :members:


Gaps
----

.. autoclass:: GapData(SeqData)
   :special-members: __init__
   :members: