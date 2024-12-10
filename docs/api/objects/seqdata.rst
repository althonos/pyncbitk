Sequence Data (``pyncbitk.objects.seqdata``)
============================================

.. currentmodule:: pyncbitk.objects.seqdata

.. automodule:: pyncbitk.objects.seqdata

Base Classes
------------

.. autoclass:: SeqData(Serial)
   :special-members:
   :members:

.. autoclass:: SeqAaData(SeqData)
   :special-members: 
   :members:

.. autoclass:: SeqNaData(SeqData)
   :special-members: 
   :members:


Nucleotide Data
---------------

.. autoclass:: IupacNaData(SeqNaData)
   :special-members: 
   :members:

.. autoclass:: Ncbi2NaData(SeqNaData)
   :special-members: 
   :members:

.. autoclass:: Ncbi4NaData(SeqNaData)
   :special-members: 
   :members:

.. autoclass:: Ncbi8NaData(SeqNaData)
   :special-members: 
   :members:

.. autoclass:: NcbiPNaData(SeqNaData)
   :special-members: 
   :members:


Protein Data
------------

.. autoclass:: IupacAaData(SeqAaData)
   :special-members: 
   :members:

.. autoclass:: Ncbi8AaData(SeqNaData)
   :special-members: 
   :members:

.. autoclass:: NcbiEAaData(IupacAaData)
   :special-members: 
   :members:

.. autoclass:: NcbiPAaData(SeqNaData)
   :special-members: 
   :members:

.. autoclass:: NcbiStdAa(SeqNaData)
   :special-members: 
   :members:


Gaps
----

.. autoclass:: GapData(SeqData)
   :special-members: 
   :members: