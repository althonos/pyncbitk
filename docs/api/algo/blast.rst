BLAST (Basic Local Alignment Search Tool)
=========================================

.. currentmodule:: pyncbitk.algo

Runners
-------

.. autoclass:: Blast
   :special-members: __init__
   :members:

.. autoclass:: NucleotideBlast(Blast)
   :special-members: __init__
   :members:

.. autoclass:: ProteinBlast(Blast)
   :special-members: __init__
   :members:

.. autoclass:: MappingBlast(Blast)
   :special-members: __init__
   :members:

.. autoclass:: BlastP(ProteinBlast)
   :special-members: __init__
   :members:

.. autoclass:: BlastN(NucleotideBlast)
   :special-members: __init__
   :members:

.. autoclass:: TBlastN(ProteinBlast)
   :special-members: __init__
   :members:


Results
-------

.. autoclass:: SearchResultsSet
   :special-members: __init__
   :members:

.. autoclass:: SearchResults
   :special-members: __init__
   :members: