BLAST (Basic Local Alignment Search Tool)
=========================================

.. currentmodule:: pyncbitk.algo.blast

Runners
-------

.. autoclass:: Blast
   :special-members:
   :members:

.. autoclass:: NucleotideBlast(Blast)
   :special-members:
   :members:

.. autoclass:: ProteinBlast(Blast)
   :special-members:
   :members:

.. autoclass:: MappingBlast(Blast)
   :special-members:
   :members:

.. autoclass:: BlastN(NucleotideBlast)
   :special-members:
   :members:

.. autoclass:: BlastP(ProteinBlast)
   :special-members:
   :members:

.. autoclass:: BlastX(NucleotideBlast)
   :special-members:
   :members:

.. autoclass:: TBlastN(ProteinBlast)
   :special-members:
   :members:


Queries
-------

.. autoclass:: SearchQuery
   :special-members:
   :members:

.. autoclass:: SearchQueryVector
   :special-members:
   :members:


Results
-------

.. autoclass:: SearchResultsSet
   :special-members:
   :members:

.. autoclass:: SearchResults
   :special-members:
   :members: