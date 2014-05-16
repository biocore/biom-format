.. _converting:

===============================
Converting between file formats
===============================

The ``convert`` command in the biom-format project can be used to convert between biom and tab-delimited table formats. This is useful for several reasons:
 - converting biom format tables to tab-delimited tables for easy viewing in programs such as Excel
 - converting between sparse and dense biom formats

 .. note:: The tab-delimited tables are commonly referred to as the `classic format` tables, while BIOM formatted tables are referred to as `biom tables`.

General usage examples
----------------------

Convert a tab-delimited table to biom format. Note that you *must* specify the type of table here::

	biom convert -i table.txt -o table.from_txt.biom --table-type="otu table"

Convert biom format to tab-delimited table format::

	biom convert -i table.biom -o table.from_biom.txt -b

Convert biom format to classic format, including the ``taxonomy`` observation metadata as the last column of the classic format table. Because the BIOM format can support an arbitrary number of observation (or sample) metadata entries, and the classic format can support only a single observation metadata entry, you must specify which of the observation metadata entries you want to include in the output table::

	biom convert -i table.biom -o table.from_biom_w_taxonomy.txt -b --header-key taxonomy

Convert biom format to classic format, including the ``taxonomy`` observation metadata as the last column of the classic format table, but renaming that column as ``ConsensusLineage``. This is useful when using legacy tools that require a specific name for the observation metadata column.::

	biom convert -i table.biom -o table.from_biom_w_consensuslineage.txt -b --header-key taxonomy --output-metadata-id "ConsensusLineage"

Special case usage examples
---------------------------

Converting QIIME 1.4.0 and earlier OTU tables to BIOM format
````````````````````````````````````````````````````````````
If you are converting a QIIME 1.4.0 or earlier OTU table to BIOM format, there are a few steps to go through. First, for convenience, you might want to rename the ``ConsensusLineage`` column ``taxonomy``. You can do this with the following command::

	sed 's/Consensus Lineage/ConsensusLineage/' < otu_table.txt | sed 's/ConsensusLineage/taxonomy/' > otu_table.taxonomy.txt

Then, you'll want to perform the conversion including a step to convert the taxonomy `string` from the classic OTU table to a taxonomy `list`, as it's represented in QIIME 1.4.0-dev and later::

	biom convert -i otu_table.taxonomy.txt -o otu_table.from_txt.biom --table-type="otu table" --process-obs-metadata taxonomy
