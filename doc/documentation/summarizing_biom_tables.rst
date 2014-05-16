.. _summarizing_biom_tables:

====================================================
Summarizing BIOM tables
====================================================

If you have an existing BIOM file and want to compile a summary of the information in that table, you can use the ``biom summarize-table`` command.

To get help with ``biom summarize-table`` you can call::

	biom summarize-table -h

This command takes a BIOM file or gzipped BIOM file as input, and will print a summary of the count information on a per-sample basis to the new file specified by the ``-o`` parameter. The example file used in the commands below can be found in the ``biom-format/examples`` directory.

Summarizing sample data
-----------------------

To summarize the per-sample data in a BIOM file, you can run::
	
	biom summarize-table -i rich_sparse_otu_table.biom -o rich_sparse_otu_table_summary.txt

The following information will be written to ``rich_sparse_otu_table_summary.txt``::

	Num samples: 6
	Num observations: 5
	Total count: 27
	Table density (fraction of non-zero values): 0.500
	
	Counts/sample summary:
	 Min: 3.0
	 Max: 7.0
	 Median: 4.000
	 Mean: 4.500
	 Std. dev.: 1.500
	 Sample Metadata Categories: LinkerPrimerSequence; BarcodeSequence; Description; BODY_SITE
	 Observation Metadata Categories: taxonomy
	
	Counts/sample detail:
	 Sample5: 3.0
	 Sample2: 3.0
	 Sample6: 4.0
	 Sample3: 4.0
	 Sample4: 6.0
	 Sample1: 7.0

As you can see, general summary information about the table is provided, including the number of samples, the number of observations, the total count (i.e., the sum of all values in the table), and so on, followed by the per-sample counts.

Summarizing sample data qualitatively
--------------------------------------

To summarize the per-sample data in a BIOM file qualitatively, where the number of unique observations per sample (rather than the total count of observations per sample) are provided, you can run::

	biom summarize-table -i rich_sparse_otu_table.biom --qualitative -o rich_sparse_otu_table_qual_summary.txt

The following information will be written to ``rich_sparse_otu_table_qual_summary.txt``::

	Num samples: 6
	Num observations: 5
	
	Observations/sample summary:
	 Min: 1
	 Max: 4
	 Median: 2.500
	 Mean: 2.500
	 Std. dev.: 0.957
	 Sample Metadata Categories: LinkerPrimerSequence; BarcodeSequence; Description; BODY_SITE
	 Observation Metadata Categories: taxonomy
	
	Observations/sample detail:
	 Sample5: 1
	 Sample4: 2
	 Sample1: 2
	 Sample6: 3
	 Sample2: 3
