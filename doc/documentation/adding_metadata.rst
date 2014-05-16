.. _adding_metadata:

====================================================
Adding sample and observation metadata to biom files
====================================================

Frequently you'll have an existing BIOM file and want to add sample and/or observation metadata to it. For samples, metadata is frequently environmental or technical details about your samples: the subject that a sample was collected from, the pH of the sample, the PCR primers used to amplify DNA from the samples, etc. For observations, metadata is frequently a categorization of the observation: the taxonomy of an OTU, or the EC hierarchy of a gene. You can use the ``biom add-metadata`` command to add this information to an existing BIOM file.

To get help with ``add-metadata`` you can call::

	biom add-metadata -h

This command takes a BIOM file, and corresponding sample and/or observation mapping files. The following examples are used in the commands below. You can find these files in the ``biom-format/examples`` directory.

Your BIOM file might look like the following::

    {
        "id":null,
        "format": "1.0.0",
        "format_url": "http://biom-format.org",
        "type": "OTU table",
        "generated_by": "some software package",
        "date": "2011-12-19T19:00:00",
        "rows":[
                {"id":"GG_OTU_1", "metadata":null},
                {"id":"GG_OTU_2", "metadata":null},
                {"id":"GG_OTU_3", "metadata":null},
                {"id":"GG_OTU_4", "metadata":null},
                {"id":"GG_OTU_5", "metadata":null}
            ],  
        "columns": [
                {"id":"Sample1", "metadata":null},
                {"id":"Sample2", "metadata":null},
                {"id":"Sample3", "metadata":null},
                {"id":"Sample4", "metadata":null},
                {"id":"Sample5", "metadata":null},
                {"id":"Sample6", "metadata":null}
            ],
        "matrix_type": "sparse",
        "matrix_element_type": "int",
        "shape": [5, 6], 
        "data":[[0,2,1],
                [1,0,5],
                [1,1,1],
                [1,3,2],
                [1,4,3],
                [1,5,1],
                [2,2,1],
                [2,3,4],
                [2,5,2],
                [3,0,2],
                [3,1,1],
                [3,2,1],
                [3,5,1],
                [4,1,1],
                [4,2,1]
               ]
    }

A sample metadata mapping file could then look like the following. Notice that there is an extra sample in here with respect to the above BIOM table. Any samples in the mapping file that are not in the BIOM file are ignored.
::

	#SampleID	BarcodeSequence	DOB
	# Some optional
	# comment lines...
	Sample1	AGCACGAGCCTA	20060805
	Sample2	AACTCGTCGATG	20060216
	Sample3	ACAGACCACTCA	20060109
	Sample4	ACCAGCGACTAG	20070530
	Sample5	AGCAGCACTTGT	20070101
	Sample6	AGCAGCACAACT	20070716

An observation metadata mapping file might look like the following. Notice that there is an extra observation in here with respect to the above BIOM table. Any observations in the mapping file that are not in the BIOM file are ignored.
::

	#OTUID	taxonomy	confidence
	# Some optional
	# comment lines
	GG_OTU_0	Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__	0.980
	GG_OTU_1	Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae	0.665
	GG_OTU_2	Root;k__Bacteria	0.980
	GG_OTU_3	Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae	1.000
	GG_OTU_4	Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae	0.842
	GG_OTU_5	Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae	1.000


Adding metadata
===============

To add sample metadata to a BIOM file, you can run the following::

	biom add-metadata -i min_sparse_otu_table.biom -o table.w_smd.biom --sample-metadata-fp sam_md.txt

To add observation metadata to a BIOM file, you can run the following::

	biom add-metadata -i min_sparse_otu_table.biom -o table.w_omd.biom --observation-metadata-fp obs_md.txt

You can also combine these in a single command to add both observation and sample metadata::

	biom add-metadata -i min_sparse_otu_table.biom -o table.w_md.biom --observation-metadata-fp obs_md.txt --sample-metadata-fp sam_md.txt

In the last case, the resulting BIOM file will look like the following::

	{
	    "columns": [
	        {
	            "id": "Sample1", 
	            "metadata": {
	                "BarcodeSequence": "AGCACGAGCCTA", 
	                "DOB": "20060805"
	            }
	        }, 
	        {
	            "id": "Sample2", 
	            "metadata": {
	                "BarcodeSequence": "AACTCGTCGATG", 
	                "DOB": "20060216"
	            }
	        }, 
	        {
	            "id": "Sample3", 
	            "metadata": {
	                "BarcodeSequence": "ACAGACCACTCA", 
	                "DOB": "20060109"
	            }
	        }, 
	        {
	            "id": "Sample4", 
	            "metadata": {
	                "BarcodeSequence": "ACCAGCGACTAG", 
	                "DOB": "20070530"
	            }
	        }, 
	        {
	            "id": "Sample5", 
	            "metadata": {
	                "BarcodeSequence": "AGCAGCACTTGT", 
	                "DOB": "20070101"
	            }
	        }, 
	        {
	            "id": "Sample6", 
	            "metadata": {
	                "BarcodeSequence": "AGCAGCACAACT", 
	                "DOB": "20070716"
	            }
	        }
	    ], 
	    "data": [
	        [0, 2, 1.0], 
	        [1, 0, 5.0], 
	        [1, 1, 1.0], 
	        [1, 3, 2.0], 
	        [1, 4, 3.0], 
	        [1, 5, 1.0], 
	        [2, 2, 1.0], 
	        [2, 3, 4.0], 
	        [2, 5, 2.0], 
	        [3, 0, 2.0], 
	        [3, 1, 1.0], 
	        [3, 2, 1.0], 
	        [3, 5, 1.0], 
	        [4, 1, 1.0], 
	        [4, 2, 1.0]
	    ], 
	    "date": "2012-12-11T07:36:15.467843", 
	    "format": "Biological Observation Matrix 1.0.0", 
	    "format_url": "http://biom-format.org", 
	    "generated_by": "some software package", 
	    "id": null, 
	    "matrix_element_type": "float", 
	    "matrix_type": "sparse", 
	    "rows": [
	        {
	            "id": "GG_OTU_1", 
	            "metadata": {
	                "confidence": "0.665", 
	                "taxonomy": "Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae"
	            }
	        }, 
	        {
	            "id": "GG_OTU_2", 
	            "metadata": {
	                "confidence": "0.980", 
	                "taxonomy": "Root;k__Bacteria"
	            }
	        }, 
	        {
	            "id": "GG_OTU_3", 
	            "metadata": {
	                "confidence": "1.000", 
	                "taxonomy": "Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae"
	            }
	        }, 
	        {
	            "id": "GG_OTU_4", 
	            "metadata": {
	                "confidence": "0.842", 
	                "taxonomy": "Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae"
	            }
	        }, 
	        {
	            "id": "GG_OTU_5", 
	            "metadata": {
	                "confidence": "1.000", 
	                "taxonomy": "Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae"
	            }
	        }
	    ], 
	    "shape": [5, 6], 
	    "type": "OTU table"
	}


Processing metadata while adding
================================

There are some additional parameters you can pass to this command for more complex processing. 

You can tell the command to process certain metadata column values as integers (``--int-fields``), floating point (i.e., decimal or real) numbers (``--float-fields``), or as hierarchical semicolon-delimited data (``--sc-separated``).

::

	biom add-metadata -i min_sparse_otu_table.biom -o table.w_md.biom --observation-metadata-fp obs_md.txt --sample-metadata-fp sam_md.txt --int-fields DOB --sc-separated taxonomy --float-fields confidence

Here your resulting BIOM file will look like the following, where ``DOB`` values are now integers (compare to the above: they're not quoted now), ``confidence`` values are now floating point numbers (again, not quoted now), and ``taxonomy`` values are now lists where each entry is a taxonomy level, opposed to above where they appear as a single semi-colon-separated string.
::

	{
	    "columns": [
	        {
	            "id": "Sample1", 
	            "metadata": {
	                "BarcodeSequence": "AGCACGAGCCTA", 
	                "DOB": 20060805
	            }
	        }, 
	        {
	            "id": "Sample2", 
	            "metadata": {
	                "BarcodeSequence": "AACTCGTCGATG", 
	                "DOB": 20060216
	            }
	        }, 
	        {
	            "id": "Sample3", 
	            "metadata": {
	                "BarcodeSequence": "ACAGACCACTCA", 
	                "DOB": 20060109
	            }
	        }, 
	        {
	            "id": "Sample4", 
	            "metadata": {
	                "BarcodeSequence": "ACCAGCGACTAG", 
	                "DOB": 20070530
	            }
	        }, 
	        {
	            "id": "Sample5", 
	            "metadata": {
	                "BarcodeSequence": "AGCAGCACTTGT", 
	                "DOB": 20070101
	            }
	        }, 
	        {
	            "id": "Sample6", 
	            "metadata": {
	                "BarcodeSequence": "AGCAGCACAACT", 
	                "DOB": 20070716
	            }
	        }
	    ], 
	    "data": [
	        [0, 2, 1.0], 
	        [1, 0, 5.0], 
	        [1, 1, 1.0], 
	        [1, 3, 2.0], 
	        [1, 4, 3.0], 
	        [1, 5, 1.0], 
	        [2, 2, 1.0], 
	        [2, 3, 4.0], 
	        [2, 5, 2.0], 
	        [3, 0, 2.0], 
	        [3, 1, 1.0], 
	        [3, 2, 1.0], 
	        [3, 5, 1.0], 
	        [4, 1, 1.0], 
	        [4, 2, 1.0]
	    ], 
	    "date": "2012-12-11T07:30:29.870689", 
	    "format": "Biological Observation Matrix 1.0.0", 
	    "format_url": "http://biom-format.org", 
	    "generated_by": "some software package", 
	    "id": null, 
	    "matrix_element_type": "float", 
	    "matrix_type": "sparse", 
	    "rows": [
	        {
	            "id": "GG_OTU_1", 
	            "metadata": {
	                "confidence": 0.665, 
	                "taxonomy": ["Root", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnospiraceae"]
	            }
	        }, 
	        {
	            "id": "GG_OTU_2", 
	            "metadata": {
	                "confidence": 0.98, 
	                "taxonomy": ["Root", "k__Bacteria"]
	            }
	        }, 
	        {
	            "id": "GG_OTU_3", 
	            "metadata": {
	                "confidence": 1.0, 
	                "taxonomy": ["Root", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnospiraceae"]
	            }
	        }, 
	        {
	            "id": "GG_OTU_4", 
	            "metadata": {
	                "confidence": 0.842, 
	                "taxonomy": ["Root", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnospiraceae"]
	            }
	        }, 
	        {
	            "id": "GG_OTU_5", 
	            "metadata": {
	                "confidence": 1.0, 
	                "taxonomy": ["Root", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnospiraceae"]
	            }
	        }
	    ], 
	    "shape": [5, 6], 
	    "type": "OTU table"
	}

If you have multiple fields that you'd like processed in one of these ways, you can pass a comma-separated list of field names (e.g., ``--float-fields confidence,pH``).

Renaming (or naming) metadata columns while adding
==================================================

You can also override the names of the metadata fields provided in the mapping files with the ``--observation-header`` and ``--sample-header`` parameters. This is useful if you want to rename metadata columns, or if metadata column headers aren't present in your metadata mapping file. If you pass either of these parameters, you must name all columns in order. If there are more columns in the metadata mapping file then there are headers, extra columns will be ignored (so this is also a useful way to select only the first n columns from your mapping file). For example, if you want to rename the ``DOB`` column in the sample metadata mapping you could do the following::
	
	biom add-metadata -i min_sparse_otu_table.biom -o table.w_smd.biom --sample-metadata-fp sam_md.txt --sample-header SampleID,BarcodeSequence,DateOfBirth

If you have a mapping file without headers such as the following::

	GG_OTU_0	Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__	0.980
	GG_OTU_1	Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae	0.665
	GG_OTU_2	Root;k__Bacteria	0.980
	GG_OTU_3	Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae	1.000
	GG_OTU_4	Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae	0.842
	GG_OTU_5	Root;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae	1.000

you could name these while adding them as follows::

	biom add-metadata -i min_sparse_otu_table.biom -o table.w_omd.biom --observation-metadata-fp obs_md.txt --observation-header OTUID,taxonomy,confidence

As a variation on the last command, if you only want to include the ``taxonomy`` column and exclude the ``confidence`` column, you could run::

	biom add-metadata -i min_sparse_otu_table.biom -o table.w_omd.biom --observation-metadata-fp obs_md.txt --observation-header OTUID,taxonomy

