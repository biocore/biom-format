#!/usr/bin/env python

from biom.parse import get_axis_indices, direct_slice_data, direct_parse_key

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, The BIOM-Format Project"
__credits__ = ["Daniel McDonald"]
__url__ = "http://biom-format.org"
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

try:
    from cogent.util.option_parsing import parse_command_line_parameters, \
            make_option
    cogent_cl_parsing = True
except ImportError:
    from sys import argv
    cogent_cl_parsing = False

if cogent_cl_parsing:
    script_info = {}
    script_info['brief_description'] = "Subset a BIOM file."
    script_info['script_description'] = "Subset a BIOM file without fully parsing it"
    script_info['script_usage'] = [("","Subset the observations in my_data.biom file.","%prog -i my_data.biom -a observations -s file_with_ids")]
    script_info['output_description']= ""
    script_info['required_options'] = [
     make_option('-i','--biom_fp',type="existing_filepath",
                 help='the BIological Observation Matrix filepath'),
     make_option('-a','--axis', type='choice',
                  choices=['observations','samples'],
                  help="The axis to subset over"),
     make_option('-s','--ids_fp',type="existing_filepath",
                 help="A file containing a single column of IDs to retain"),
     make_option('-o','--output_fp',type="new_filepath",
                 help="A file to write the result to")
    ]
    script_info['version'] = __version__
else:
    from optparse import OptionParser, make_option
    options = [
     make_option('-i','--biom_fp',type="string",
                 help='the BIological Observation Matrix filepath'),
     make_option('-a','--axis', type='string',
                  help="The axis to subset over"),
     make_option('-s','--ids_fp',type="string",
                 help="A file containing a single column of IDs to retain"),
     make_option('-o','--output_fp',type="string",
                 help="A file to write the result to")
    ]
    
if __name__ == '__main__':
    if cogent_cl_parsing:
        option_parser, opts, args =\
                     parse_command_line_parameters(**script_info)
    else:
        parser = OptionParser(option_list=options)
        opts, args = parser.parse_args()

    ids = [l.strip() for l in open(opts.ids_fp)]
    biom_str = open(opts.biom_fp).read()

    idxs, new_axis_md = get_axis_indices(biom_str, ids, opts.axis)
    new_data = direct_slice_data(biom_str, idxs, opts.axis)
    output = open(opts.output_fp,'w')

    # multiple walks over the file. bad form, but easy right now
    # ...should add a yield_and_ignore parser or something.
    output.write('{')
    output.write(direct_parse_key(biom_str, "id"))
    output.write(",")
    output.write(direct_parse_key(biom_str, "format"))
    output.write(",")
    output.write(direct_parse_key(biom_str, "format_url"))
    output.write(",")
    output.write(direct_parse_key(biom_str, "type"))
    output.write(",")
    output.write(direct_parse_key(biom_str, "generated_by"))
    output.write(",")
    output.write(direct_parse_key(biom_str, "date"))
    output.write(",")
    output.write(direct_parse_key(biom_str, "matrix_type"))
    output.write(",")
    output.write(direct_parse_key(biom_str, "matrix_element_type"))
    output.write(",")
    output.write(new_data)
    output.write(",")
    output.write(new_axis_md)
    output.write(",")
    if opts.axis == "observations":
        output.write(direct_parse_key(biom_str, "columns"))
    else:
        output.write(direct_parse_key(biom_str, "rows"))
    output.write("}")
    output.close()
    


