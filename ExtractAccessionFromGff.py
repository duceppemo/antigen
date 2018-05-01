#!/usr/local/env python3.5

import os

__author__ = 'duceppemo'


class ExtractAccessionFromGff(object):

    def __init__(self, args):
        # Define variables based on supplied arguments
        self.args = args
        self.input = args.input
        self.output = args.output
        self.db = self
        self.database = self.input + '.db'

        # run the script
        self.run()

    def run(self):
        self.parse_input()
        self.write_output()

    def parse_input(self):
        import gffutils

        # Check if database file already exists
        if os.path.exists(self.database):
            self.db = gffutils.FeatureDB(self.database, keep_order=True)
        else:
            self.db = gffutils.create_db(self.input, dbfn=self.database, force=True, keep_order=True,
                                         sort_attribute_values=True)

    def write_output(self):
        with open(self.output, 'w') as out:
            for cds in self.db.features_of_type('CDS'):
                try:  # In case some entries don't have 'Name' attributes
                    cds_name = cds['Name'][0]
                except KeyError:
                    cds_name = cds['ID'][0]

                for g in self.db.parents(self.db[cds], level=1):
                    gene_name = g['Name'][0]
                    # print(gene_name + ' -> ' + cds_name)
                    out.write(gene_name + "\t" + cds_name + "\n")

        # Remove database file
        os.remove(self.database)


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Read sequence information file and output in different format')
    parser.add_argument('-i', '--input', metavar='input.gff',
                        required=True,
                        help='Input file')
    parser.add_argument('-o', '--output', metavar='output.txt',
                        required=True,
                        help='Output file')

    # Get the arguments into an object
    arguments = parser.parse_args()

    ExtractAccessionFromGff(arguments)


