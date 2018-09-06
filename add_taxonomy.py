from __future__ import print_function
import sys
import csv
from BioSQL import BioSeqDatabase
from biosqlx.taxon_tree import TaxonTree
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--database", help="path to the sqlite3 database prepopulated with the NCBI taxonomy")
    parser.add_argument("-i", "--input", help="csv file, with header row containing 'taxonomy' and 'seqfeature' headers, containing the taxonomies of genes/proteins that are not in NCBI. The taxonomy should be a semicolon (';') separated list.")
    args = parser.parse_args()

    conn = BioSeqDatabase.open_database(db=args.database, driver='sqlite3')

    tax_tree = TaxonTree(conn.adaptor)

    taxon2id = {}
    with open(args.input) as fp:
        parser = csv.DictReader(fp)
        for row in parser:
            try:
                taxon2id[row['taxonomy']].add(row['seqfeature'])
            except KeyError:
                taxon2id[row['taxonomy']] = set([row['seqfeature']])

    rows = conn.adaptor.execute_and_fetchall("select max(ncbi_taxon_id) from taxon")
    max_tax_id = int(rows[0][0])


    for taxon, seqfeatures in taxon2id.items():
        parts = taxon.split(';')
        parent = None
        for part in parts:
            nodes = tax_tree.find_elements(name=part)
            nodes2 = []
            for n in nodes:
                if n.__dict__['scientific name'][0] == part:
                    nodes2.append(n)

            nodes = nodes2

            if len(nodes) > 1:
                raise ValueError("There are more than one taxon associated with {}\n{}".format(part, [x.name for x in nodes]))
            elif len(nodes) == 0:
                # not found. Everything after this must be added in
                print(parent, part)
                max_tax_id += 1
                parent = tax_tree.add(part, 'scientific name', ncbi_taxon_id=max_tax_id, parent=parent)
            else:
                # it found. add in the parent
                parent = nodes[0]

        inserts = [ (x, parent._id) for x in seqfeatures ]
        conn.adaptor.executemany('INSERT INTO acc_taxon VALUES (?, ?)', inserts)

    conn.commit()


