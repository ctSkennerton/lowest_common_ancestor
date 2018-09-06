from __future__ import print_function
import sys
import csv
from BioSQL import BioSeqDatabase
from biosqlx.taxon_tree import TaxonTree

conn = BioSeqDatabase.open_database(db='taxdb', driver='sqlite3')

tax_tree = TaxonTree(conn.adaptor)

taxon2id = {}
with open('biosqldb_annotations.csv') as fp:
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


