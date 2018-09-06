from __future__ import print_function
import sqlite3

conn = sqlite3.connect('taxdb')
c = conn.cursor()
c.execute('PRAGMA foreign_keys = ON')
c.execute('PRAGMA journal_mode=WAL')

problems = []
with open('prot.accession2taxid.subset') as fp:
    next(fp)
    inserts = []
    line_number = 0
    for line in fp:
        fields = line.split()
        acc = fields[1]
        taxid = int(fields[2])
        c.execute('SELECT taxon_id FROM taxon WHERE ncbi_taxon_id = ?', (taxid,))
        tax_info = c.fetchone()
        if tax_info is None:
            problems.append(line)
            #raise RuntimeError("Could not find %i for %s in database, this likely means that the taxonomy db is out of date" % (taxid, acc) )
        else:
            inserts.append((acc, tax_info[0]))

        if len(inserts) >= 100000:
            line_number += 100000
            try:
                c.executemany('INSERT INTO acc_taxon VALUES (?, ?)', inserts)
            except sqlite3.IntegrityError:
                print( acc, taxid)
                raise
            else:
                conn.commit()
                del inserts[:]
                print(line_number)

    if len(inserts):
        try:
            c.executemany('INSERT INTO acc_taxon VALUES (?, ?)', inserts)
        except sqlite3.IntegrityError:
            print(acc, taxid)
            raise

conn.commit()

if len(problems) > 0:
    print("identified %i accessions that are missing in the taxonomy information" % len(problems))
    print("they will be saved to the file 'problems.txt'")
    with open("problems.txt", 'w') as fp:
        for i in problems:
            fp.write(i)
