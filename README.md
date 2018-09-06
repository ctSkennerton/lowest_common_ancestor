# Lowest Common Ancestor
This is a set of code for finding the lowest common ancestor for taxonomic
predictions. It was made because I had a number of organisms that weren't
in NCBI and wanted to include their novel taxonomies in my predictions. To
do this we first require an integrated hierarchy with both NCBI sequences
and any new taxonomies.

# Dependancies
You'll need a c++ compiler that handles c++11, sqlite3, python3, biopython and biosql-extensions to be installed

# Compile `lca`
```
g++ -std=c++11 -o lca lca.cpp
```

# Workflow
This isn't a complete stand alone script but below are the steps to get things going.

1. Make an sqlite3 database using the provided file `taxdb.sql`
   ```   
   sqlite3 taxdb < taxdb.sql
   ```
2. Load the NCBI taxonomy using the `load_ncbi_taxonomy.pl` script. This script is part of bioSQL and 
   can automatically download the taxonomy database from NCBI if you don't have it already
   ```
	perl load_ncbi_taxonomy.pl --driver SQLite --dbname taxdb
   ```
3. Run blast or diamond or any other program that creates a blast-style tabular output.
4. Get just the reference sequences that had a match from your query sample
   ```
   cut -f 2 <blast_style_tabular_file.tsv> | sort -u > accessions.txt
   ```
5. Get the taxonomic information for any NCBI sequence that had a hit. The file "prot.accession2taxid" can be downloaded from NCBI and contains the full set but we can speed things up by just taking the subset of proteins that actually matched our query.
   ```
   fgrep -wf accessions.txt prot.accession2taxid > prot.accession2taxid.subset
   ```
6. Load in the NCBI accessions
   ```
   python load_acc_taxon.py -d taxdb -i prot.accession2taxid.subset
   ```
7. Load in novel taxonomies using `add_taxonomy.py`. Provide a csv file containing protein IDs and their taxonomy if they don't exist in NCBI.
   ```	
   python add_taxonomy.py -d taxdb -i <file_with_novel_proteins_and_taxonomy.csv>
   ```
8. Now with the integrated set of taxonomy information from NCBI and your personal genomes, make the files for running the LCA algorithm
   ```
   sqlite3 taxdb "select taxon_id, parent_taxon_id from taxon;" | tr "|" "\t" > nodes.tsv
   sqlite3 taxdb "select taxon_id, name from taxon_name where name_class = 'scientific name';" | tr "|" "\t" > names.tsv
   sqlite3 taxdb "select * from acc_taxon;" | tr "|" "\t" > acc2tax.tsv
   ```
9. Run LCA script
   ```
   lca -b <blast_style_tabular_file.tsv> -n nodes.tsv -N names.tsv -a acc2tax.tsv > lca.txt
   ```

