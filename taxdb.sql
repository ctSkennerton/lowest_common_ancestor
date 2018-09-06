CREATE TABLE taxon (
       taxon_id		INTEGER PRIMARY KEY,
       ncbi_taxon_id 	INT(10),
       parent_taxon_id	INT(10) ,
       node_rank	VARCHAR(32),
       genetic_code	TINYINT ,
       mito_genetic_code TINYINT ,
       left_value	INT(10) ,
       right_value	INT(10) ,
       UNIQUE (ncbi_taxon_id),
       UNIQUE (left_value),
       UNIQUE (right_value)
);

CREATE INDEX taxparent ON taxon(parent_taxon_id);

CREATE TABLE taxon_name (
       taxon_id		INTEGER,
       name		VARCHAR(255)  NOT NULL,
       name_class	VARCHAR(32)  NOT NULL,
    UNIQUE (taxon_id,name,name_class),
    FOREIGN KEY ( taxon_id ) REFERENCES taxon ( taxon_id ) ON DELETE CASCADE
);

CREATE INDEX taxnametaxonid ON taxon_name(taxon_id);
CREATE INDEX taxnamename    ON taxon_name(name);

CREATE TABLE acc_taxon (
	accession  TEXT,
	taxon_id   INTEGER,
	UNIQUE(accession, taxon_id),
	FOREIGN KEY( taxon_id ) REFERENCES taxon (taxon_id )
);

CREATE INDEX acctaxonaccession ON acc_taxon(accession);
