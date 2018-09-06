//
// much of this code is copied from:
// https://www.geeksforgeeks.org/lca-for-general-or-n-ary-trees-sparse-matrix-dp-approach-onlogn-ologn/
// Sparse Matrix DP approach to find LCA of two nodes
#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include <sstream>

using namespace std;
#define MAXN 3000000
#define level 24
 
vector <int> tree[MAXN];
int depth[MAXN];
int parent[MAXN][level];
 
// pre-compute the depth for each node and their
// first parent(2^0th parent)
// time complexity : O(n)
void dfs(int cur, int prev)
{
    depth[cur] = depth[prev] + 1;
    parent[cur][0] = prev;
    for (int i=0; i<tree[cur].size(); i++)
    {
        if (tree[cur][i] != prev)
            dfs(tree[cur][i], cur);
    }
}
 
// Dynamic Programming Sparse Matrix Approach
// populating 2^i parent for each node
// Time complexity : O(nlogn)
void precomputeSparseMatrix(int n)
{
    for (int i=1; i<level; i++)
    {
        for (int node = 1; node <= n; node++)
        {
            if (parent[node][i-1] != -1) {
                parent[node][i] =
                    parent[parent[node][i-1]][i-1];
                //cout << i << " " << node << " " << parent[node][i] << endl;
            }
        }
    }
}
 
// Returning the LCA of u and v
// Time complexity : O(log n)
int lca(int u, int v)
{
    if (depth[v] < depth[u])
        swap(u, v);
 
    int diff = depth[v] - depth[u];
 
    // Step 1 of the pseudocode
    for (int i=0; i<level; i++)
        if ((diff>>i)&1)
            v = parent[v][i];
 
    // now depth[u] == depth[v]
    if (u == v)
        return u;
 
    // Step 2 of the pseudocode
    for (int i=level-1; i>=0; i--)
        if (parent[u][i] != parent[v][i])
        {
            u = parent[u][i];
            v = parent[v][i];
        }
 
    return parent[u][0];
}
 
void addEdge(int u,int v)
{
    tree[u].push_back(v);
    tree[v].push_back(u);
}

void help() {
    cout << "lca -n nodes.tsv -b file.blast6 -a acc2taxid.tsv\n";
    cout << "Options:\n";
    cout << "-n FILE        file name containing 2 tab separated columns of numbers\n";
    cout << "               the first column is the node name and the second column is\n";
    cout << "               the parent node name\n";
    cout << "-b FILE        Tab separated 12 column file containing the blast results\n";
    cout << "               use -outfmt 6 in blast to get the right format\n";
    cout << "-a FILE        tab separated file containing 2 columns. The first column\n";
    cout << "               must be match the second column of the blast file, the\n";
    cout << "               second column must match to the first column of the nodes\n";
    cout << "               file (the taxon id)\n";
    cout << "-N FILE        Tab separated 2 column file containing the taxonomy id in the\n";
    cout << "               first column and the name of that taxonomy in the second column\n";
    exit(1);
}
struct options {
    string blastFile;
    string accessionFile;
    string nodeFile;
    string namesFile;
};

options parse_cmdline(int argc, char * argv[]) {
    int opt;
    options opts;

    while ((opt = getopt(argc, argv, "hn:b:a:N:")) != -1) {
        switch (opt) {
            case 'n':
                opts.nodeFile = optarg;
                break;
            case 'b':
                opts.blastFile = optarg;
                break;
            case 'a':
                opts.accessionFile = optarg;
                break;
            case 'N':
                opts.namesFile = optarg;
                break;
            case 'h':
            default: /* '?' */
                help();
                break;
        }
    }
    return opts;
}

struct blastResults {
    string query;
    vector<int> hitTaxids;
};

vector<blastResults> allBlastRes;

void process_blast(string& blastFile, unordered_map<string, int>& accessionMap){
    ifstream blastfs(blastFile);
    bool first = true;
    string query, subject, prev_query;
    float perc_id, bitscore, bitscore_limit;
    double evalue, evalue_limit;
    int length, qstart, qend, sstart, send, gap, mismatch;
    blastResults blast_res;

    // process first line
    blastfs >> query >> subject >> perc_id >> length >> mismatch >> gap >> qstart >> qend >> sstart >> send >>evalue >> bitscore;
    blast_res.query = query;
    blast_res.hitTaxids.push_back(accessionMap[subject]);
    evalue_limit = evalue * 0.9;
    bitscore_limit = bitscore * 0.9;

    while(blastfs >> query >> subject >> perc_id >> length >> mismatch >> gap >> qstart >> qend >> sstart >> send >>evalue >> bitscore ) {

        if (blast_res.query != query) {
            // we have reached a new query
            // prepare a new struct
            allBlastRes.push_back(blast_res);

            blast_res.query = query;
            blast_res.hitTaxids.erase(blast_res.hitTaxids.begin(), blast_res.hitTaxids.end());
            blast_res.hitTaxids.push_back(accessionMap[subject]);
            evalue_limit = evalue * 0.9;
            bitscore_limit = bitscore * 0.9;

        } else if (bitscore >= bitscore_limit) {
            blast_res.hitTaxids.push_back(accessionMap[subject]);
        }
    }
    // push back the last query
    allBlastRes.push_back(blast_res);
}

void process_accessions(string& accessionFile, unordered_map<string, int>& accessionMap) {
    string accession;
    int tax_id;
    ifstream accessionfs(accessionFile);
    while(accessionfs >> accession >> tax_id) {
        accessionMap[accession] = tax_id;
    }
}

void process_names(string& namesFile, unordered_map<int, string>& namesMap) {
    string line;
    ifstream namefs(namesFile);
    while(getline(namefs, line)) {
        stringstream ss;
        ss.str(line);
        string tax_id;
        string name;
        getline(ss, tax_id, '\t');
        getline(ss, name, '\t');
        int taxId = stoi(tax_id);
        namesMap[taxId] = name;
    }
}
 
// driver function
int main(int argc, char * argv[])
{
    memset(parent,-1,sizeof(parent));

    options opts;
    opts = parse_cmdline(argc, argv);
    ifstream nodesfs(opts.nodeFile);
    int child_id, parent_id, n;
    n = 1;
    while(nodesfs >> child_id >> parent_id) {
        addEdge(child_id, parent_id);
        n++; 
    }
    depth[0] = 0;
 
    // running dfs and precalculating depth
    // of each node.
    dfs(1,0);
 
    // Precomputing the 2^i th ancestor for evey node
    precomputeSparseMatrix(n);


    unordered_map<string, int> accession_map;
    accession_map.reserve(3000000);
    process_accessions(opts.accessionFile, accession_map);

    unordered_map<int, string> name_map;
    process_names(opts.namesFile, name_map);

    process_blast(opts.blastFile, accession_map);
 
    vector<blastResults>::iterator iter;
    for(iter = allBlastRes.begin(); iter != allBlastRes.end(); ++iter) {
        int lca_node = 0;
        vector<int>::iterator tax_iter;
        if (iter->hitTaxids.size() >= 2) {
            // get the lca of the first two nodes
            lca_node = lca(iter->hitTaxids[0], iter->hitTaxids[1]);
            //cerr << "first lca node is " << lca_node << "based on subjects " << iter->hitTaxids[0] << " "<< iter->hitTaxids[1]<<endl;
        }

        if (iter->hitTaxids.size() > 2) {
            for (tax_iter = iter->hitTaxids.begin()+3; tax_iter != iter->hitTaxids.end(); ++tax_iter) {
                //cerr << "refining lca node based on "<< *tax_iter << endl;
                lca_node = lca(lca_node, *tax_iter);
                //cerr << "new lca node is " << lca_node <<endl;
            }
        }
        cout << iter->query << '\t';
        if (lca_node > 1) {
            cout << name_map[lca_node];
            int curr_node = lca_node;
            while (parent[curr_node][0] > 1) {
                curr_node = parent[curr_node][0];
                cout << "\t" << name_map[curr_node];
            }
            cout <<"\t";
        }
        cout << "root" <<endl;
    }

    return 0;
}
