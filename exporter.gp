\\ ====================================================================
\\ CONFIGURATION: Set your number field defining polynomial here
\\ ====================================================================
\\ Example: x^3 - 3*x - 1 (a totally real cubic field)
P = x^3 - 3*x - 1; 

\\ ====================================================================
\\ PHASE 1: Heavy Number Theory Work (Done once)
\\ ====================================================================

nf = bnfinit(P);           \\ Computes unit group information

degree = poldegree(P);

r1     = nf.sign[1];
r2     = nf.sign[2];

\\ Extract the Integral Basis elements (as algebraic objects)
basis_elements = nf.zk;

\\ Extract Fundamental Units (expressed as elements in the field)
fundamental_units = nf.fu;
num_units = length(fundamental_units);

\\ ====================================================================
\\ PHASE 2: Geometric Embedding Transformation
\\ ====================================================================

\\ Helper function to extract and format embeddings of a given element
get_embeddings_matrix(element, nf_struct) = {
    my(embeds = nfeltembed(nf_struct, element));
    my(r1 = nf_struct.sign[1]);
    my(r2 = nf_struct.sign[2]);
    my(out_mat = matrix(r1 + r2, 2)); \\ Rows: each embedding, Columns: [Real, Imag]
    
    \\ Populate Real Embeddings
    for(i = 1, r1,
        out_mat[i, 1] = precision(real(embeds[i]), 12);
        out_mat[i, 2] = 0.0;
    );
    \\ Populate Complex Embeddings 
    for(j = 1, r2,
        my(c_val = embeds[r1 + j]);
        out_mat[r1 + j, 1] = precision(real(c_val), 12);
        out_mat[r1 + j, 2] = precision(imag(c_val), 12);
    );
    return(out_mat);
};

\\ ====================================================================
\\ PHASE 3: File Output Generation
\\ ====================================================================
filename = "field_data.txt";
file = fileopen(filename, "w");

filewrite(file, "# ==========================================");
filewrite(file, "# FIELD METADATA");
filewrite(file, "# ==========================================");
filewrite(file, Str("degree: ", degree));
filewrite(file, Str("r1: ", r1));
filewrite(file, Str("r2: ", r2));
filewrite(file, "");

filewrite(file, "# ==========================================");
filewrite(file, "# INTEGRAL BASIS EMBEDDINGS");
filewrite(file, "# Columns represent generators of O_K.");
filewrite(file, "# Row format: [Real_Part, Imag_Part]");
filewrite(file, "# (r1 rows have Imag=0, followed by r2 complex rows)");
filewrite(file, "# ==========================================");
filewrite(file, "INTEGRAL_BASIS_EMBEDDINGS:");

\\ Loop through the integral basis vectors and evaluate them at all embeddings
{
for(k=1, degree, 
    my(b_element = nf.zk[k]);
    my(emb_mat = get_embeddings_matrix(b_element, nf));
    filewrite(file, Str("# Generator ", k));
    for(row = 1, r1 + r2,
        filewrite(file, Str(strprintf("%.12f\t%.12f", emb_mat[row,1], emb_mat[row,2])));
    );
);
}

filewrite(file, "");

filewrite(file, "# ==========================================");
filewrite(file, "# FUNDAMENTAL UNITS");
filewrite(file, "# Each unit blocks out its layout identically");
filewrite(file, "# ==========================================");
filewrite(file, Str("FUNDAMENTAL_UNITS_COUNT: ", num_units));
filewrite(file, "");

\\ Loop through fundamental units and write out their raw embeddings
{
for(u = 1, num_units,
    my(unit_element = fundamental_units[u]);
    my(emb_mat = get_embeddings_matrix(unit_element, nf));
    filewrite(file, Str("UNIT_", u, "_EMBEDDINGS:"));
    for(row = 1, r1 + r2,
        filewrite(file, Str(strprintf("%.12f\t%.12f", emb_mat[row,1], emb_mat[row,2])));
    );
    filewrite(file, "");
);
}

fileclose(file);
print("Success! field_data.txt generated cleanly.");
quit;
