\\ Get polynomial P from environment variable
env_str = getenv("POLY"); \\ 1. Read the text representation from the system environment

\\ 2. Safety check: make sure the variable actually exists
{
if (env_str == "0", 
    print("Error: The POLY environment variable is not set."); quit;
);
}
P = eval(env_str); \\ 3. Convert the string into a live algebraic polynomial object

nf = bnfinit(P);
degree = poldegree(P);
r1 = nf.sign[1];
r2 = nf.sign[2];
basis_elements = nf.zk;
fundamental_units = nf.fu;
num_units = length(fundamental_units);

get_embeddings_matrix(element, nf_struct) = {
    my(embeds = nfeltembed(nf_struct, element));
    my(r1 = nf_struct.sign[1]);
    my(r2 = nf_struct.sign[2]);
    my(out_mat = matrix(r1 + r2, 2));
    for(i = 1, r1,
        out_mat[i, 1] = precision(real(embeds[i]), 12);
        out_mat[i, 2] = 0.0;
    );
    for(j = 1, r2,
        my(c_val = embeds[r1 + j]);
        out_mat[r1 + j, 1] = precision(real(c_val), 12);
        out_mat[r1 + j, 2] = precision(imag(c_val), 12);
    );
    return(out_mat);
};

filename = "field_data.txt";
file = fileopen(filename, "w");
filewrite(file, "# FIELD METADATA");
filewrite(file, Str("polynomial: ", P));
filewrite(file, Str("degree: ", degree));
filewrite(file, Str("r1: ", r1));
filewrite(file, Str("r2: ", r2));
filewrite(file, Str("discriminant: ", nf.disc));
filewrite(file, "");
filewrite(file, "INTEGRAL_BASIS_EMBEDDINGS:");
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
filewrite(file, Str("FUNDAMENTAL_UNITS_COUNT: ", num_units));
filewrite(file, "");
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
print("Generated field_data.txt cleanly.");
quit;
