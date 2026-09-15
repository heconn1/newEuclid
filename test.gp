{
\\ 1. Read the text representation from the system environment
env_str = getenv("POLY");

\\ 2. Safety check: make sure the variable actually exists
if (env_str == "0", 
    print("Error: The POLY environment variable is not set."); quit;
);

\\ 3. Convert the string into a live algebraic polynomial object
P = eval(env_str);

\\ 4. Print your result
print("poly = ", P);
quit;
}
