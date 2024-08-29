data {
    int<lower=1> dim;
}
parameters {
    vector[dim] x;
}
model {
    x ~ std_normal();
}
