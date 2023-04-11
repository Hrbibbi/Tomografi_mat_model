function noised_vector = add_noise(vec,sigma)
    n=normrnd(0,sigma,[size(vec)]);
    noised_vector=n+vec;
