q_in = [0.1,0.1,0.1];
q3 = q_in(1);
q4 = q_in(2);
q5 = q_in(3);

q_guess = [0, 0, q3, q4, q5, pi/4, pi/4, pi/4, 0.1, 0, 0, 0, 0, 0];

% initialize variables for keeping track of error
tol = 1e-10;
err = 2*tol;
max_it = 10;

% run no more than max_it iterations of updating the solution for m_qp
% exit loop once the error is below the input tolerance
it = 0;
while (it < max_it && err > tol) {
    psi = get_psi(q_guess_);

    VectorXd to_psi_bar(14);
    to_psi_bar << 0, 0, 0, 0, 0, 0, 0, 0, 0, q_guess(0), q_guess(1), q_guess(2), q_guess(3), q_guess(4);

    VectorXd psi_bar = psi - to_psi_bar;

    q_guess -= get_psi_dq(q_guess_).fullPivLu().solve(psi_bar);

    // update the error. this ends up being sqrt of sum of squares
    err = 0;
    for (size_t i = 0; i < 14; ++i) {
        err += psi(i)*psi(i);
    }
    err = sqrt(err);

    // while iterator
    it++;
}