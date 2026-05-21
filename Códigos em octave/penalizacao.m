function pen = penalizacao(params)

    Vinf = params(1);
    A = params(2);
    tau = params(3);
    p = params(4);

    pen = 0;

    if Vinf <= 0
        pen = pen + 1e6 + 1e6*abs(Vinf);
    end

    if A <= 0
        pen = pen + 1e6 + 1e6*abs(A);
    end

    if tau <= 0
        pen = pen + 1e6 + 1e6*abs(tau);
    end

    if p <= 0
        pen = pen + 1e6 + 1e6*abs(p);
    end

end
