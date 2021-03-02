function Pas_De_Cauchy(g,H,delta)
    
    t = 1/2 * g' * H * g
    norm_g = norm(g,2)
    b = -norm_g^2
    
    if (t > 0)
        l = -b/(2*t)
        s = -l * g ;
        e = 1
        if norm(s,2) > delta
            e = -1
        end
    else
        l = 1
        s = -l * g
        e = 0
    end
    return s,e
end
