@doc doc"""
Minimise le problème : ``min_{||s||< \delta_{k}} q_k(s) = s^{t}g + (1/2)s^{t}Hs``
                        pour la ``k^{ème}`` itération de l'algorithme des régions de confiance

# Syntaxe
```julia
sk = Gradient_Conjugue_Tronque(fk,gradfk,hessfk,option)
```

# Entrées :   
   * **gradfk**           : (Array{Float,1}) le gradient de la fonction f appliqué au point xk
   * **hessfk**           : (Array{Float,2}) la Hessienne de la fonction f appliqué au point xk
   * **options**          : (Array{Float,1})
      - **delta**    : le rayon de la région de confiance
      - **max_iter** : le nombre maximal d'iterations
      - **tol**      : la tolérance pour la condition d'arrêt sur le gradient


# Sorties:
   * **s** : (Array{Float,1}) le pas s qui approche la solution du problème : ``min_{||s||< \delta_{k}} q(s)``

# Exemple d'appel:
```julia
gradf(x)=[-400*x[1]*(x[2]-x[1]^2)-2*(1-x[1]) ; 200*(x[2]-x[1]^2)]
hessf(x)=[-400*(x[2]-3*x[1]^2)+2  -400*x[1];-400*x[1]  200]
xk = [1; 0]
options = []
s = Gradient_Conjugue_Tronque(gradf(xk),hessf(xk),options)
```
"""
function Gradient_Conjugue_Tronque(gradfk,hessfk,options)

    "# Si option est vide on initialise les 3 paramètres par défaut"
    if options == []
        deltak = 2
        max_iter = 100
        tol = 1e-6
    else
        deltak = options[1]
        max_iter = options[2]
        tol = options[3]
    end
    
    n = length(gradfk)
    s_k = zeros(n)
    g_k = gradfk
    p_k = -gradfk
    
    for k in 1:max_iter
        k_k = (p_k)'*hessfk*p_k
        if k_k <= 0
            a = norm(p_k,2)^2
            b = 2 * s_k' * p_k
            c = norm(s_k,2)^2 - deltak^2
            
            if a == 0
                break
            end
            d = sqrt(b^2 - 4a*c)
            #q
            function q(s)
                return (gradfk)'*(s) + 1/2*(s)'*hessfk*(s)
            end
            
            rho_1,rho_2 = (-b - d)/(2*a),(-b + d)/(2*a)
            if q(s_k + rho_1*p_k) > q(s_k + rho_2*p_k) 
                rho_k = rho_2
            else
                rho_k = rho_1
            end
            
            s_k += rho_k*p_k
            break
        end
        
        alpha_k = (g_k)'*(g_k)/k_k
        
        if norm(s_k + alpha_k*p_k,2) >= deltak
            a = norm(p_k,2)^2
            b = 2 * s_k' * p_k
            c = norm(s_k,2)^2 - deltak^2
            if a == 0
                break
            end
            d = sqrt(b^2 - 4a*c)
            rho_1,rho_2 = (-b - d)/(2*a),(-b + d)/(2*a)
            if rho_1 > 0
                rho_k = rho_1
            elseif rho_2 > 0
                rho_k = rho_2
                s_k += rho_k * p_k
                break
            end
        end
        
        s_k += alpha_k * p_k
        g_k += alpha_k * hessfk * p_k
        beta_k = g_k'*g_k/((g_k - alpha_k * hessfk * p_k)'*(g_k - alpha_k * hessfk * p_k))
        p_k = -g_k + beta_k *p_k
        
        if norm(g_k,2) <= tol*norm(gradfk)
            break
        end
    end
   
   return s_k
end