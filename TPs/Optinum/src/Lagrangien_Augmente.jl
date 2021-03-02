function Lagrangien_Augmente(algo,fonc::Function,contrainte::Function,gradfonc::Function,
	hessfonc::Function,grad_contrainte::Function,hess_contrainte::Function,x0,options)

	if options == []
		epsilon = 1e-8
		tol = 1e-5
		itermax = 1000
		lambda0 = 2
		mu0 = 100
		tho = 2
	else
		epsilon = options[1]
		tol = options[2]
		itermax = options[3]
		lambda0 = options[4]
		mu0 = options[5]
		tho = options[6]
	end

    n = length(x0)
    xmin = zeros(n)
	fxmin = 0
	flag = 0
	iter = 0
    
    eta_cha = 0.1258925
    alpha = 0.1
    beta = 0.9
    epsilon_0 = 1
    epsilon_k = epsilon ###
    eta_k = eta_cha / (mu0^alpha)
    lambda_k = lambda0
    mu_k = mu0
    x_k = x0
    flag = 0
    iter = 0
    gradL0 = (gradfonc(x0) + lambda0'*grad_contrainte(x0) + mu0*grad_contrainte(x0)*contrainte(x0))
    
    while ( norm(gradfonc(x_k),2) > tol*(norm(gradfonc(x0),2)+ epsilon) )||  (norm(contrainte(x_k),2) > norm(contrainte(x0),2)*tol + epsilon ) && (iter < itermax)
                    
        function L(x)
            return fonc(x) + lambda_k'*contrainte(x)+mu_k/2*norm(contrainte(x),2)^2
        end
        
        function gradL(x)
            return gradfonc(x) + lambda_k'*grad_contrainte(x) + mu_k*grad_contrainte(x)*contrainte(x)
        end
        
        function hessL(x)
            return hessfonc(x) + lambda_k'*hess_contrainte(x) + mu_k*(hess_contrainte(x)*contrainte(x) + grad_contrainte(x)*grad_contrainte(x)')
        end
        
        # Calculer Approximation un minimiseur x_k du problème sans contraintes.
        if algo == "newton"
            x_k,~ = Algorithme_De_Newton(L,gradL,hessL,x_k,[])
        elseif algo =="cauchy"
            x_k,~ = Regions_De_Confiance("cauchy",L,gradL,hessL,x_k,[])
        elseif algo == "gct"
            x_k,~= Regions_De_Confiance("gct",L,gradL,hessL,x_k,[])
        else
            flag = -1
            break
        end
        
        
        
        if norm(gradL(x_k),2) <= tol*(norm(gradL0,2)+epsilon) && norm(contrainte(x_k)) <= (norm(contrainte(x0))*tol +epsilon)
            
            xmin = x_k
            break
        elseif norm(contrainte(x_k),2) <= eta_k
            lambda_k += mu_k*contrainte(x_k)
            epsilon_k = epsilon_k/mu_k
            eta_k = eta_k/(mu_k^beta)
        else
            mu_k =tho*mu_k
            epsilon_k = epsilon_0 / mu_k
            eta_k = eta_cha/(mu_k^alpha)
        end
        
        iter += 1
       if iter == itermax
            flag = 3
            break
        end
        
    end
    fxmin = fonc(xmin)
    return xmin,fxmin,flag,iter
end
