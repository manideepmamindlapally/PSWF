function [psi, lambda] = pswf(c,T,tt,Fs,N)
    t = -tt:1/Fs:tt;
    mid_t_idx = floor(Fs*(tt-T/2))+1:floor(Fs*(tt+T/2));
    
    no = 0:N-1;
    even_no = 0:2:N-1;
    odd_no = 1:2:N-1;
    FACT = factorial(0:N+2);
    %=========================================================================
    %===========EVALUATION of 'd's for even n ================================
    alpha = c^2 * ( (no+2).*(no+1) ) ./ ( (2*no+5).*(2*no+3) );
    beta  = no.*(no+1) + c^2 * ( 2*no.*(no+1)-1 )./( (2*no-1).*(2*no+3) ); 
    gamma = c^2 * ( no.*(no-1) ) ./ ( (2*no-3).*(2*no-1) );

    MAT = zeros(N/2,N/2);
    for k = 0:2:N-3
        MAT(k/2+1,k/2+2) = alpha(k+1);
    end
    for k = 0:2:N-1
        MAT(k/2+1,k/2+1) = beta(k+1);
    end
    for k = 2:2:N-1
        MAT(k/2+1,k/2) = gamma(k+1);
    end

    [D_even, lam_even] = eig(MAT, 'vector');

    X1 = ( power(-1,even_no/2) .* FACT(even_no+1) ./ power(2,even_no) ./ FACT(even_no/2+1).^2  ) * D_even;
    X2 = power(-1,even_no/2) .* FACT(even_no+1) ./ power(2,even_no) ./ FACT(even_no/2+1).^2;
    multifact = X2 ./ X1;
    D_even = D_even .* multifact;

    %=========================================================================
    %===========EVALUATION of 'd's for odd n ================================
    MAT = zeros(N/2,N/2);
    for k = 1:2:N-3
        MAT((k+1)/2,(k+1)/2+1) = alpha(k+1);
    end
    for k = 1:2:N-1
        MAT((k+1)/2,(k+1)/2) = beta(k+1);
    end
    for k = 3:2:N-1
        MAT((k+1)/2,(k+1)/2-1) = gamma(k+1);
    end

    [D_odd, lam_odd] = eig(MAT, 'vector');

    X1 = ( power(-1,(odd_no-1)/2) .* FACT(odd_no+1+1) ./ power(2,odd_no) ./ FACT((odd_no-1)/2+1) ./ FACT((odd_no+1)/2+1) ) * D_odd;
    X2 = power(-1,(odd_no-1)/2) .* FACT(odd_no+1+1) ./ power(2,odd_no) ./ FACT((odd_no-1)/2+1) ./ FACT((odd_no+1)/2+1);
    multifact = X2 ./ X1;
    D_odd = D_odd .* multifact;

    %=========================================================================
    %===========EVALUATION of LAMBDA =========================================
    lambda = zeros(1,N);
    for n=0:2:N-1
        lambda(n+1) = 2*c/pi * ( 2^n * D_even(1,n/2+1) * (FACT(n/2+1))^2 / FACT(n+1) )^2;
    end
    for n=1:2:N-1
        lambda(n+1) = 2*c/pi * ( 2^n * D_odd(1,(n+1)/2) * c * FACT((n-1)/2+1) * FACT((n+1)/2+1) / (3*FACT(n+1+1)) )^2;
    end

    Nn_even = (2 ./ (2*even_no + 1)) * (D_even.^2);
    Nn_odd = (2 ./ (2*odd_no + 1)) * (D_odd.^2);

    Kappan_even = FACT(even_no+1) ./ ( 2.^even_no .* FACT(even_no/2+1).^2 .* D_even(1,:));
    Kappan_odd  = 3*FACT(odd_no+2) ./ ( c * 2.^odd_no .* FACT((odd_no-1)/2+1) .* FACT((odd_no+1)/2+1) .* D_odd(1,:));

    %========================================================================
    %========= psi calculation ==============================================
    
    %====== spherical bessel functions =======
    spJ = zeros(N,size(t,2));
    for n=0:N-1
        spJ(n+1,:) = sign(t) .* sphbes(n,c*t/T);
    end
    spJ(:,floor(tt*Fs)+1) = spJ(:,floor(tt*Fs)+2);

    %====== psi in interval outside |T| =======
    psi_even = diag( sqrt( (lambda(even_no+1)./Nn_even)/T) .* Kappan_even .* ((-1).^(even_no/2))) * ( diag((-1).^(even_no/2)) .* D_even )' * spJ(even_no+1,:);
    psi_odd  = diag( sqrt( (lambda(odd_no+1)./Nn_odd)/T) .* Kappan_odd .* ((-1).^((odd_no-1)/2))) * ( diag((-1).^((odd_no-1)/2)) .* D_odd )' * spJ(odd_no+1,:);

    psi(even_no+1,:) = real(psi_even);
    psi(odd_no+1,:) = real(psi_odd);
end

function js = sphbes(nu, x)
% returns the spherical Bessel functions jnu(x)
% x is a vector or it may be a matrix if nu is a scalar
% if nu is a row and x a column vector, the output js is a matrix

[nnu lnu] = size(nu);
[nx lx] = size(x);
xm = repmat(x, 1, lnu);
js = sqrt(pi ./(2* xm)) .* (besselj(nu + 0.5, x));

end