function [index_opt, Jmin, x_reset]=FindJCostmin_reset(cnt, indexes, numSalti, numStati_a, numStati_c,numeroCampioni, timeMesh, Mc1, Mc2, Mc3, Aav, x0, Jnow, index_opt, x_reset)
Jmin=Jnow;

if(cnt>numSalti)
    %calcola Jmin
    invalid=0;
    
    dynamics=expm(Aav(1:numStati_a,1:numStati_a)*timeMesh(1,indexes(1)));
    P1=Pv(1:numStati,1:numStati);
    T=Tv(1:numStati,1:numStati_a);
    newP=T'*P1*T-dynamics'*T'*P1*T*dynamics;
    Rc=newP(1:numStati_c,1:numStati_c);
    Rpr=newP(1:numStati_c,numStati_c+1:end);
    xnc=x0(numStati_c+1:end,1);
    xc_o=-(Rc^-1)*Rpr*xnc;
    x_o=[xc_o; xnc];
    xo(1:numStati_a,1)=x_o;
    % x0* computed. I have the initial dynamic
    for i=1:numSalti
        if(invalid)
            break;
        end
        for j=1:i
            %evolve the dynamic computed beginning for the optimal
            xi=dynamics*x_o;
            % computed P2
            Pnext=Pv(1:numStati,numStati*(j)+1:numStati*(j+1));
            % computed T1
            T=Tv(1:numStati,numStati_a*(j-1)+1:numStati_a*j);
            % computed T2
            Tnext=Tv(1:numStati,numStati_a*(j)+1:numStati_a*(j+1));
            % P2-eAt*P3*eAt
            % computes new dynamic if needed
            newP=T'*Pnext*T;
            if((j+2)<=(numSalti+1))
                if(timeMesh(j, indexes(j))-timeMesh(j-1, indexes(j-1))<0)
                    invalid=1;
                    break;
                end
                dynamics =expm(Aav(1:numStati_a,numStati_a*(j)+1:(j+1)*numStati_a)*(timeMesh(j+1,indexes(j+1))-timeMesh(j,indexes(j))));
                newP=newP-dynamics'*Tnext'*Pnext*Tnext*dynamics;
            end
            %compute new reset
            Rc=newP(1:numStati_c,1:numStati_c);
            Rpr=newP(1:numStati_c,numStati_c+1:end);
            xnc=xi(numStati_c+1:end,1);
            xc_o=-(Rc^-1)*Rpr*xnc;
            x_o=[xc_o; xnc];
        end
        xx(numStati_a*(i-1)+1:i*numStati_a,1)=xi;
        xo(numStati_a*(i)+1:i*numStati_a*(i+1),1)=x_o;
    end
    if(invalid)
        J=Inf;
    else
        J=xx'*Mc1*xx+xo'*Mc2*xo;
    end
    if(J<Jmin)
        Jmin=J;
        index_opt=indexes;
        x_reset=xo;
    end
else
    for i=1:numeroCampioni
        indexes(cnt)=i;
        cnt=cnt+1;
        [index_opt, Jmin, x_reset]=FindJCostmin_reset(cnt, indexes, numSalti, numStati_a, numStati_c,numeroCampioni, timeMesh, Mc1, Mc2, Mc3, Aav, x0, Jmin, index_opt, x_reset);
        cnt=cnt-1;
    end
end
end