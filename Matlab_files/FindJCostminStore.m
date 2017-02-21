function [index_opt, Jmin, Jstore]=FindJCostminStore(cnt, indexes, numSalti, numStati,numeroCampioni, timeMesh, Mc, Aev, x0, Jnow, index_opt, JstoreN)
persistent globCnt;
if(isempty(globCnt))
    globCnt=1;
end

Jmin=Jnow;
Jstore=JstoreN;
if(cnt>numSalti)
    %calcola Jmin
    invalid=0;
 %   for i=1:numSalti
        xi=x0;
        for j=1:numSalti
            if(invalid)
                break;
            end
            if(j==1)
                xi=expm(Aev(1:numStati,numStati*(j-1)+1:j*numStati)*timeMesh(j,indexes(j)))*xi;
            else
                if(timeMesh(j, indexes(j))-timeMesh(j-1, indexes(j-1))>=0)
                    xi=expm(Aev(1:numStati,numStati*(j-1)+1:j*numStati)*(timeMesh(j,indexes(j))-timeMesh(j-1,indexes(j-1))))*xi;
                else 
                    invalid=1;
                    break;
                end
            end
            xx(numStati*(j-1)+1:j*numStati,1)=xi;
        end
  %  end
    if(invalid)
        J=Inf;
    else
        J=xx'*Mc*xx;
    end
    if(J<Jmin)
        Jmin=J;
        index_opt=indexes;
    end
  
    Jstore(globCnt)=J;
    globCnt=globCnt+1;
else
    for i=1:numeroCampioni
        indexes(cnt)=i;
        cnt=cnt+1;
        [index_opt, Jmin, Jstore]=FindJCostminStore(cnt, indexes, numSalti, numStati, numeroCampioni, timeMesh, Mc, Aev, x0, Jmin, index_opt,Jstore);
        cnt=cnt-1;
    end
end
end