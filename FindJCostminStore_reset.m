function [index_opt, Jmin, Jstore, x_reset]=FindJCostminStore_reset(cnt, indexes, numSalti, numStati_a, numStati_c,numStati_p,numStati_r, numeroCampioni, timeMesh, Mc1, Mc2, Aav, x0, Jnow, index_opt, JstoreN, x_reset, Tv, Pv, numStati)
persistent globCnt;
if(isempty(globCnt))
    globCnt=1;
end
Jmin=Jnow;
Jstore=JstoreN;
if(cnt>numSalti)
    %calcola Jmin
    invalid=0;
    
    % the vector is the reset states plus the lagrange multiplicators
    index_lambda=(numSalti+1)*(numStati_c);
    MO=zeros(index_lambda, index_lambda);
    BO=zeros(index_lambda, 1);
    
    
    % x0* computed. I have the initial dynamic
    
    for j=1:numSalti
        if j==1
            ev(:,1:numStati_a)=expm(Aav(1:numStati_a,1:numStati_a)*timeMesh(1,indexes(1)));    
        else
            if(timeMesh(j, indexes(j))-timeMesh(j-1, indexes(j-1))<0)
                invalid=1;
                break;
            end
            ev(:,(j-1)*numStati_a+1:j*numStati_a)=expm(Aav(1:numStati_a,numStati_a*(j)+1:(j+1)*numStati_a)*(timeMesh(j+1,indexes(j+1))-timeMesh(j,indexes(j))));  
        end
        
    end
    
    
    if(invalid)
        J=Inf;
    else
        for j=0:numSalti
            MOtemp=MO;
            BOtemp=BO;
            [MO,BO]=SolvePiece(ev,  Pv, Tv, j, numSalti, x0, MOtemp,BOtemp, numStati, numStati_a, numStati_c, numStati_p, numStati_r,1);
        end

%         MO
%         BO
%         
%         T1=Tv(:,1:numStati_a);
%         T2=Tv(:,numStati_a+1:2*numStati_a);
%         P1=Pv(:,1:numStati);
%         P2=Pv(:,numStati+1:2*numStati);
% 
%         newP1=T1'*P1*T1-ev'*T1'*P1*T1*ev;
%         newP2=T2'*P2*T2;
%         ec=ev(numStati_c+1:numStati_c+numStati_p, 1:numStati_c);
%         ep=ev(numStati_c+1:numStati_c+numStati_p, numStati_c+1:numStati_c+numStati_p);
%         er=ev(numStati_c+1:numStati_c+numStati_p, numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
%         
%         x0p=x0(numStati_c+1:numStati_c+numStati_p,1);        
%         x0r=x0(numStati_c+numStati_p+1:end,1);
%         
%         x1=ev*x0;
%         x1p=x1(numStati_c+1:numStati_c+numStati_p,1);        
%         x1r=x1(numStati_c+numStati_p+1:end,1);
%         
%         
%         
%         P1cc=newP1(1:numStati_c,1:numStati_c);        
%         P1cp=newP1(1:numStati_c,numStati_c+1:numStati_c+numStati_p);
%         P1cr=newP1(1:numStati_c,numStati_c+numStati_p+1:end);        
%         P1pp=newP1(numStati_c+1:numStati_c+numStati_p,numStati_c+1:numStati_c+numStati_p);
%         P1pr=newP1(numStati_c+1:numStati_c+numStati_p,numStati_c+numStati_p+1:end);
%         P2cc=newP2(1:numStati_c,1:numStati_c);
%         P2cp=newP2(1:numStati_c,numStati_c+1:numStati_c+numStati_p);
%         P2cr=newP2(1:numStati_c,numStati_c+numStati_p+1:end);        
%         P2pp=newP2(numStati_c+1:numStati_c+numStati_p,numStati_c+1:numStati_c+numStati_p);
%         P2pr=newP2(numStati_c+1:numStati_c+numStati_p,numStati_c+numStati_p+1:end);
%         MO=[P1cc+ec'*P2pp*ec, ec'*P2cp'; P2cp*ec, P2cc]
%         BO=[-P1cp*x0p-P1cr*x0r-ec'*P2pr*x1r-ec'*P2pp*ep*x0p-ec'*P2pp*er*x0r;
%             -P2cr*x1r-P2cp*ep*x0p-P2cp*er*x0r]
        
        
        %det(MO)
        if(det(MO)==0)
            J=Inf;
        else
            resets=(MO^-1)*BO;
            %resets
            
            xo=zeros(numStati_a*(numSalti+1),1);
            xx=zeros(numStati_a*numSalti,1);
            xo(1:numStati_a,1)=[resets(1:numStati_c,1); x0(numStati_c+1:end,1)];
            dynamics=expm(Aav(1:numStati_a,1:numStati_a)*timeMesh(1,indexes(1)));
            
            for j=1:numSalti
                xi=dynamics*xo((j-1)*numStati_a+1:j*numStati_a,1);
                xx(numStati_a*(j-1)+1:j*numStati_a,1)=xi;
                x_res=resets(j*numStati_c+1:(j+1)*numStati_c,1);
                xo(numStati_a*(j)+1:numStati_a*(j+1),1)=[x_res(1:numStati_c,1); xi(numStati_c+1:end,1)];
                if((j+2)<=(numSalti+1))
                    dynamics =expm(Aav(1:numStati_a,numStati_a*(j)+1:(j+1)*numStati_a)*(timeMesh(j+1,indexes(j+1))-timeMesh(j,indexes(j))));
                end
            end
            % xx
            % xo
            J=xx'*Mc1*xx+xo'*Mc2*xo;
            
        end
    end
    if(J<Jmin)
        Jmin=J;
        index_opt=indexes;
        x_reset=xo;
    end
    Jstore(globCnt)=J;
    globCnt=globCnt+1;
else
    for i=1:numeroCampioni
        indexes(cnt)=i;
        cnt=cnt+1;
        [index_opt, Jmin, Jstore, x_reset]=FindJCostminStore_reset(cnt, indexes, numSalti, numStati_a, numStati_c,numStati_p,numStati_r, numeroCampioni, timeMesh, Mc1, Mc2, Aav, x0, Jmin, index_opt, Jstore, x_reset, Tv, Pv, numStati);
        cnt=cnt-1;
    end
end
end

% function [index_opt, Jmin, Jstore, x_reset]=FindJCostminStore_reset(cnt, indexes, numSalti, numStati_a, numStati_c,numeroCampioni, timeMesh, Mc1, Mc2, Aav, x0, Jnow, index_opt, JstoreN, x_reset, Tv, Pv, numStati)
% persistent globCnt;
% if(isempty(globCnt))
%     globCnt=1;
% end
% Jmin=Jnow;
% Jstore=JstoreN;
% if(cnt>numSalti)
%     %calcola Jmin
%     invalid=0;
%     
%     dynamics=expm(Aav(1:numStati_a,1:numStati_a)*timeMesh(1,indexes(1)));
%     P1=Pv(1:numStati,1:numStati);
%     T=Tv(1:numStati,1:numStati_a);
%     newP=T'*P1*T-dynamics'*T'*P1*T*dynamics;
%     Rc=newP(1:numStati_c,1:numStati_c);
%     Rpr=newP(1:numStati_c,numStati_c+1:end);
%     xnc=x0(numStati_c+1:end,1);
%     xc_o=-(Rc^-1)*Rpr*xnc;
%     x_o=[xc_o; xnc];
%     xo(1:numStati_a,1)=x_o;
%     % x0* computed. I have the initial dynamic
%     
%     for j=1:numSalti
%         %evolve the dynamic computed beginning for the optimal
%         xi=dynamics*x_o;
%         % computed P2
%         Pnext=Pv(1:numStati,numStati*(j)+1:numStati*(j+1));
%         % computed T1
%         T=Tv(1:numStati,numStati_a*(j-1)+1:numStati_a*j);
%         % computed T2
%         Tnext=Tv(1:numStati,numStati_a*(j)+1:numStati_a*(j+1));
%         % P2-eAt*P3*eAt
%         % computes new dynamic if needed
%         newP=Tnext'*Pnext*Tnext;
%         if((j+2)<=(numSalti+1))
%             if(timeMesh(j, indexes(j))-timeMesh(j-1, indexes(j-1))<0)
%                 invalid=1;
%                 break;
%             end
%             dynamics =expm(Aav(1:numStati_a,numStati_a*(j)+1:(j+1)*numStati_a)*(timeMesh(j+1,indexes(j+1))-timeMesh(j,indexes(j))));
%             newP=newP-dynamics'*Tnext'*Pnext*Tnext*dynamics;
%         end
%         %compute new reset
%         Rc=newP(1:numStati_c,1:numStati_c);
%         Rpr=newP(1:numStati_c,numStati_c+1:end);
%         xnc=xi(numStati_c+1:end,1);
%         xc_o=-(Rc^-1)*Rpr*xnc;
%         x_o=[xc_o; xnc];
%         xx(numStati_a*(j-1)+1:j*numStati_a,1)=xi;
%         xo(numStati_a*(j)+1:numStati_a*(j+1),1)=x_o;
%     end
%     if(invalid)
%         J=Inf;
%     else
%         xx
%         xo
%         J=xx'*Mc1*xx+xo'*Mc2*xo
%     end
%     if(J<Jmin)
%         Jmin=J;
%         index_opt=indexes;
%         x_reset=xo;
%     end
%     
%     Jstore(globCnt)=J;
%     globCnt=globCnt+1;
% else
%     for i=1:numeroCampioni
%         indexes(cnt)=i;
%         cnt=cnt+1;
%         [index_opt, Jmin, Jstore, x_reset]=FindJCostminStore_reset(cnt, indexes, numSalti, numStati_a, numStati_c,numeroCampioni, timeMesh, Mc1, Mc2, Aav, x0, Jmin, index_opt, Jstore, x_reset, Tv, Pv, numStati);
%         cnt=cnt-1;
%     end
% end
% end