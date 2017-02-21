function [index_opt, Jmin, Jstore, x_reset]=FindJCostminStore_resetFFBumpless(cnt, indexes, numSalti, numStati_a, numStati_c,numStati_p,numStati_r, numStati_w, numeroCampioni, timeMesh, Mc1, Mc2, Aav, x0, Jnow, index_opt, JstoreN, x_reset, Tv, Pv, numStati,E)
persistent globCnt;
if(isempty(globCnt))
    globCnt=1;
end
Jmin=Jnow;
Jstore=JstoreN;
if(cnt>numSalti)
    %calcola Jmin
    invalid=0;
    numStati_cw=numStati_c+numStati_w;

    % the vector is the reset states plus the lagrange multiplicators
    index_lambda=(numSalti+1)*(numStati_cw);
    MO=zeros(index_lambda+numSalti, index_lambda+numSalti);
    BO=zeros(index_lambda+numSalti, 1);
    
    
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
            [MO,BO]=SolvePieceFFBL(ev,  Pv, Tv, j, numSalti, x0, MOtemp,BOtemp, numStati, numStati_a, numStati_c, numStati_p, numStati_r, numStati_w,E,1);
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
%         ec=[ev(numStati_c+1:numStati_c+numStati_p, 1:numStati_c), ev(numStati_c+1:numStati_c+numStati_p, end-numStati_w+1:end)];
%         ep=ev(numStati_c+1:numStati_c+numStati_p, numStati_c+1:numStati_c+numStati_p);
%         er=ev(numStati_c+1:numStati_c+numStati_p, numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
%         
%         x0p=x0(numStati_c+1:numStati_c+numStati_p,1);        
%         x0r=x0(numStati_c+numStati_p+1:end-numStati_w,1);
%         
%         x1=ev*x0;
%         x1p=x1(numStati_c+1:numStati_c+numStati_p,1);        
%         x1r=x1(numStati_c+numStati_p+1:end-numStati_w,1);
%         
%         
%         
%         P1cc=[newP1(1:numStati_c,1:numStati_c),         newP1(1:numStati_c,end-numStati_w+1:end);
%               newP1(end-numStati_w+1:end,1:numStati_c), newP1(end-numStati_w+1:end,end-numStati_w+1:end)];        
%         P1cp=[newP1(1:numStati_c,numStati_c+1:numStati_c+numStati_p);
%               newP1(end-numStati_w+1:end,numStati_c+1:numStati_c+numStati_p)];
%         P1cr=[newP1(1:numStati_c,numStati_c+numStati_p+1:end-numStati_w);
%               newP1(end-numStati_w+1:end,numStati_c+numStati_p+1:end-numStati_w)];
%         P1pp=newP1(numStati_c+1:numStati_c+numStati_p,numStati_c+1:numStati_c+numStati_p);
%         P1pr=newP1(numStati_c+1:numStati_c+numStati_p,numStati_c+numStati_p+1:end-numStati_w);
%         
%         P2cc=[newP2(1:numStati_c,1:numStati_c),         newP2(1:numStati_c,end-numStati_w+1:end);
%               newP2(end-numStati_w+1:end,1:numStati_c), newP2(end-numStati_w+1:end,end-numStati_w+1:end)];        
%         P2cp=[newP2(1:numStati_c,numStati_c+1:numStati_c+numStati_p);
%               newP2(end-numStati_w+1:end,numStati_c+1:numStati_c+numStati_p)];
%         P2cr=[newP2(1:numStati_c,numStati_c+numStati_p+1:end-numStati_w);
%               newP2(end-numStati_w+1:end,numStati_c+numStati_p+1:end-numStati_w)];
%         P2pp=newP2(numStati_c+1:numStati_c+numStati_p,numStati_c+1:numStati_c+numStati_p);
%         P2pr=newP2(numStati_c+1:numStati_c+numStati_p,numStati_c+numStati_p+1:end-numStati_w);
%         
%         E1=E*T1*ev;
%         E2=E*T2;
%         
%         E1c=[E1(:,1:numStati_c), E1(:, end-numStati_w+1:end)];
%         E1p=E1(:,numStati_c+1:numStati_c+numStati_p);
%         E1r=E1(:,numStati_c+numStati_p+1:end-numStati_w);
%         
%         E2c=[E2(:,1:numStati_c), E2(:, end-numStati_w+1:end)];
%         E2p=E2(:,numStati_c+1:numStati_c+numStati_p);
%         E2r=E2(:,numStati_c+numStati_p+1:end-numStati_w);
% 
%     
% 
%         MO=[P1cc+ec'*P2pp*ec, ec'*P2cp', (E2p*ec-E1c)'; 
%             P2cp*ec,          P2cc,       E2c';
%             E2p*ec-E1c        E2c         0]  
%         
% 
%         BO=[-P1cp*x0p-P1cr*x0r-ec'*P2pr*x1r-ec'*P2pp*ep*x0p-ec'*P2pp*er*x0r;
%             -P2cr*x1r-P2cp*ep*x0p-P2cp*er*x0r;
%             -E2p*ep*x0p-E2p*er*x0r-E2r*x1r+E1p*x0p+E1r*x0r]
%         
%         
        %det(MO)
        if(det(MO)==0)
            J=Inf;
        else
            resets=(MO^-1)*BO;
            %resets
            
            xo=zeros(numStati_a*(numSalti+1),1);
            xx=zeros(numStati_a*numSalti,1);
            xo(1:numStati_a,1)=[resets(1:numStati_c,1); x0(numStati_c+1:end-numStati_w,1); resets(numStati_c+1:numStati_c+numStati_w,1)];
            dynamics=expm(Aav(1:numStati_a,1:numStati_a)*timeMesh(1,indexes(1)));
            
            for j=1:numSalti
                xi=dynamics*xo((j-1)*numStati_a+1:j*numStati_a,1);
                xx(numStati_a*(j-1)+1:j*numStati_a,1)=xi;
                x_res=resets(j*numStati_cw+1:(j+1)*numStati_cw,1);
                xo(numStati_a*(j)+1:numStati_a*(j+1),1)=[x_res(1:numStati_c,1); xi(numStati_c+1:end-numStati_w,1);x_res(numStati_c+1:numStati_c+numStati_w,1)];
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
        [index_opt, Jmin, Jstore, x_reset]=FindJCostminStore_resetFFBumpless(cnt, indexes, numSalti, numStati_a, numStati_c,numStati_p,numStati_r, numStati_w, numeroCampioni, timeMesh, Mc1, Mc2, Aav, x0, Jmin, index_opt, Jstore, x_reset, Tv, Pv, numStati,E);
        cnt=cnt-1;
    end
end
end


% 
% function [index_opt, Jmin, Jstore, x_reset]=FindJCostminStore_resetFFBumpless(cnt, indexes, numSalti, numStati_a, numStati_c, numStati_p, numStati_r, numStati_w, numeroCampioni, timeMesh, Mc1, Mc2, Aav, x0, Jnow, index_opt, JstoreN, x_reset, Tv, Pv, numStati, C)
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
%     
%     numStati_cw=numStati_c+numStati_w;
%     dynamics=expm(Aav(1:numStati_a,1:numStati_a)*timeMesh(1,indexes(1)));
%     P1=Pv(1:numStati,1:numStati);
%     T=Tv(1:numStati,1:numStati_a);
%     Tcw=[T(:,1:numStati_c), T(:,end-numStati_w+1:end)];
%     newP=T'*P1*T-dynamics'*T'*P1*T*dynamics;
%     Rc_c=newP(1:numStati_c,1:numStati_c);
%     Rc_w=newP(1:numStati_c,end-numStati_w+1:end);
%     Rc_pr=newP(1:numStati_c,numStati_c+1:end-numStati_w);
%     Rw_w=newP(end-numStati_w+1:end,end-numStati_w+1:end);
%     Rw_c=newP(end-numStati_w+1:end,1:numStati_c);
%     Rw_pr=newP(end-numStati_w+1:end,numStati_c+1:end-numStati_w);
%     Rc=[Rc_c,    Rc_w;  Rw_c,    Rw_w];
%     Rpr=[Rc_pr;  Rw_pr];
%     xnc=x0(numStati_c+1:end-numStati_w,1);
%     
%     index_lambda=(numSalti+1)*(numStati_cw);
%     % the vector is the reset states plus the lagrange multiplicators
%     MO=zeros(index_lambda+numSalti, index_lambda+numSalti);
%     BO=zeros(index_lambda+numSalti, 1);
%     dynamics_cw=[dynamics(:,1:numStati_c), dynamics(:, end-numStati_w+1, end)];
%     
%     MO(1:numStati_cw, 1:numStati_cw)=Rc;
%     
%     BO(1:numStati_cw,1)=-Rpr*xnc;
%     chain_xp=xnc(1:numStati_p,1);
%     row_M=zeros(numStati_p, index_lambda);
%     xr_prev=x0(numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r,1);
%     
%     % x0* computed. I have the initial dynamic
%     
%     for j=1:numSalti
%         
%         Ep=dynamics(numStati_c+1:numStati_c+numStati_p, numStati_c+1:numStati_c+numStati_p);
%         Er=dynamics(numStati_c+1:numStati_c+numStati_p, numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
%         
%         Ec=dynamics(numStati_c+1:numStati_c+numStati_p, 1:numStati_c);
%         Ew=dynamics(numStati_c+1:numStati_c+numStati_p,end-numStati_w+1: end);
%         %evolve the dynamic computed beginning for the optimal
%         xr=dynamics(numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r,numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r)*xr_prev;
%         
%         % computed P2
%         Pnext=Pv(1:numStati,numStati*(j)+1:numStati*(j+1));
%         % computed T1
%         T=Tv(1:numStati,numStati_a*(j-1)+1:numStati_a*(j));
%         Tcw=[T(:,1:numStati_c), T(:,end-numStati_w+1:end)];
%         
%         row_M=Ep*row_M;
%         row_M(:,(j-1)*numStati_cw+1: j*numStati_cw)=[Ec, Ew];
%         % computed T2
%         Tnext=Tv(1:numStati,numStati_a*(j)+1:numStati_a*(j+1));        
%                 
%         lambM=Tnext'*C';
% 
%         lambM=C*T*dynamics;
%         lambM1=C*Tnext;
%         MO(index_lambda+j, (j-1)*numStati_cw+1:j*numStati_cw)=-[lambM(1,1:numStati_c) lambM(1,end-numStati_w+1:end)]+lambM1(1,numStati_c+1:numStati_c+numStati_p)*[Ec Ew];
%         MO(index_lambda+j, (j)*numStati_cw+1:(j+1)*numStati_cw)=[lambM1(1,1:numStati_c) lambM1(1,end-numStati_w+1:end)];
%         
%         MO(j*numStati_cw+1:(j+1)*numStati_cw, index_lambda+j)=MO(index_lambda+j, (j)*numStati_cw+1:(j+1)*numStati_cw)';
%         MO((j-1)*numStati_cw+1:j*numStati_cw,index_lambda+j)=MO(index_lambda+j, (j-1)*numStati_cw+1:j*numStati_cw)';
%         chain_xp=Ep*chain_xp+Er*xr_prev;
% 
%         BO(index_lambda+j,1)=-lambM1(1,numStati_c+1:numStati_c+numStati_p)*chain_xp-lambM1(1,numStati_c+numStati_p+1:end-numStati_w)*xr+lambM(1,numStati_c+1:end-numStati_w)*xnc;
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
%         dynamics_cw=[dynamics(:,1:numStati_c), dynamics(:, end-numStati_w+1, end)];
% 
%         Rc_c=newP(1:numStati_c,1:numStati_c);
%         Rc_w=newP(1:numStati_c,end-numStati_w+1:end);
%         Rc_pr=newP(1:numStati_c,numStati_c+1:end-numStati_w);
%         Rw_w=newP(end-numStati_w+1:end,end-numStati_w+1:end);
%         Rw_c=newP(end-numStati_w+1:end,1:numStati_c);
%         Rw_pr=newP(end-numStati_w+1:end,numStati_c+1:end-numStati_w);
%         Rc=[Rc_c,    Rc_w;  Rw_c,    Rw_w];
%         Rpr=[Rc_pr;  Rw_pr];
%    
%         MO(j*numStati_cw+1:(j+1)*numStati_cw,1:index_lambda)=Rpr(:,1:numStati_p)*row_M;
%         MO(1:index_lambda,j*numStati_cw+1:(j+1)*numStati_cw)= MO(j*numStati_cw+1:(j+1)*numStati_cw,1:index_lambda)';
%         
%         MO(j*numStati_cw+1:(j+1)*numStati_cw,j*numStati_cw+1:(j+1)*numStati_cw)=Rc;
%         
% 
%         BO(j*numStati_cw+1:(j+1)*numStati_cw,1)=-Rpr(:,1:numStati_p)*chain_xp-Rpr(:,numStati_p+1:end)*xr;
%         
%         xr_prev=xr;
%     end
%     
%     
%     if(invalid)
%         J=Inf;
%     else
%         MO 
%         BO
%         %MO
%         %det(MO)
%         if(det(MO)==0)
%             J=Inf;
%         else
%             resets=(MO^-1)*BO;
%             %resets
%             for(kk=1:numSalti)
%                 
%                 %resets(index_lambda+1,1)
%                 if(resets(index_lambda+kk,1)==0)
%                     invalid=1;
%                     break;
%                 end
%             end
%             if(invalid)
%                 J=Inf;
%             else
%                 xo=zeros(numStati_a*(numSalti+1),1);
%                 xx=zeros(numStati_a*numSalti,1);
%                 xo(1:numStati_a,1)=[resets(1:numStati_c,1); x0(numStati_c+1:numStati_c+numStati_p+numStati_r,1); resets(numStati_c+1:numStati_c+numStati_w,1)];
%                 dynamics=expm(Aav(1:numStati_a,1:numStati_a)*timeMesh(1,indexes(1)));
%                 
%                 for j=1:numSalti
%                     xi=dynamics*xo((j-1)*numStati_a+1:j*numStati_a,1);
%                     xx(numStati_a*(j-1)+1:j*numStati_a,1)=xi;
%                     x_res=resets(j*numStati_cw+1:(j+1)*numStati_cw,1);
%                     xo(numStati_a*(j)+1:numStati_a*(j+1),1)=[x_res(1:numStati_c,1); xi(numStati_c+1:numStati_c+numStati_p+numStati_r,1); x_res(numStati_c+1:end,1)];
%                     if((j+2)<=(numSalti+1))
%                         dynamics =expm(Aav(1:numStati_a,numStati_a*(j)+1:(j+1)*numStati_a)*(timeMesh(j+1,indexes(j+1))-timeMesh(j,indexes(j))));
%                     end
%                 end
%                % xx
%                % xo
%                 J=xx'*Mc1*xx+xo'*Mc2*xo;
%             end
%         end
%     end
%     if(J<Jmin)
%         Jmin=J;
%         index_opt=indexes;
%         x_reset=xo;
%     end
%     Jstore(globCnt)=J;
%     globCnt=globCnt+1;
% else
%     for i=1:numeroCampioni
%         indexes(cnt)=i;
%         cnt=cnt+1;
%         [index_opt, Jmin, Jstore, x_reset]=FindJCostminStore_resetFFBumpless(cnt, indexes, numSalti, numStati_a, numStati_c,numStati_p, numStati_r, numStati_w, numeroCampioni, timeMesh, Mc1, Mc2, Aav, x0, Jmin, index_opt, Jstore, x_reset, Tv, Pv, numStati, C);
%         cnt=cnt-1;
%     end
% end
% end