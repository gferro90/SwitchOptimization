function[MO, BO]= SolvePieceFFBL(ev,  Pv, Tv, j, numSalti, x0, MOtemp, BOtemp, numStati, numStati_a, numStati_c, numStati_p, numStati_r, numStati_w, Em, flag)
  MO=MOtemp;
  BO=BOtemp;
  numStati_cw=numStati_c+numStati_w;

  if(j==0)
      E1=ev(:,1:numStati_a);
      P1=Pv(1:numStati,1:numStati);
      T=Tv(1:numStati,1:numStati_a);
      if(flag)
        newP=T'*P1*T-E1'*T'*P1*T*E1;
      else
        newP=Pv;
      end
      Rcc=newP(1:numStati_c,1:numStati_c);
      Rcpr=newP(1:numStati_c,numStati_c+1:end-numStati_w);
      Rcw=newP(1:numStati_c,end-numStati_w+1:end);
      Rwc=newP(end-numStati_w+1:end,1:numStati_c);
      Rwpr=newP(end-numStati_w+1,numStati_c+1:end-numStati_w);
      Rww=newP(end-numStati_w+1,end-numStati_w+1);
      Rcw=[Rcc, Rcw; Rwc, Rww];
      Rpr=[Rcpr; Rwpr];
      xnc=x0(numStati_c+1:end-numStati_w,1);
    
      
      MO(1:numStati_cw, 1:numStati_cw)=MO(1:numStati_cw, 1:numStati_cw)+Rcw;
    
      BO(1:numStati_cw,1)=BO(1:numStati_cw,1)-Rpr*xnc;
  else
        %extract the last e^At
        Ej=ev(:,(j-1)*numStati_a+1:j*numStati_a);
        % computed P2
        Pnext=Pv(1:numStati,numStati*(j)+1:numStati*(j+1));
        
        %extract Ec, Ep and Er pieces from Ej
        Ec=Ej(numStati_c+1:numStati_c+numStati_p, 1:numStati_c);
        Ep=Ej(numStati_c+1:numStati_c+numStati_p, numStati_c+1:numStati_c+numStati_p);
        Er=Ej(numStati_c+1:numStati_c+numStati_p, numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
        Ew=Ej(numStati_c+1:numStati_c+numStati_p, end-numStati_w+1:end);

        % computed T2
        Tnext=Tv(1:numStati,numStati_a*(j)+1:numStati_a*(j+1));        

        % compute the piece 2xcq P(q+q, cpp)epq xp(q-1)
        index_lambda=(numSalti+1)*(numStati_cw);

        row_M=zeros(numStati_p, index_lambda);
        xnc=x0(numStati_c+1:end-numStati_w,1);
        chain_xp=xnc(1:numStati_p,1);
        xr_prev=xnc(numStati_p+1:end,1);
        
        if(flag)
            T=Tv(1:numStati,numStati_a*(j-1)+1:numStati_a*j); 
            Lm=Em*Tnext;
            Lm1=Em*T*Ej;
            chainLm=Lm;
            chainLb=x0(numStati_c+1:numStati_c+numStati_p,1);
            chainLm_1=Lm1;
            chainLb_1=x0(numStati_c+1:numStati_c+numStati_p,1);
            xr_pp=0;
        end
        for k=1:j
            % compute the piece 2xcq P(q+1, cpp)epq xp(q-1)
            Ek=ev(:,(k-1)*numStati_a+1:(k)*numStati_a);
            
            %extract Ec, Ep and Er pieces from Ek
            Eck=Ek(numStati_c+1:numStati_c+numStati_p, 1:numStati_c);
            Epk=Ek(numStati_c+1:numStati_c+numStati_p, numStati_c+1:numStati_c+numStati_p);
            Erk=Ek(numStati_c+1:numStati_c+numStati_p, numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
            Ewk=Ek(numStati_c+1:numStati_c+numStati_p, end-numStati_w+1:end);
            
            Ecwk=[Eck, Ewk];
            
            row_M=Epk*row_M;
            row_M(:,(k-1)*numStati_cw+1: k*numStati_cw)=Ecwk;
            chain_xp=Ep*chain_xp+Er*xr_prev;
            if(flag)
                chainLb=(Lm(:,numStati_c+1:numStati_c+numStati_p))*Epk*chainLb+(Lm(:,numStati_c+1:numStati_c+numStati_p))*Erk*xr_prev;
                if(j>1)
                    Ek_1=ev(:,(k-2)*numStati_a+1:(k-1)*numStati_a);
                    
                    %extract Ec, Ep and Er pieces from Ek
                    Epk_1=Ek_1(numStati_c+1:numStati_c+numStati_p, numStati_c+1:numStati_c+numStati_p);
                    Erk_1=Ek_1(numStati_c+1:numStati_c+numStati_p, numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
                    
                    chainLb_1=(Lm1(:,numStati_c+1:numStati_c+numStati_p))*Epk_1*chainLb_1+(Lm1(:,numStati_c+1:numStati_c+numStati_p))*Erk_1*xr_pp;
                else
                    chainLb_1=(Lm1(:,numStati_c+1:numStati_c+numStati_p))*chainLb_1;
                end
                
            end
            xr_pp=xr_prev;
            xr_prev=Ek(numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r, numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r)*xr_prev;
        end
        if(flag)
            BO(index_lambda+j,1)=BO(index_lambda+j,1)-chainLb-(Lm(:,numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r))*xr_prev;
            BO(index_lambda+j,1)=BO(index_lambda+j,1)+chainLb_1+(Lm1(:,numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r))*xr_pp;
        end
        % P2-eAt*P3*eAt
        % computes new dynamic if needed
        if(flag)
            newP=Tnext'*Pnext*Tnext;
            if((j+2)<=(numSalti+1))
                dynamics =ev(:,(k)*numStati_a+1:(k+1)*numStati_a);
                newP=newP-dynamics'*Tnext'*Pnext*Tnext*dynamics;
            end
        else
            newP=Pv;
        end
        
        
        chainP=newP;

        for k=j:-1:1
            %compute the piect 2xp(q-1)epq P(q+1, prqp) xrq
            
            % compute the piece 2xcq P(q+1, cpp)epq xp(q-1)
            Ek=ev(:,(k-1)*numStati_a+1:(k)*numStati_a);

            %extract Ec, Ep and Er pieces from Ek
            Eck=Ek(numStati_c+1:numStati_c+numStati_p, 1:numStati_c);
            Epk=Ek(numStati_c+1:numStati_c+numStati_p, numStati_c+1:numStati_c+numStati_p);
            Erk=Ek(numStati_c+1:numStati_c+numStati_p, numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
            Ewk=Ek(numStati_c+1:numStati_c+numStati_p, end-numStati_w+1:end);
            
            Ecwk=[Eck, Ewk];
            BO(numStati_cw*(k-1)+1:numStati_cw*k,1)=BO(numStati_cw*(k-1)+1: numStati_cw*k,1)-Ecwk'*chainP(numStati_c+1:numStati_c+numStati_p,numStati_c+numStati_p+1:end-numStati_w)*xr_prev;

            chainP=Epk'*chainP(numStati_c+1:numStati_c+numStati_p,numStati_c+numStati_p+1:end-numStati_w);
            if(flag)
                % lambda chains
                MO(index_lambda+j,numStati_cw*(k-1)+1:numStati_cw*k)=MO(index_lambda+j,numStati_cw*(k-1)+1:numStati_cw*k)+chainLm(:,numStati_c+1:numStati_c+numStati_p)*Ecwk;
                chainLm=chainLm(:,numStati_c+1:numStati_c+numStati_p)*Epk;
                if(j>1)
                    Ek_1=ev(:,(k-2)*numStati_a+1:(k-1)*numStati_a);
                    Eck_1=Ek_1(numStati_c+1:numStati_c+numStati_p, 1:numStati_c);
                    Epk_1=Ek_1(numStati_c+1:numStati_c+numStati_p, numStati_c+1:numStati_c+numStati_p);
                    Ewk_1=Ek_1(numStati_c+1:numStati_c+numStati_p, end-numStati_w+1:end);

                    Ecwk_1=[Eck_1, Ewk_1];
                    MO(index_lambda+j,numStati_cw*(k-2)+1:numStati_cw*(k-1))=MO(index_lambda+j,numStati_cw*(k-2)+1:numStati_cw*(k-1))-chainLm_1(:,numStati_c+1:numStati_c+numStati_p)*Ecwk_1;
                    chainLm_1=chainLm_1(:,numStati_c+1:numStati_c+numStati_p)*Epk_1;
                end
                
                MO(index_lambda+j, numStati_cw*k+1:numStati_cw*(k+1))=MO(index_lambda+j, numStati_cw*k+1:numStati_cw*(k+1))+[Lm(1,1:numStati_c), Lm(1,end-numStati_w+1:end)];
                MO(index_lambda+j, numStati_cw*(k-1)+1:numStati_cw*(k))=MO(index_lambda+j, numStati_cw*(k-1)+1:numStati_cw*(k))-[Lm1(1,1:numStati_c), Lm1(1,end-numStati_w+1:end)];
                
                
                MO(numStati_cw*(k-1)+1:numStati_cw*k,index_lambda+j)=MO(index_lambda+j,numStati_cw*(k-1)+1:numStati_cw*k)';
                MO(numStati_cw*k+1:numStati_cw*(k+1),index_lambda+j)=MO(index_lambda+j, numStati_cw*k+1:numStati_cw*(k+1));

            end
        end
        

        
        %build this piece
        Rcc=newP(1:numStati_c,1:numStati_c);
        Rcpr=newP(1:numStati_c,numStati_c+1:end-numStati_w);
        Rcw=newP(1:numStati_c,end-numStati_w+1:end);
        Rwc=newP(end-numStati_w+1:end,1:numStati_c);
        Rwpr=newP(end-numStati_w+1,numStati_c+1:end-numStati_w);
        Rww=newP(end-numStati_w+1,end-numStati_w+1);
        Rcw=[Rcc, Rcw; Rwc, Rww];
        Rpr=[Rcpr; Rwpr];

        MO(j*numStati_cw+1:(j+1)*numStati_cw,1:index_lambda)=MO(j*numStati_cw+1:(j+1)*numStati_cw,1:index_lambda)+Rpr(:,1:numStati_p)*row_M;
        MO(1:index_lambda,j*numStati_cw+1:(j+1)*numStati_cw)=MO(1:index_lambda,j*numStati_cw+1:(j+1)*numStati_cw)+ MO(j*numStati_cw+1:(j+1)*numStati_cw,1:index_lambda)';
        
        MO(j*numStati_cw+1:(j+1)*numStati_cw,j*numStati_cw+1:(j+1)*numStati_cw)=MO(j*numStati_cw+1:(j+1)*numStati_cw,j*numStati_cw+1:(j+1)*numStati_cw)+Rcw;
        
        
        BO(j*numStati_cw+1:(j+1)*numStati_cw,1)=BO(j*numStati_cw+1:(j+1)*numStati_cw,1)-Rpr(:,1:numStati_p)*chain_xp-Rpr(:,numStati_p+1:end)*xr_prev;
    
        %build the matrix for the recursive
        Px=[Ec, Ep, Er, Ew];
        Prec=Px'*newP(numStati_c+1:numStati_c+numStati_p,numStati_c+1:numStati_c+numStati_p)*Px;
        %go recursively on this
        [MO, BO]=SolvePieceFFBL(ev,  Prec, Tv, j-1, numSalti, x0, MO, BO, numStati, numStati_a, numStati_c, numStati_p, numStati_r,numStati_w,Em,0);
    
  end
end